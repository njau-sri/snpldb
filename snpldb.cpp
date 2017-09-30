#include <cmath>
#include <limits>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "vcfio.h"
#include "cmdline.h"
#include "split.h"
#include "util.h"

using std::size_t;
using std::log;

namespace {

    struct Par
    {
        std::string vcf;
        std::string blk;
        std::string out;
        double maf = 0.01;
        double inform = 0.95;
        int maxlen = 200000;
        int llim = 70;
        int ulim = 98;
        int recomb = 90;
        int batch = 20000;
    };

    struct Block
    {
        int first;
        int last;
        int length;
        float inform;
    };

    Par par;

    // Gabriel, S.B. et al. The structure of haplotype blocks in the human genome. Science, 2002, 296(5576): 2225-9.
    // Barrett, J.C. et al. Haploview: analysis and visualization of LD and haplotype maps. Bioinformatics, 2005, 21: 263-5.

    // Estimates haplotype frequencies via the EM algorithm
    // AB, Ab, aB, ab, AaBb
    void calc_hap_prob_EM(int n11, int n12, int n21, int n22, int nhet, double &p11, double &p12, double &p21, double &p22)
    {
        static const int maxit = 1000;
        static const double tol = 1e-7;

        double n = n11 + n12 + n21 + n22;
        p11 = n11 / n;
        p12 = n12 / n;
        p21 = n21 / n;
        p22 = n22 / n;

        if (nhet == 0)
            return;

        auto cp11 = p11;
        auto cp12 = p12;
        auto cp21 = p21;
        auto cp22 = p22;

        auto h = nhet / n;
        auto x = h / 2;
        auto y = h - x;

        for (int i = 0; i < maxit; ++i) {
            p11 = cp11 + x;
            p12 = cp12 + y;
            p21 = cp21 + y;
            p22 = cp22 + x;
            auto z = h * p11 * p22 / (p11 * p22 + p12 * p21);
            if (std::fabs(x - z) < tol)
                break;
            x = z;
            y = h - x;
        }
    }

    // D' 95% confidence interval estimate
    int calc_dprime_CI(int ploidy, const std::vector<allele_t> &x, const std::vector<allele_t> &y, int &lower, int &upper)
    {
        size_t n = x.size() / ploidy;
        int f[3][3] = { { 0,0,0 },{ 0,0,0 },{ 0,0,0 } };

        if (ploidy == 1) {
            for (size_t i = 0; i < n; ++i) {
                if (x[i] && y[i]) {
                    auto a = x[i] * 2 - 2;
                    auto b = y[i] * 2 - 2;
                    f[a][b] += 1;
                }
            }
        }
        else if (ploidy == 2) {
            for (size_t i = 0; i < n; ++i) {
                auto i1 = i * 2, i2 = i * 2 + 1;
                if (x[i1] && x[i2] && y[i1] && y[i2]) {
                    auto a = x[i1] + x[i2] - 2;
                    auto b = y[i1] + y[i2] - 2;
                    f[a][b] += 1;
                }
            }
        }

        int n11 = f[0][0] * 2 + f[0][1] + f[1][0];
        int n12 = f[0][2] * 2 + f[0][1] + f[1][2];
        int n21 = f[2][0] * 2 + f[1][0] + f[2][1];
        int n22 = f[2][2] * 2 + f[2][1] + f[1][2];
        int nhet = f[1][1];

        double nn = n11 + n12 + n21 + n22 + nhet * 2;
        if (nn < 4)
            return 1;

        double p11 = 0.0;
        double p12 = 0.0;
        double p21 = 0.0;
        double p22 = 0.0;

        calc_hap_prob_EM(n11, n12, n21, n22, nhet, p11, p12, p21, p22);

        if (n11 > n22) {
            std::swap(n11, n22);
            std::swap(n12, n21);
            std::swap(p11, p22);
            std::swap(p12, p21);
        }

        if (n11 > n12 || n11 > n21) {
            if (n12 < n21) {
                std::swap(n11, n12);
                std::swap(n21, n22);
                std::swap(p11, p12);
                std::swap(p21, p22);
            }
            else {
                std::swap(n11, n21);
                std::swap(n12, n22);
                std::swap(p11, p21);
                std::swap(p12, p22);
            }
        }

        auto p1x = (n11 + n12 + nhet) / nn;
        auto p2x = 1 - p1x;
        auto px1 = (n11 + n21 + nhet) / nn;
        auto px2 = 1 - px1;

        auto D = p11 - p1x*px1;
        auto Dmax = D < 0.0 ? std::min(p1x*px1, p2x*px2) : std::min(p1x*px2, p2x*px1);

        if (p11 < 1e-10) p11 = 1e-10;
        if (p12 < 1e-10) p12 = 1e-10;
        if (p21 < 1e-10) p21 = 1e-10;
        if (p22 < 1e-10) p22 = 1e-10;
        auto LL1 = n11*log(p11) + n12*log(p12) + n21*log(p21) + n22*log(p22) + nhet*log(p11*p22 + p12*p21);

        if (D < 0.0) {
            std::swap(p1x, p2x);
            std::swap(n11, n21);
            std::swap(n12, n22);
        }

        double tp = 0.0;
        double ls[101];
        auto Dstep = Dmax / 100;

        for (int i = 0; i < 101; ++i) {
            auto q11 = i*Dstep + p1x*px1;
            auto q12 = p1x - q11;
            auto q21 = px1 - q11;
            auto q22 = p2x - q21;
            if (i == 100) {
                if (q11 < 1e-10) q11 = 1e-10;
                if (q12 < 1e-10) q12 = 1e-10;
                if (q21 < 1e-10) q21 = 1e-10;
                if (q22 < 1e-10) q22 = 1e-10;
            }
            auto LL2 = n11*log(q11) + n12*log(q12) + n21*log(q21) + n22*log(q22) + nhet*log(q11*q22 + q12*q21);
            auto prob = std::exp(LL2 - LL1);
            ls[i] = prob;
            tp += prob;
        }

        double sp = 0.0;
        auto tp5 = tp * 0.05;
        for (int i = 0; i < 101; ++i) {
            sp += ls[i];
            if (sp > tp5 && sp - ls[i] < tp5) {
                lower = i - 1;
                break;
            }
        }

        sp = 0.0;
        for (int i = 100; i >= 0; --i) {
            sp += ls[i];
            if (sp > tp5 && sp - ls[i] < tp5) {
                upper = i + 1;
                break;
            }
        }

        return 0;
    }

    void classify_dprime_CI(int llim, int ulim, bool &strong, bool &recomb, bool &strong2, bool &strong34)
    {
        if (ulim >= par.ulim) {
            if (llim > par.llim)
                strong = true;
            if (llim > 80)
                strong2 = strong34 = true;
            else if (llim > 50)
                strong34 = true;
        }

        if (ulim < par.recomb)
            recomb = true;
    }

    double test_block_Gabriel(int x, int y, int n, const std::vector< std::vector<bool> > &sr, const std::vector< std::vector<bool> > &ss)
    {
        int s = 0, r = 0;
        for (int i = x; i <= y; ++i) {
            for (int j = i + 1; j <= y; ++j) {
                r += sr[j][i];
                if (n > 4)
                    s += sr[i][j];
                else if (n == 2)
                    s += ss[i][j];
                else if (n == 3 || n == 4)
                    s += ss[j][i];
            }
        }

        int t = s + r;
        if (n == 2) {
            if (t < 1)
                return -1;
        }
        else if (n == 3) {
            if (t < 3)
                return -2;
        }
        else {
            if (t < 6)
                return -3;
        }

        static const double eps = std::numeric_limits<double>::epsilon();
        double inform = (double)s / t;
        if (inform - par.inform > eps)
            return inform;

        return -4;
    }

    void find_block_Gabriel(const Genotype &gt, const std::vector<int> &snps, std::vector< std::pair<int, int> > &ppos)
    {
        // !!! polymorphic & bi-allelic only !!!

        int n = snps.size();
        int w = par.batch;

        for (int i = 0; i < n; ++i) {
            auto pos1 = gt.pos[snps[i]];
            for (int j = i + 1; j < n; ++j) {
                auto pos2 = gt.pos[snps[j]];
                if (pos2 <= pos1) {
                    std::cerr << "ERROR: chromosome positions must be in ascending order: " << pos1 << " " << pos2 << "\n";
                    return;
                }
                auto dist = pos2 - pos1;
                if (dist <= par.maxlen && w < (j - i + 1))
                    w = j - i + 1;
            }
        }

        if (w > n)
            w = n;

        std::vector<std::vector<int>> ci(w, std::vector<int>(w, -1));
        std::vector<std::vector<bool>> sr(w, std::vector<bool>(w, false));
        std::vector<std::vector<bool>> ss(w, std::vector<bool>(w, false));

        int llim = -1, ulim = -1;
        bool strong = false, recomb = false, strong2 = false, strong34 = false;

        int ploidy = gt.ploidy;
        std::vector<Block> blks;

        int start = 0;
        std::cerr << "INFO: " << start + 1 << " - " << start + w << "\n";

        for (;;) {
            for (int i = 0; i < w; ++i) {
                auto x = snps[start + i];
                auto pos1 = gt.pos[x];
                for (int j = i + 1; j < w; ++j) {
                    auto y = snps[start + j];
                    auto pos2 = gt.pos[y];
                    if (pos2 - pos1 > par.maxlen)
                        break;

                    llim = ulim = -1;
                    calc_dprime_CI(ploidy, gt.dat[x], gt.dat[y], llim, ulim);

                    strong = recomb = strong2 = strong34 = false;
                    classify_dprime_CI(llim, ulim, strong, recomb, strong2, strong34);

                    ci[j][i] = llim;
                    ci[i][j] = ulim;
                    sr[i][j] = strong;
                    sr[j][i] = recomb;
                    ss[i][j] = strong2;
                    ss[j][i] = strong34;
                }
            }

            for (int i = 0; i < w; ++i) {
                auto pos1 = gt.pos[snps[start + i]];
                for (int j = i + 1; j < w; ++j) {
                    auto pos2 = gt.pos[snps[start + j]];
                    auto dist = pos2 - pos1;
                    if (dist > par.maxlen)
                        break;
                    if (ci[j][i] < par.llim || ci[i][j] < par.ulim)
                        continue;
                    int q = j - i + 1;
                    if ((q == 2 && dist > 20000) || (q == 3 && dist > 30000))
                        continue;
                    float inform = test_block_Gabriel(i, j, q, sr, ss);
                    if (inform > 0.0)
                        blks.emplace_back(Block{ start + i, start + j, dist, inform });
                }
            }

            if (start + w >= n)
                break;

            auto pos2 = gt.pos[snps[start + w]];
            for (int i = 1; i < w; ++i) {
                auto pos1 = gt.pos[snps[++start]];
                auto dist = pos2 - pos1;
                if (pos2 - pos1 <= par.maxlen)
                    break;
            }

            if (start + w >= n)
                w = n - start;

            std::cerr << "INFO: " << start + 1 << " - " << start + w << "\n";
        }

        auto cmp = [&](const Block &a, const Block &b) {
            if (a.length > b.length) return true;
            if (a.length < b.length) return false;
            if ((b.first > a.first && b.first < a.last) || (b.last > a.first && b.last < a.last)) {
                if (a.inform > b.inform) return true;
                if (a.inform < b.inform) return false;
                int a1 = -1, a2 = -1, b1 = -1, b2 = -1;
                calc_dprime_CI(ploidy, gt.dat[snps[a.first]], gt.dat[snps[a.last]], a1, a2);
                calc_dprime_CI(ploidy, gt.dat[snps[b.first]], gt.dat[snps[b.last]], b1, b2);
                if (a1 > b1) return true;
                if (a1 < b1) return false;
            }
            return a.first < b.first;
        };

        std::sort(blks.begin(), blks.end(), cmp);

        std::vector<bool> inblock(n + 1, false);

        for (auto &e : blks) {
            auto first = e.first, last = e.last;
            if (inblock[first] || inblock[last])
                continue;
            ppos.emplace_back(gt.pos[snps[first]], gt.pos[snps[last]]);
            for (auto i = first; i <= last; ++i)
                inblock[i] = true;
        }

        std::sort(ppos.begin(), ppos.end());
    }

    int read_block(const std::string &filename, std::vector<std::string> &chr, std::vector<int> &pos1, std::vector<int> &pos2)
    {
        std::ifstream ifs(filename);
        if (!ifs) {
            std::cerr << "ERROR: can't open file: " << filename << "\n";
            return 1;
        }

        size_t ln = 0;
        for (std::string line; std::getline(ifs, line); ) {
            ++ln;
            line.erase(line.find_last_not_of("\r\n") + 1);

            std::vector<std::string> vs;
            split(line, " \t", vs);
            if (vs.empty())
                continue;

            if (vs.size() < 3) {
                std::cerr << "ERROR: expected at least 3 columns for each line of block file\n";
                return 1;
            }

            auto start = std::stoi(vs[1]);
            auto stop = std::stoi(vs[2]);

            if (stop <= start) {
                std::cerr << "ERROR: invalid block: " << vs[0] << " " << vs[1] << " " << vs[2] << "\n";
                return 1;
            }

            chr.push_back(vs[0]);
            pos1.push_back(start);
            pos2.push_back(stop);
        }

        return 0;
    }

    std::vector< std::vector<int> > index_snps(const Genotype &gt, const std::vector<std::string> &chr)
    {
        int nchr = chr.size();

        std::map<std::string, int> mapchr;
        for (int i = 0; i < nchr; ++i)
            mapchr[chr[i]] = i;

        size_t m = gt.loc.size();

        std::vector< std::vector<int> > idx(nchr);
        for (size_t i = 0; i < m; ++i) {
            auto j = mapchr[gt.chr[i]];
            idx[j].push_back(i);
        }

        return idx;
    }

    size_t count_match(const std::vector<allele_t> &x, const std::vector<allele_t> &y)
    {
        size_t n = x.size();
        size_t c = 0;
        for (size_t i = 0; i < n; ++i) {
            if (x[i] == y[i])
                ++c;
        }
        return c;
    }

    void group_snps(const Genotype &gt, std::vector<int> &snps, Genotype &ogt)
    {
        int n = gt.ind.size();
        bool haploid = gt.ploidy == 1;

        std::vector< std::vector<allele_t> > dat;
        for (int i = 0; i < n; ++i) {
            std::vector<allele_t> v1, v2;
            for (auto j : snps) {
                if (haploid) {
                    v1.push_back(gt.dat[j][i]);
                }
                else {
                    v1.push_back(gt.dat[j][i * 2]);
                    v2.push_back(gt.dat[j][i * 2 + 1]);
                }
            }
            dat.push_back(v1);
            if (!v2.empty())
                dat.push_back(v2);
        }

        auto haps = unique(dat);
        int m = haps.size();

        std::vector<int> freq;
        for (auto &e : haps) {
            auto c = std::count(dat.begin(), dat.end(), e);
            freq.push_back(c);
        }

        auto ord = order(freq);
        std::reverse(ord.begin(), ord.end());
        subset(haps, ord).swap(haps);
        subset(freq, ord).swap(freq);

        auto mac = static_cast<int>(std::ceil(par.maf * n * gt.ploidy));
        int na = std::count_if(freq.begin(), freq.end(), [mac](int a) { return a >= mac; });

        std::vector<int> hc(m);
        std::iota(hc.begin(), hc.end(), 0);
        for (int i = na; i < m; ++i) {
            std::vector<int> ns;
            for (int j = 0; j < na; ++j) {
                int c = count_match(haps[i], haps[j]);
                ns.push_back(c);
            }
            int wh = index(ns, *std::max_element(ns.begin(), ns.end()));
            hc[i] = wh;
        }

        auto start = gt.pos[snps[0]];
        auto stop = start;
        for (auto j : snps) {
            auto pos = gt.pos[j];
            if (pos < start)
                start = pos;
            else if (pos > stop)
                stop = pos;
        }

        auto loc = gt.chr[snps[0]] + "_LDB_" + std::to_string(start) + "_" + std::to_string(stop);
        ogt.loc.push_back(loc);
        ogt.chr.push_back(gt.chr[snps[0]]);
        ogt.pos.push_back(start);

        std::vector<allele_t> v;
        if (haploid) {
            for (int i = 0; i < n; ++i) {
                auto wh = index(haps, dat[i]);
                auto code = static_cast<allele_t>(hc[wh] + 1);
                v.push_back(code);
            }
        }
        else {
            for (int i = 0; i < n; ++i) {
                auto wh = index(haps, dat[i * 2]);
                auto code = static_cast<allele_t>(hc[wh] + 1);
                v.push_back(code);
                wh = index(haps, dat[i * 2 + 1]);
                code = static_cast<allele_t>(hc[wh] + 1);
                v.push_back(code);
            }
        }
        if (na == 0)
            std::fill(v.begin(), v.end(), 0);
        ogt.dat.push_back(v);

        std::vector<std::string> allele;
        for (int i = 0; i < na; ++i) {
            std::string si;
            auto p = snps.size();
            for (size_t j = 0; j < p; ++j) {
                auto jj = snps[j];
                auto a = haps[i][j];
                if (a)
                    si.append(gt.allele[jj][a - 1]);
                else
                    si.push_back('N');
            }
            allele.push_back(si);
        }
        ogt.allele.push_back(allele);
    }

    void recode_save_allele(Genotype &gt)
    {
        std::ofstream ofs(par.out + ".allele");

        if (!ofs)
            std::cerr << "ERROR: can't open file for writing: " << par.out << ".allele\n";
        else
            ofs << "Locus\tCode\tAllele\n";

        auto m = gt.loc.size();
        for (size_t i = 0; i < m; ++i) {
            auto n = gt.allele[i].size();
            for (size_t j = 0; j < n; ++j) {
                if (ofs)
                    ofs << gt.loc[i] << "\t" << j << "\t" << gt.allele[i][j] << "\n";
                gt.allele[i][j] = std::to_string(j);
            }
        }
    }

} // namespace

int snpldb(int argc, char *argv[])
{
    CmdLine cmd("snpldb [options]");

    cmd.add("--vcf", "VCF file", "");
    cmd.add("--blk", "predefined block file", "");
    cmd.add("--out", "output file", "snpldb.out");
    cmd.add("--maf", "minimum minor haplotype frequency", "0.01");
    cmd.add("--maxlen", "maximum length of blocks", "200000");

    cmd.add("--llim", "lower limit CI for strong LD", "70");
    cmd.add("--ulim", "upper limit CI for string LD", "98");
    cmd.add("--recomb", "upper limit CI for strong recombination", "90");
    cmd.add("--inform", "minimum fraction of informative strong LD", "0.95");

    cmd.add("--batch", "number of SNPs in a batch", "20000");

    if (argc < 2) {
        cmd.help();
        return 1;
    }

    cmd.parse(argc, argv);

    par.vcf = cmd.get("--vcf");
    par.blk = cmd.get("--blk");
    par.out = cmd.get("--out");

    par.maf = std::stod(cmd.get("--maf"));
    par.maxlen = std::stoi(cmd.get("--maxlen"));

    par.llim = std::stoi(cmd.get("--llim"));
    par.ulim = std::stoi(cmd.get("--ulim"));
    par.recomb = std::stoi(cmd.get("--recomb"));
    par.inform = std::stod(cmd.get("--inform"));

    par.batch = std::stoi(cmd.get("--batch"));

    std::cerr << "INFO: reading genotype file...\n";
    Genotype gt;
    if (read_vcf(par.vcf, gt) != 0)
        return 2;
    std::cerr << "INFO: " << gt.ind.size() << " individuals, " << gt.loc.size() << " loci\n";

    std::vector<std::string> blk_chr;
    std::vector<int> blk_pos1, blk_pos2;

    if (!par.blk.empty())
        read_block(par.blk, blk_chr, blk_pos1, blk_pos2);

    auto chrid = stable_unique(gt.chr);
    int nchr = chrid.size();
    auto snps = index_snps(gt, chrid);

    if (blk_chr.empty()) {
        for (int i = 0; i < nchr; ++i) {
            std::cerr << "INFO: finding blocks on chromosome: " << chrid[i] << "\n";
            std::vector<std::pair<int, int>> ppos;
            find_block_Gabriel(gt, snps[i], ppos);
            blk_chr.insert(blk_chr.end(), ppos.size(), chrid[i]);
            for (auto &e : ppos) {
                blk_pos1.push_back(e.first);
                blk_pos2.push_back(e.second);
            }
        }
    }

    auto m = gt.loc.size();
    auto nb = blk_chr.size();

    Genotype ldb;
    std::vector<bool> inblock(m, false);

    std::vector<int> blk_length(nb, 0);
    std::vector<int> blk_size(nb, 0);

    for (int i = 0; i < nchr; ++i) {
        Genotype gtchr;
        for (size_t k = 0; k < nb; ++k) {
            if (blk_chr[k] != chrid[i])
                continue;
            std::vector<int> idx;
            for (auto j : snps[i]) {
                if (gt.pos[j] < blk_pos1[k] || gt.pos[j] > blk_pos2[k])
                    continue;
                inblock[j] = true;
                idx.push_back(j);
            }
            blk_length[k] = blk_pos2[k] - blk_pos1[k];
            blk_size[k] = idx.size();
            if (idx.empty()) {
                std::cerr << "WARNING: no SNPs were found in block: " << chrid[i] << " " << blk_pos1[k] << " " << blk_pos2[k] << "\n";
                continue;
            }
            group_snps(gt, idx, gtchr);
        }

        for (auto j : snps[i]) {
            if (inblock[j]) {
                auto found = std::find(gtchr.pos.begin(), gtchr.pos.end(), gt.pos[j]);
                if (found != gtchr.pos.end()) {
                    auto wh = found - gtchr.pos.begin();
                    ldb.loc.push_back(gtchr.loc[wh]);
                    ldb.chr.push_back(gtchr.chr[wh]);
                    ldb.pos.push_back(gtchr.pos[wh]);
                    ldb.dat.push_back(gtchr.dat[wh]);
                    ldb.allele.push_back(gtchr.allele[wh]);
                }
            }
            else {
                ldb.loc.push_back(gt.loc[j]);
                ldb.chr.push_back(gt.chr[j]);
                ldb.pos.push_back(gt.pos[j]);
                ldb.dat.push_back(gt.dat[j]);
                ldb.allele.push_back(gt.allele[j]);
            }
        }
    }

    std::ofstream ofs(par.out + ".block");

    if (!ofs)
        std::cerr << "ERROR: can't open file for writing: " << par.out << ".block\n";
    else {
        ofs << "Chromosome\tStart\tStop\tLength\tSNPs\n";
        for (size_t i = 0; i < nb; ++i)
            ofs << blk_chr[i] << "\t" << blk_pos1[i] << "\t" << blk_pos2[i] << "\t" << blk_length[i] << "\t" << blk_size[i] << "\n";
    }

    ldb.ind = gt.ind;
    ldb.ploidy = gt.ploidy;

    recode_save_allele(ldb);

    if (write_vcf(ldb, par.out + ".vcf") != 0)
        return 3;

    return 0;
}
