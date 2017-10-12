#!/bin/bash

rm -rf snpldb-$1
mkdir snpldb-$1

if [ $1 == "lnx64" ]; then

    g++ *.cpp -o snpldb-$1/snpldb -s -O2 -std=c++11 -static

elif [ $1 == "win32" ]; then

    i686-w64-mingw32-g++ *.cpp -o snpldb-$1/snpldb.exe -s -O2 -std=c++11 -static

elif [ $1 == "win64" ]; then

    x86_64-w64-mingw32-g++ *.cpp -o snpldb-$1/snpldb.exe -s -O2 -std=c++11 -static

fi
