#!/bin/bash
set -e

rm -rf build deploy

cmake -B build 
cmake --build build --parallel

mkdir deploy
cp ./build/freight_graphs deploy/
rm -rf build
