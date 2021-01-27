#!/bin/bash

echo "0.01"
cp in.lj-data01 in.lj-data
./a.out
rm in.lj-data
mv histo.dat histo01p.dat

echo "0.001"
cp in.lj-data001 in.lj-data
./a.out
rm in.lj-data
mv histo.dat histo001p.dat
