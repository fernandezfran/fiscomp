#!/bin/bash

echo "0.1"
cp in.lj-data.1 in.lj-data
./a.out
rm in.lj-dat
mv thermo.dat thermo.1.dat

echo "0.2"
cp in.lj-data.2 in.lj-data
./a.out
rm in.lj-dat
mv thermo.dat thermo.2.dat

echo "0.3"
cp in.lj-data.3 in.lj-data
./a.out
rm in.lj-dat
mv thermo.dat thermo.3.dat

echo "0.4"
cp in.lj-data.4 in.lj-data
./a.out
rm in.lj-dat
mv thermo.dat thermo.4.dat

echo "0.5"
cp in.lj-data.5 in.lj-data
./a.out
rm in.lj-dat
mv thermo.dat thermo.5.dat

echo "0.6"
cp in.lj-data.6 in.lj-data
./a.out
rm in.lj-dat
mv thermo.dat thermo.6.dat

echo "0.7"
cp in.lj-data.7 in.lj-data
./a.out
rm in.lj-dat
mv thermo.dat thermo.7.dat

echo "0.8"
cp in.lj-data.8 in.lj-data
./a.out
rm in.lj-dat
mv thermo.dat thermo.8.dat

echo "0.9"
cp in.lj-data.9 in.lj-data
./a.out
rm in.lj-dat
mv thermo.dat thermo.9.dat
