#!/bin/bash
cd ../src/
ls -l *f90 | awk '{if($NF~/f90/) print $NF}' | while read frqcase
do
  echo ${frqcase}
  ../postprocess/a.out ../tmp/${frqcase} ${frqcase}
done
