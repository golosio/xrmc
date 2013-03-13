#!/bin/sh

export RUNDIR=$1

cd $RUNDIR

TH1=45
IRUN=$2
TH2=$3

mkdir d1_${IRUN}_${TH2}
for fn in $(ls template1/*.dat); do
    fn1=$(basename $fn)
    cat $fn | sed "s/TH1/${TH1}/g;s/TH2/${TH2}/g;s/IRUN/${IRUN}/g" > d1_${IRUN}_${TH2}/$fn1
done
cd d1_${IRUN}_${TH2}
xrmc $RUNDIR/d1_${IRUN}_${TH2}/input.dat
