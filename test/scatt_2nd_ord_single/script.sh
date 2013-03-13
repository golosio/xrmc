#!/bin/sh

TH1=45
TH2=30
IRUN=1

for fn in $(ls ../scatt_2nd_ord_cluster/template/*.dat); do
    fn1=$(basename $fn)
    cat $fn | sed "s/TH1/${TH1}/g;s/TH2/${TH2}/g;s/IRUN/${IRUN}/g" > $fn1
done
