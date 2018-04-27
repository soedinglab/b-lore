#!/bin/bash

BLORE='../bin/blore'
POP1_IN="input/pop1"
POP2_IN="input/pop2"
POP3_IN="input/pop3"
POP1_OUT="outputn/pop1reg"
POP2_OUT="outputn/pop2reg"
POP3_OUT="outputn/pop3reg"
META_OUT="outputn/meta"

#$BLORE --summary --gen ${POP1_IN}/Locus.*.gen --sample ${POP1_IN}/pop1.sample --out ${POP1_OUT} --pca 3 --reg 0.01 --regoptim
#$BLORE --summary --gen ${POP2_IN}/Locus.*.gen --sample ${POP2_IN}/pop2.sample --out ${POP2_OUT} --pca 3 --reg 0.01 --regoptim
#$BLORE --summary --gen ${POP3_IN}/Locus.*.gen --sample ${POP3_IN}/pop3.sample --out ${POP3_OUT} --pca 3 --reg 0.01 --regoptim


$BLORE --meta --input ${POP1_OUT}/pop1reg_summary.locusnames \
              --statdir ${POP1_OUT} ${POP2_OUT} ${POP3_OUT} \
              --out ${META_OUT} \
              --prefix blore \
              --zmax 2
