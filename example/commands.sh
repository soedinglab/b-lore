#!/bin/bash

BLORE='../bin/blore'
POP1_IN="input/pop1"
POP2_IN="input/pop2"
POP3_IN="input/pop3"
POP1_OUT="output/pop1reg"
POP2_OUT="output/pop2reg"
POP3_OUT="output/pop3reg"

$BLORE --summary --gen ${POP1_IN}/Locus.*.gen --sample ${POP1_IN}/pop1.sample --out ${POP1_OUT} --pca 3 --regoptim
$BLORE --summary --gen ${POP2_IN}/Locus.*.gen --sample ${POP2_IN}/pop2.sample --out ${POP2_OUT} --pca 3 --regoptim
$BLORE --summary --gen ${POP3_IN}/Locus.*.gen --sample ${POP3_IN}/pop3.sample --out ${POP3_OUT} --pca 3 --regoptim

$BLORE --meta --statinfo ${POP1_OUT}/pop1reg_summary ${POP2_OUT}/pop2reg_summary ${POP3_OUT}/pop3reg_summary \
              --feature input/functional_annotation/Locus.*.feat \
              --out output/meta \
              --prefix bvslr \
              --zmax 2 \
              --muvar \
              --params 0.0005 0 0.001 0.001
