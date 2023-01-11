#!/bin/zsh

TRUNK=$1

if [[ -z $TRUNK ]]; then
  echo "Usage: $0 foo"
  exit 1
fi

sed "/\[ atomtypes \]/i#include \"joung-cheatham-atomtype.itp\"
/\[ atomtypes \]/i#include \"tip3p-atomtype.itp\"
/\[ system \]/i#include \"joung-cheatham-molecule.itp\"
/\[ system \]/i#include \"tip3p.itp\"
" < ${TRUNK}.top > ${TRUNK}_solvated.top

