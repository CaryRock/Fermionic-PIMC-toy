#!/bin/sh

name=$1
file=$2

pimcave.py $file | awk '{ print $1,'\t',$2,'\t',$3 }' > "$name".dat
