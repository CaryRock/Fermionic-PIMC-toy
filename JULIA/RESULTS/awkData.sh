#! /bin/sh
# Shell script to do the pimcave.py + awk thing I'd otherwise manually do
# Writes Temperature in first row: [$T   -1  -1] (since awks three columns)
#
# This will probably need to be updated to handle multiple files at once, or to
# account for some later requirement in another file

# $1 ce-lineardensity-00.500-0001-01.000-1.00000-86e3e49a-5082-11ec-3984-ffca6b93bce8.dat 
temp=$(basename $1 .dat | egrep [0-9]{2}[.][0-9]{3} -o | head -n1)
echo "$temp"

$(pimcave.py $1 | awk ' NR>2 { print $1,"\t",$2,"\t",$3 }' > awkTest.dat)
#pimcave.py ce-lineardensity-00.500-0001-01.000-1.00000-86e3e49a-5082-11ec-3984-ffca6b93bce8.dat | awk 'NR>3 { print ,t,,t, }'
