#!/bin/bash

echo 0 > counter.tmp
declare -i num=0
offset=0
for halo; do
	if [ -f $halo ] ; then
		if [ $num -eq 0 ] ; then
			head -n20 $halo > head.tmp
		fi
		head $halo | awk 'NR == 9 { print $4; }' >> counter.tmp
		(( num++ ))
	else
		offset=$halo
	fi
done
head -n2 head.tmp
echo \#Number of halo output files: $num
head -n8 head.tmp | tail -n5
awk '{ s += $1; } END { print "#Total particles processed:", s; }' counter.tmp
tail -n11 head.tmp
rm counter.tmp head.tmp
for halo; do
	if [ -f $halo ] ; then
		awk 'NR > 20 { print $1, $2, $3, $4, $5, $6, $7, $8, $9-'$offset', $10-'$offset', $11-'$offset', $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44, $45, $46, $47, $48, $49, $50, $51, $52, $53, $54; }' $halo
	fi
done
exit 0
