#!/bin/bash

head1="MSC_MS_210526082000_100000120_"
head2="MSC_MS_210526082000_100000120_"
head3="0"
apped1="*.fits"
apped2=".dat"
lpath="../TEST300/MSC_0000009/"

fget=`ls ${lpath}${apped1}`
#echo $fget
farr=(${fget})
for var in ${farr[@]}
do 
	#echo ${var}
	svar=(${var##*/})
	#split string
	echo ${svar}
	./cut ${lpath} ${svar}${apped2} ${lpath} ${svar} ../TEST300/150/
done


:<<!
for((i=6;i<=9;i++));
do
	tmp="$i"
	orders1=${head3}${tmp}${apped2}
	orders2=${head1}${head3}${tmp}${apped1}
	./cut ../${lpath}/ ${orders1} ../${lpath}/ ${orders2} ../${lpath}/08/
done
for((i=11;i<=20;i++));
do
	tmp="$i"
	orders1=${tmp}${apped2}
	orders2=${head1}${tmp}${apped1}
	./cut ../${lpath}/ ${orders1} ../${lpath}/ ${orders2} ../${lpath}/08/
done
for((i=22;i<=25;i++));
do
	tmp="$i"
	orders1=${tmp}${apped2}
	orders2=${head1}${tmp}${apped1}
	./cut ../${lpath}/ ${orders1} ../${lpath}/ ${orders2} ../${lpath}/08/
done

!