#!/bin/bash
dir="/Volumes/SHIN/Nie/MCI/"
snr50="50/"
snr100="100/"
head1="MSC_"
apped1=".fits"
apped2=".cat"
left="-14.3"
right="-8"

for((i=6;i<=25;i++));
do
	tmp="$i"
	if [ $i -le 9 ];then
		tmp="0"${tmp}
	fi
	if [ $i -ne 10  -a  $i -ne 21 ];then
		echo $tmp
		#tmp=$dir$head1${tmp}$apped1
		#echo $tmp 
		#sex $tmp -CATALOG_NAME $tmp$apped2
		./cut ${dir} $head1${tmp}$apped1$apped2 ${dir} $head1${tmp}$apped1 $dir${snr50} ${left} ${right}
	fi
done


