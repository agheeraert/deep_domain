#!/bin/bash

for file in `ls /home/aria/Stage/deep_domain_a/dataset/aln/*.sto`; do
	cat $file | awk '{print $2}' | head -n -1 >> $file-aln
done

cd /home/aria/Stage/deep_domain_a/dataset 
for file in `ls *.sto-aln`; do
	echo $file | awk -F . '{print $1}' | xargs -I {} mv {}.sto-aln ../{}.aln
done
cd ..
for file in `ls *.aln`; do
	sed -i -e '1,4d' $file
done
