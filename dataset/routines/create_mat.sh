#!/bin/bash

cd ../aln/
for file in `ls *.aln`; do
	noext=`echo $file | awk -F . '{print $1}'`
	ccmpred $file ../mat/$noext.mat 
done
