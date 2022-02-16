#!/bin/bash

echo "Enter the filename with .txt at the end"
read filename

echo "Check Sum" >> $filename

for i in *fastq.gz;
do 
echo $i
echo `md5sum $i` >> $filename
done

