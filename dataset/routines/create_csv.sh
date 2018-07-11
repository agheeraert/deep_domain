#!/bin/bash
echo fas, mat, cmap > /home/aria/Stage/deep_domain/dataset/index.csv
ls -1 /home/aria/Stage/deep_domain/dataset/*fas | awk -F . '{print $1".fas," $1".mat," $1".cmap"}' >> ../index.csv
