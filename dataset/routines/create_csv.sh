#!/bin/bash
echo fas, mat, cmap > ~/Stage/deep_domain_a/dataset/index.csv
ls -1 ~/Stage/deep_domain_a/dataset/*.fas | awk -F . '{print $1".fas""," $1".mat," $1".cmap"}' >> ~/Stage/deep_domain_a/dataset/index.csv
