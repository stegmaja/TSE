#!/bin/bash
mkdir plots
mkdir data
cd mobse/src
make mobse
cd ../..
conda env create -f environment.yml
conda activate tse