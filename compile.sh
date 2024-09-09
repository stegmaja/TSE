#!/bin/bash
mkdir plots
mkdir data
rm -rf mobse
git clone git@gitlab.com:mobse/source-code.git mobse
cd mobse
patch -p1 < ../mobse.patch
cd src
make mobse
cd ../..
pip install -r requirements.txt
