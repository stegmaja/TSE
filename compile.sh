#!/bin/bash
mkdir plots
mkdir data
pip install -r requirements.txt
git clone git@gitlab.com:mobse/source-code.git mobse
cd mobse
patch -p1 < ../mobse.patch
cd src
make mobse
cd ../..
