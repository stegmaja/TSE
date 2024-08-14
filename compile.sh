#!/bin/bash
mkdir plots
mkdir data
cd mobse/src
make mobse
cd ../..
pip install -r requirements.txt