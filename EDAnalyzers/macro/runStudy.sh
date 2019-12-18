#!/bin/env bash

input='output.root'

echo "Run PV study"
python pvStudy.py \
--input=${input}

echo "Run IP study"
python ipStudy.py \
--input=${input}
