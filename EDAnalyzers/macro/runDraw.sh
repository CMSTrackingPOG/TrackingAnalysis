#!/bin/env bash

#mode='pv'
mode='ippt'
#input='pv.json'
input='ip.json'

python draw.py \
--mode=${mode} \
--input=${input}
