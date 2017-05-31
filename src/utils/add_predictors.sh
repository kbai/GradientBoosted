#!/bin/bash
paste $1 $2| awk -v x=$3 -v y=$4 '{print x*$1+y*$2}'
