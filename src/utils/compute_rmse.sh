#!/bin/sh
awk '{print $1+2*$2+3*$3+4*$4+5*$5}' $1 > tmp1asdf
awk '{print $6}' $2 >tmp2asdf
paste tmp1asdf tmp2asdf > asdfgh
awk '{print $1-$2}'  asdfgh > lllll
awk 'BEGIN{n=0;s=0} {s = s+ $1*$1; n++} END{ print sqrt(s/n)}' ./lllll
