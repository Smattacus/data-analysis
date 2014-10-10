#!/bin/sh
#Call this script with parameters:
#make_webm.sh <prefix> <outputwebmname>

num=`ls $1*.png | wc -l`

png2yuv -I p -f 1 -b 0 -n $num -j $1_%03d.png > temp.yuv
vpxenc --good --end-usage=vbr --passes=2 --best --target-bitrate=3000 -o $2 temp.yuv
rm temp.yuv
