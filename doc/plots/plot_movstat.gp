#!/usr/bin/env gnuplot

set term pngcairo enh col size 1600,800

set out "../images/movstat2.png"
file = "../examples/movstat2.txt"

load 'grid.cfg'
set key left

set multiplot layout 2,1

plot file us 1:2 w li lc rgb "black" ti "x(t)"

set xlabel "time"
#set yrange [0:10]

plot file us 1:3 w li ti "True sigma", \
     file us 1:4 w li ti "MAD", \
     file us 1:5 w li ti "IQR", \
     file us 1:6 w li ti "S_n", \
     file us 1:7 w li ti "Q_n", \
     file us 1:8 w li lc rgb "gray" ti "sigma"

unset multiplot
