#!/bin/sh

./plot_cheb.gp
./plot_dwt.gp
./plot_fft.gp
./plot_histogram.gp
./plot_interp_compare.gp
./plot_min.gp
./plot_ntuples.gp
./plot_ode.gp
./plot_qrng.gp
./plot_randist.gp
./plot_siman.gp

graph -T png < ../examples/interp.txt > ../images/interp.png
graph -T png < ../examples/interpp.txt > ../images/interpp.png
