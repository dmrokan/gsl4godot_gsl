#!/usr/bin/env gnuplot

set term pngcairo enh mono

unset key
set pointsize 1

## exponential fitting example

# best fit curve
A = 5.17379
lambda = 0.11104
b = 1.05283
f(x) = A * exp(-lambda*x) + b

set out "../images/fit-exp.png"
plot f(x), \
     '../examples/nlfit.txt' us 2:3:4 ls 1 w errorbars
