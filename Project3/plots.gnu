set term postscript eps enhanced color 14
set grid
set xlabel 'x'
set ylabel 'Wave function'
f(x)=exp(-4*x)
set xrange [0:3]
set output 'exp.eps'
plot f(x) title 'exp(-4x)'
set autoscale
set xlabel 'n_{points}'
set ylabel 'Integral value'
set output 'quadrature.eps'
g(x)=0.19276571095
plot 'quadrature.dat' u 1:2 w lp title 'Gauss-Legendre', 'quadrature.dat' u 1:3 w lp title 'Gauss-Laguerre', g(x) title 'Closed form value'
set logscale x
set format x "10^{%L}"
set xlabel 'MC cycles'
set output 'mc_values.eps'
plot 'mc.dat' u 1:2 w lp title 'Brute force', 'mc.dat' u 1:3 w lp title 'Importance sampling', g(x) title 'Closed form value'
set ylabel 'Standard deviation'
set format y "10^{%L}"
set logscale y
set output 'mc_sd.eps'
plot 'mc.dat' u 1:4 w lp title 'Brute force', 'mc.dat' u 1:5 w lp title 'Importance sampling'
set ylabel 'Run time'
set output 'mc_time.eps'
plot 'mc.dat' u 1:6 w lp title 'Brute force', 'mc.dat' u 1:7 w lp title 'Importance sampling'
unset logscale x
set format x "%1.0f"
set xlabel 'n_{points}'
set output 'quadrature_time.eps'
plot 'quadrature.dat' u 1:4 w lp title 'Gauss-Legendre', 'quadrature.dat' u 1:5 w lp title 'Gauss-Laguerre'