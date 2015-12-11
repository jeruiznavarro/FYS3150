set term postscript eps enhanced color 14
set grid
unset key
set output 'diffusion.eps'
set xlabel 'Temperature'
set ylabel 'Diffusion constant'
plot 'diffusion.dat' w lp
set key bottom right
set logscale
set xlabel 'Timestep'
set ylabel '{/Symbol s}_E'
set output 'sigma_e.eps'
set format x "%l*10^{%L}"
plot 'sigma_e.dat' u 1:4 w lp title 'Velocity-Verlet at ~200K', 'sigma_e.dat' u 1:2 w lp title 'Velocity-Verlet at ~350K', 'sigma_e.dat' u 1:5 w lp title 'Euler-Cromer at ~200K', 'sigma_e.dat' u 1:3 w lp title 'Euler-Cromer at ~350K'