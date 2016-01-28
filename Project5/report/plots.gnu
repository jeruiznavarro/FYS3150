set term postscript eps enhanced color 14
set grid
unset key
set xrange [0.8:2]
set yrange [-1:5]
set output 'lennard_jones.eps'
set xlabel 'Inter-particle distance (r/{/Symbol s})'
set ylabel 'Potential (V/{/Symbol e})'
f(x)=4*(1/x)**6*((1/x)**6-1)
plot f(x)
set autoscale
set key top center
set xrange [0:1000]
set output 'energy_evolution.eps'
set xlabel 'Time (fs)'
set ylabel 'Energy per particle (hartree)'
plot 'time_evolution.dat' u ($1*1e15):2 w lp title '{/Symbol D}t=1', 'time_evolution_fast.dat' u ($1*1e15):2 w lp title '{/Symbol D}t=10'
set output 'temperature_evolution.eps'
set xlabel 'Time (fs)'
set ylabel 'Temperature (K)'
plot 'time_evolution.dat' u ($1*1e15):5 w lp title '{/Symbol D}t=1', 'time_evolution_fast.dat' u ($1*1e15):5 w lp title '{/Symbol D}t=10'
unset key
set autoscale
set output 'diffusion.eps'
set xlabel 'Temperature (K)'
set ylabel 'Diffusion constant (m^2/s)'
plot 'diffusion.dat' w lp
set key bottom right
set logscale
set xlabel 'Timestep (fs)'
set ylabel '{/Symbol s}_E (hartree)'
set output 'sigma_e.eps'
plot 'sigma_e.dat' u 1:4 w lp title 'Velocity-Verlet at 200K', 'sigma_e.dat' u 1:2 w lp title 'Velocity-Verlet at 350K', 'sigma_e.dat' u 1:5 w lp title 'Euler-Cromer at 200K', 'sigma_e.dat' u 1:3 w lp title 'Euler-Cromer at 350K'