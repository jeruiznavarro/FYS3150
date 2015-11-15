set term postscript eps enhanced color 14
set grid
set xlabel 'MC cycles'
set ylabel '<E>'
set output 'energy_T1.eps'
plot 'energy_T1_disordered.dat' every 100 w lp title 'Disordered configuration', 'energy_T1_ordered.dat' every 100 w lp title 'Ordered configuration'
set output 'energy_T2.4.eps'
plot 'energy_T2.4_disordered.dat' every 100 w lp title 'Disordered configuration', 'energy_T2.4_ordered.dat' every 100 w lp title 'Ordered configuration'
set key right center
set ylabel '<|M|>'
set output 'magnetisation_T1.eps'
plot 'magnetisation_T1_disordered.dat' every 100 w lp title 'Disordered configuration', 'magnetisation_T1_ordered.dat' every 100 w lp title 'Ordered configuration'
set output 'magnetisation_T2.4.eps'
plot 'magnetisation_T2.4_disordered.dat' every 100 w lp title 'Disordered configuration', 'magnetisation_T2.4_ordered.dat' every 100 w lp title 'Ordered configuration'
set key left top
set xlabel 'Temperature'
set ylabel 'Acceptance rate'
set output 'acceptance.eps'
plot 'acceptance.dat' w lp title 'Ratio of moves accepted per spin'
set ylabel '<E>'
set output 'temperature_energy.eps'
plot '20L.dat' u 1:2 w lp title 'L=20', '40L.dat' u 1:2 w lp title 'L=40', '60L.dat' u 1:2 w lp title 'L=60', '80L.dat' u 1:2 w lp title 'L=80', '100L.dat' u 1:2 w lp title 'L=100'
set ylabel '<|M|>'
set output 'temperature_magnetisation.eps'
plot '20L.dat' u 1:5 w lp title 'L=20', '40L.dat' u 1:5 w lp title 'L=40', '60L.dat' u 1:5 w lp title 'L=60', '80L.dat' u 1:5 w lp title 'L=80', '100L.dat' u 1:5 w lp title 'L=100'
set key right top
set ylabel 'C_V'
set output 'temperature_cv.eps'
plot '20L.dat' u 1:3 w lp title 'L=20', '40L.dat' u 1:3 w lp title 'L=40', '60L.dat' u 1:3 w lp title 'L=60', '80L.dat' u 1:3 w lp title 'L=80', '100L.dat' u 1:3 w lp title 'L=100'
set ylabel '{/Symbol c}'
set output 'temperature_chi.eps'
plot '20L.dat' u 1:4 w lp title 'L=20', '40L.dat' u 1:4 w lp title 'L=40', '60L.dat' u 1:4 w lp title 'L=60', '80L.dat' u 1:4 w lp title 'L=80', '100L.dat' u 1:4 w lp title 'L=100'
set xlabel 'Energy'
set ylabel '# times energy value appears'
set output 'energy_histogram_T1.eps'
binwidth=0.01
set xrange[-794.6:-794.2]
bin(x,width)=width*floor(x/width)
plot 'energy_histogram_T1.dat' u (bin($1,binwidth)):(1.0) smooth freq with boxes title 'T=1'
set output 'energy_histogram_T24.eps'
binwidth=0.01
set xrange[-487:-480]
plot 'energy_histogram_T24.dat' u (bin($1,binwidth)):(1.0) smooth freq with boxes title 'T=2.4'