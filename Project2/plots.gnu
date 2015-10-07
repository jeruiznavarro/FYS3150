set term postscript eps enhanced color 14
set grid
unset key
set xlabel 'n_{steps}'
set ylabel 'Transformations'
set output 'trans.eps'
plot 'transformations.dat' w lp
set key top left
set ylabel 'Time'
set logscale y
set output 'time.eps'
plot 'time.dat' w lp title 'Jacobi', 'time.dat' u 1:3 w lp title 'Householder'
unset logscale y
set key top right
set xlabel '{/Symbol r}'
set ylabel 'Probability'
set output 'first_probability.eps'
plot '0.01_no-coulomb.dat' every 2 w lp title '{/Symbol w}_r=0.01 with no coulomb interaction', '0.5_no-coulomb.dat' every 2 w lp title '{/Symbol w}_r=0.5 with no coulomb interaction', '1_no-coulomb.dat' every 2 w lp title '{/Symbol w}_r=1 with no coulomb interaction', '5_no-coulomb.dat' every 2 w lp title '{/Symbol w}_r=5 with no coulomb interaction', '0.01_coulomb.dat' every 2 w lp title '{/Symbol w}_r=0.01 with coulomb interaction', '0.5_coulomb.dat' every 2 w lp title '{/Symbol w}_r=0.5 with coulomb interaction', '1_coulomb.dat' every 2 w lp title '{/Symbol w}_r=1 with coulomb interaction', '5_coulomb.dat' every 2 w lp title '{/Symbol w}_r=5 with coulomb interaction'
set output 'second_probability.eps'
plot '0.01_no-coulomb.dat' every 2::1 w lp title '{/Symbol w}_r=0.01 with no coulomb interaction', '0.5_no-coulomb.dat' every 2::1 w lp title '{/Symbol w}_r=0.5 with no coulomb interaction', '1_no-coulomb.dat' every 2::1 w lp title '{/Symbol w}_r=1 with no coulomb interaction', '5_no-coulomb.dat' every 2::1 w lp title '{/Symbol w}_r=5 with no coulomb interaction', '0.01_coulomb.dat' every 2::1 w lp title '{/Symbol w}_r=0.01 with coulomb interaction', '0.5_coulomb.dat' every 2::1 w lp title '{/Symbol w}_r=0.5 with coulomb interaction', '1_coulomb.dat' every 2::1 w lp title '{/Symbol w}_r=1 with coulomb interaction', '5_coulomb.dat' every 2::1 w lp title '{/Symbol w}_r=5 with coulomb interaction'