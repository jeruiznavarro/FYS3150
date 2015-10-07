set term postscript eps enhanced color 14
set grid
set key top right
set xlabel "x"
set ylabel "u(x)"
set output "10reg.eps"
plot '10_reg.dat' u 1:2 w lp title "Numerical solution", '10_reg.dat' u 1:3 w lp title "Analytical solution"
set output "100reg.eps"
plot '100_reg.dat' u 1:2 w lp title "Numerical solution", '100_reg.dat' u 1:3 w lp title "Analytical solution"
set output "1000reg.eps"
plot '1000_reg.dat' u 1:2 w lp title "Numerical solution", '1000_reg.dat' u 1:3 w lp title "Analytical solution"
set output "10lu.eps"
plot '10_reg.dat' u 1:2 w lp title "Substitution solution", '10_lu.dat' u 1:2 w lp title "LU decomposition solution"
set output "100lu.eps"
plot '100_reg.dat' u 1:2 w lp title "Substitution solution", '100_lu.dat' u 1:2 w lp title "LU decomposition solution"
set output "1000lu.eps"
plot '1000_reg.dat' u 1:2 w lp title "Substitution solution", '1000_lu.dat' u 1:2 w lp title "LU decomposition solution"
set key top left
set xlabel "log_{10}(h)"
set ylabel "{/Symbol e}"
set output "error.eps"
plot "error.dat" u (log10($1)):($2) w lp title "Relative error in both methods"
