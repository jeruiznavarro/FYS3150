#!/bin/bash
rm *.eps
gnuplot plots.gnu
latex2pdf project3.tex