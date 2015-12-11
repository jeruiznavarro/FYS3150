#!/bin/bash
rm *.eps
gnuplot plots.gnu
latex2pdf project4.tex