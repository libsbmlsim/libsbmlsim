# use this script with gnuplot (ex. gnuplot sample.plt)
plot 'sample.result' using 1:2 w l lw 4 ti"CDK1"
replot 'sample.result' using 1:4 w l lw 4 ti"APC"
replot 'sample.result' using 1:3 w l lw 4 ti"Plk1"
