#!/bin/bash

# Input file
LOGFILE="logfile"
# Extract step and outflow_height values
paste <(grep "Step" "$LOGFILE" | awk '{print $2}' | tr -d ',') \
      <(grep "outflow_height" "$LOGFILE" | awk '{print $2}') > outflow_step.dat

# Display plot
gnuplot -persist <<EOF
set terminal qt size 1200,800 enhanced font 'Arial,12'
set title "Outflow Height Over Time"
set xlabel "Step"
set ylabel "Outflow Height"
set grid
plot "outflow_step.dat" using 1:2 with linespoints title "Outflow Height"
EOF
