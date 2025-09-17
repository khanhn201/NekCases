#!/bin/bash

# Input file
LOGFILE="logfile"
# Extract step and outflow_height values
# paste <(grep "Step" "$LOGFILE"  | head- 7000 | awk '{print $2}' | tr -d ',') \
#       <(grep "outflow_height" | head- 7000 | "$LOGFILE" | awk '{print $2}') > outflow_step.dat
paste <(grep "Step" "$LOGFILE" | awk '{for(i=1;i<=NF;i++) if($i=="t=") {print $(i+1)}}' | tr -d ',') \
      <(grep "outflow_height" "$LOGFILE" | awk '{print $2}') | head -n 7000 > outflow.dat

# Display plot
gnuplot -persist <<EOF
set terminal qt size 1200,800 enhanced font 'Arial,12'
set title "Outflow Height Over Time"
set xlabel "t"
set ylabel "Outflow Height"
set grid
plot "outflow.dat" using 1:2 with linespoints title "Outflow Height"
EOF
