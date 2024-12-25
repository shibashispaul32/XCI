rm fort.* -f
wait
ifx xci-combined.f90 -o a.out
wait
./a.out &
wait
gnuplot pm3d.sh
