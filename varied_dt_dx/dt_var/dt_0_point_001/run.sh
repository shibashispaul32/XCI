#rm output.png fort.600
#wait
ifx wave-pinning-pde-1d.f90 -o a.out
wait
./a.out &
wait
gnuplot pm3d.sh
#wait
#ristretto output.png
