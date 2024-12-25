rm fort.* -f
wait
ifx combined.f90 -o a.out 
wait
./a.out &
wait
gnuplot 2-anim.sh
wait
