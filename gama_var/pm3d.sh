load '/home/shibu/gnuplot-palettes-master/viridis.pal'
#load '../../../gnuplot-palettes-master/rdgy.pal'
set palette negative
set xlabel 'Time' font 'Arial,30'
set ylabel 'Space' font 'Arial,30'
set yrange [0.2:]
unset key
#set title "Spatio-temoral profile of the Xist-PRC2 complex in 1 dimension for single stimulation \n n_{hill}= 4, {/Symbol d} = 0.4,  {/Symbol a}= 0.1, {/Symbol b}= 0.1, n_{stim}=2" offset 0,1
#set title "Spatio-temoral profile of the Xist-PRC2 complex in 1 dimension \n n_{hill}= 4, {/Symbol d} = 1.10,  {/Symbol a}= 0.05, {/Symbol b}= 0.1, n_{stim}=1, stim_{amp}=1.0, D_{u_1}=0.025, D_{u_2}=10.0" offset 0,1
set label '' center font 'Arial' offset 2,-2.7
set style fill transparent solid 1.0
#---------------------
# 2125 pixel = 90 mm at 600 dpi
# 1062 pixel = 90 mm at 300 dpi
# 354 pixel = 90 mm at 100 dpi
#unset colorbox
set xtics rotate by 30 right
wid=1000
hei=1000
set terminal pngcairo size wid,hei enhanced font 'Arial,30'
set output 'output.png'
set pm3d map interpolate 10,10
sp 'fort.600'

#---------------------

#set terminal svg size 2400,2400 enhanced mouse font 'Arial,60'
#set output 'output.svg'
#set pm3d map interpolate 2,2
#sp 'fort.600'

