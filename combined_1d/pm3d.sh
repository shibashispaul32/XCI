load '../../gnuplot-palettes-master/viridis.pal'
set terminal pngcairo size 400*6,400*2.4 enhanced font 'Arial,25'
#set terminal svg size 380*3,290*3 enhanced font 'Arial,25'
#et terminal epslatex size 350*3,290*3 font 'Arial,14'
set autoscale y
set sample 100,100
unset key
set style fill transparent solid 1.0
#set pm3d lighting specular 0.5
#set title 'compaction of chromatin geometry in one dimension'
#set label '{/Symbol=25 d} = 0.95' center font 'Arial-Bold,18' offset 0,-4.2
#set zrange [0:1.4]
set output 'output.png'
NOYTICS = "set format y ''; unset ylabel"
LMARGIN = "set lmargin at screen 0.10; set rmargin at screen 0.37"
MMARGIN = "set lmargin at screen 0.375; set rmargin at screen 0.645"
RMARGIN = "set lmargin at screen 0.65; set rmargin at screen 0.915"
TMARGIN = "set tmargin at screen 0.83; set bmargin at screen 0.2"
#--------------------------------------------------------------
set multiplot layout 1,3
set xtics rotate by 30 right
set ytics rotate by 30 right
#--------------------------------------------------------------
set pm3d map 
#set colorbox horizontal origin 0.5,0.5 size 0.2,0.05
set xlabel 'Time' font ',35' offset 0,-1
set ylabel 'Space' font ',35' offset -2,0
#set title 'Spatio-temoral profile of free xist molecules' font ',25'
set cbtics rotate by 30 right
set colorbox horizontal user origin 0.12, 0.95 size 0.23, 0.05
@LMARGIN
@TMARGIN
sp 'fort.800'
#--------------------------------------------------------------
set pm3d map interpolate 2,2
set zrange [0:]
set xlabel 'Time' font ',35' offset 0,-1
set colorbox horizontal user origin 0.395, 0.95 size 0.23, 0.05
@NOYTICS
@MMARGIN
@TMARGIN
set xtics 200,200,1000 rotate by 30 right
#set title 'Spatio-temoral profile of tethered xist molecules' font ',25'
sp 'fort.700'
#--------------------------------------------------------------
set palette negative
set xlabel 'Time' font ',35' offset 0,-1
set colorbox horizontal user origin 0.665, 0.95 size 0.23, 0.05
@NOYTICS
@RMARGIN
@TMARGIN
#set title 'Spatio-temoral profile of the APRC2' font ',25'
sp 'fort.600'
#--------------------------------------------------------------
unset multiplot
