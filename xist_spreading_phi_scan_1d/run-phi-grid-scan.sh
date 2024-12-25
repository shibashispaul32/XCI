rm -rf phi1-* *.png *.csv
ifort xist-spreading.f90 -o a.out
concurr=15;k=0
phi1_seq=(0.01 0.1 0.5 1.0 5.0 10.0 100.0)
phi2_seq=(0.01 0.1 0.5 1.0 5.0 10.0 100.0)
for i in  "${phi1_seq[@]}"
do
for j in  "${phi2_seq[@]}"
	 #$(seq 0.0 0.10 0.10)
do
	mkdir phi1-$i-phi2-$j
	cd phi1-$i-phi2-$j 
	cp ../a.out ../input.in ../xist-spreading.f90 ../pm3d-xist-spreading.sh .
	echo $i >> phi1.in
	echo $j >> phi2.in
	./a.out &
	k=$(expr $k + 1)
	if [ $k == $concurr ]
	then
		wait
		k=0
	fi
        cd ../
done
done
wait
#------------------------------------------------
#              creating csv file
#------------------------------------------------
echo X > temp
cut -c 1-24 phi1-"${phi1_seq[1]}"-phi2-"${phi2_seq[1]}"/fort.600 > temp-x
cat temp temp-x > temp-xx
grep "\S" temp-xx > x
rm temp*
echo Y > temp
cut -c 25-48 phi1-"${phi1_seq[1]}"-phi2-"${phi2_seq[1]}"/fort.600 > temp-y
cat temp temp-y > temp-yy
grep "\S" temp-yy > y
rm temp*
paste -d , x y > combined.csv
rm x y
for i in  "${phi1_seq[@]}"
do
for j in  "${phi2_seq[@]}"
	 #$(seq 0.0 0.10 0.10)
do
	cd phi1-$i-phi2-$j
        cut -c 49-72 fort.600 > ../temp1
	cd ../
	echo phi1-$i-phi2-$j > temp2
	cat temp2 temp1 > temp3
	grep "\S" temp3 > temp4
	paste -d , combined.csv temp4 > temp5
	mv temp5 combined.csv
	rm temp* -f
done
done
wait
#------------------------------------------------
#           grid plotting in python
#------------------------------------------------
python grid-plot.py
#------------------------------------------------
#     individual  plotting using gnuplot
#------------------------------------------------
#rm -rf plots*
#mkdir plots-phi-var
#for i in 0.10 0.5 0.10 5.0.100.0 #$(seq 0.0 0.10 0.10)
#do
#for j in 0.10 0.5 0.10 5.0.100.0 #$(seq 0.0 0.10 0.10)
#do
#	cd phi1-$i-phi2-$j
#	gnuplot pm3d-xist-spreading.sh
#	mv output-u4.png  phi1-$i-phi2-$j.png
#	cp  phi1-$i-phi2-$j.png ../plots-phi-var/.
#        cd ../
#done
#done
#------------------------------------------------
