



for k in 50 100 200 300 400 500 600 700 800 900 1000
#for k in 20 30 50
do
	echo "#define NUMEROVERTICES $k" >> param_NSGAII.h
	#g++ transgenetico.cpp -ansi -Wall -O2 -o ../instancias/concave/executaveis/transgenetico$k
	#g++ transgenetico.cpp -ansi -Wall -O2 -o ../instancias/correlated/executaveis/transgenetico$k
	g++ main.cpp -std=c++11 -Wall -O3 -o Executaveis/main_$k #$15
	

done

