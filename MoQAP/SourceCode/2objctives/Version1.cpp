#ifndef VERSION1_CPP
#define VERSION1_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2019)
	
This code implements the first PathRelinking version for Multi-objective Quadratic Assignment Problem

=====================================================================================
*/

/// Version1: Todas as soluções da sequência PR(s_o, s_d) 
// são obtidas escolhendo elementos radomicamente



#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <sys/times.h>
#include <sys/types.h>
#include <iostream>
#include <math.h>       /* acos, cos */
#include <climits>


#include "PathRelinking.cpp"
#include "BoundedParetoSet.cpp"


class Version1: public PathRelinking{

private:
	Solution *aux;
	
	void PR(BoundedParetoSet *retorno, Solution *s11, Solution *s2, int &contAval){
		*aux = *s11;
		std::vector<pair<int, int> > paresTroca; // troca o first com o second


		do{
			paresTroca.clear();

			// ============================= Computa pares O(n^2) ============================= //
			for (int i=0; i<n_locations; i++){ // varre a origem (aux)
				// para cada facility que está na posicao i em aux e na posicao j em s2, tal que i!=j, entao cria-se o par (i,j) em paresTroca
				if (aux->getAssignment(i)!=s2->getAssignment(i)) { // isso assegura que a posicao de aux->getAssignment(i) em s2 não é i
					// procura-se a posicao j de aux->getAssignment(i) em s2
					for (int j = 0; j < n_locations; ++j){
						if (aux->getAssignment(i)==s2->getAssignment(j)){
							// devido o if (aux->getAssignment(i)!=s2->getAssignment(i)), pode-se assegurar que i!=j
							paresTroca.push_back(make_pair(i,j)); // O(n) insercoes
						}
					}
				}
			}

			// ============================= Critério de escolha ============================= //
			if (paresTroca.size()>0){



				int index_random = pick_a_number(0, paresTroca.size()-1);
				int i = paresTroca[index_random].first;
				int j = paresTroca[index_random].second;

				int a_i =  aux->getAssignment(i);
				aux->setAssignment_(i,aux->getAssignment(j));
				aux->setAssignment_(j,a_i);
				aux->calculaObjetivos(contAval); // contAval++
				retorno->adicionarSol(aux);


			}

		} while (paresTroca.size()>0);
		
	}

public:
	Version1(){
		aux = new Solution(n_locations,n_objetivos);

	}

};


#endif
