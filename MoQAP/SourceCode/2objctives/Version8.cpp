#ifndef VERSION8_CPP
#define VERSION8_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2019)
	
This code implements a PathRelinking version for Multi-objective Quadratic Assignment Problem

=====================================================================================
*/

/// Version8: hypervolume -- apenas para 2 obj


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



class Version8: public PathRelinking{

private:
	Solution *aux;

	double ohvc(double obj1_s1, double obj2_s1, double obj1_s2, double obj2_s2){
		return (obj1_s1 - obj1_s2)*(obj2_s2 - obj2_s1);
	}

	bool A_Domina_B(std::vector<double> &obj1, std::vector<double> &obj2){
		bool greater = false;
	    bool lesser = false;
	    // Check if objetivos_ is lesser or greater than d.objetivos_ in any dimension.
	    for (size_t i = 0; i < n_objetivos; ++i)
	    {
	      if (obj1[i] > obj2[i]) greater = true;   // obj1 perde em pelo menos um objetivo
	      else if (obj1[i] < obj2[i]) lesser = true; // obj2 ganha em pelo menos um objetivo
	    }

	    if (lesser) /// obj1 ganha em pelo menos um objetivo
	    {
	      if (greater==false) // obj1 nao perde em nenhum objetivo
	        return true; // obj1 domina obj2
	      else return false;
	    } else return false; // obj1 nao ganha em nenhum objetivo
	}

	bool insertPareto(auxUnidadePareto &nova, std::list<auxUnidadePareto> &possibilidades){
		list<auxUnidadePareto>::iterator it1 = possibilidades.begin();
		list<list<auxUnidadePareto>::iterator> remover; 
		while (it1!=possibilidades.end()){
			if (A_Domina_B(nova.objetivos, (*it1).objetivos))
				remover.push_back(it1);
			if (A_Domina_B((*it1).objetivos, nova.objetivos) || ((*it1).objetivos[0]==nova.objetivos[0] && (*it1).objetivos[1]==nova.objetivos[1]))
				return false;
			it1++;
		}

		list< list<auxUnidadePareto>::iterator >::iterator j = remover.begin();
		while (j != remover.end()) {
			possibilidades.erase( *j );
			j++;
		}
		possibilidades.push_front(nova);
		return true;

	}


  	double getOHI(auxUnidadePareto &nova, std::list<auxUnidadePareto> &possibilidades){
		list<auxUnidadePareto>::iterator it1 = possibilidades.begin();
		auxUnidadePareto s_sup,  s_inf;
		bool mod_sup= false, mod_inf = false;
		double min = 1e9;
		double max = -1;
		if (possibilidades.size()==1){
			return  1e9;
		}

		while (it1!=possibilidades.end()){
			//OHI
			if ((*it1).objetivos[1]>nova.objetivos[1]){
				if ((*it1).objetivos[1]<min){
					min = (*it1).objetivos[1];
					s_sup = *it1;
					mod_sup = true;
				}
			} else if ((*it1).objetivos[1] < nova.objetivos[1]){
				if ((*it1).objetivos[1]>max){
					max = (*it1).objetivos[1];
					s_inf = *it1;
					mod_inf = true;
				}
			}
			it1++;
		}
		if (mod_inf == false){
			return 2*ohvc(s_sup.objetivos[0],s_sup.objetivos[1], nova.objetivos[0], nova.objetivos[1]);
		} else if (mod_sup == false){
			return 2*ohvc(nova.objetivos[0], nova.objetivos[1], s_inf.objetivos[0], s_inf.objetivos[1]);
		} else {
			return ohvc(s_sup.objetivos[0],s_sup.objetivos[1], nova.objetivos[0], nova.objetivos[1])+ohvc(nova.objetivos[0], nova.objetivos[1], s_inf.objetivos[0], s_inf.objetivos[1]);
		}

	}



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


				
				std::list<auxUnidadePareto> possibilidades;

				for (int opcao = 0; opcao < paresTroca.size(); ++opcao){ // Essa linha é O(n) (no maximo, n valores em aux estarao fora do lugar). Total do laço O(qn^2), pois o delta eh O(n) 
					// cout<<"troca facilidades = "<<auxiliar->getAssignment(i)<<" "<<auxiliar->getAssignment(j)<<endl;
					
					std::vector<double> objetivosss(n_objetivos);

					for (int obj = 0; obj < n_objetivos; ++obj){
						objetivosss[obj] = aux->getObj(obj)+delta(aux->getAssignment_(),obj,paresTroca[opcao].first,paresTroca[opcao].second); // i j ou j i
					}

					auxUnidadePareto nova = auxUnidadePareto{paresTroca[opcao].first,paresTroca[opcao].second,objetivosss};
					insertPareto(nova, possibilidades);	
				
				}



				list<auxUnidadePareto>::iterator it1 = possibilidades.begin();
				double maxOHI = getOHI(*it1, possibilidades);
				auxUnidadePareto maxPossi = *it1;
				it1++;
				
				while (it1!=possibilidades.end()){

					double newOHI = getOHI(*it1, possibilidades);

					if (newOHI > maxOHI){ // MAXIMIZA OHI
						maxOHI = newOHI;
						maxPossi = *it1;
					}
					it1++;
				}

				int i = maxPossi.i;
				int j = maxPossi.j;

				int a_i =  aux->getAssignment(i);
				aux->setAssignment_(i,aux->getAssignment(j));
				aux->setAssignment_(j,a_i);
				aux->calculaObjetivos(contAval); // contAval++
				retorno->adicionarSol(aux);

			}

		} while (paresTroca.size()>0);
		
	}




public:
	Version8(){
		aux = new Solution(n_locations,n_objetivos);

	}

};


#endif
