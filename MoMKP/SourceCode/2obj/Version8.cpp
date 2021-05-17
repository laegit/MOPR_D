#ifndef VERSION8_CPP
#define VERSION8_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2018)
	
This code implements the first PathRelinking version for Multiobjective Multidimensional Knapsack problem 

hypervolume-based selection 

baseado no trabalho do anyTimePLS

apenas para 2 objetivos

=====================================================================================
*/





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
	Knapsack *aux, *aux2;
	

	double ohvc(double obj1_s1, double obj2_s1, double obj1_s2, double obj2_s2){
		return (obj1_s1 - obj1_s2)*(obj2_s2 - obj2_s1);
	}

	bool insertPareto(Knapsack *aux, std::list<int> &indicesEmAux, int i){
		list<int>::iterator it1 = indicesEmAux.begin();
		list<list<int>::iterator> remover; 
		*aux2 = *aux;
		bool domina = false;

		
		if (aux2->at(i)==false) // insere
			aux2->insert_item(i, my_momkp->profit_.at(i), my_momkp->weight_.at(i));
		else  // remove
			aux2->remove_item(i, my_momkp->profit_.at(i), my_momkp->weight_.at(i));
		
		while (it1!=indicesEmAux.end()){
			if (aux->at(*it1)==false) // insere
				aux->insert_item(*it1, my_momkp->profit_.at(*it1), my_momkp->weight_.at(*it1));
			else  // remove
				aux->remove_item(*it1, my_momkp->profit_.at(*it1), my_momkp->weight_.at(*it1));
			

			if (*aux2 >> *aux)
				remover.push_back(it1);
			if (*aux >> *aux2)
				domina = true; // return false;
			


			// volta ao que era antes
			if (aux->at(*it1)==false) // insere (antes era true)
				aux->insert_item(*it1, my_momkp->profit_.at(*it1), my_momkp->weight_.at(*it1));
			else  // remove (antes era false)
				aux->remove_item(*it1, my_momkp->profit_.at(*it1), my_momkp->weight_.at(*it1));
			
			if (domina) return false;

			it1++;
		}

		list< list<int>::iterator >::iterator j = remover.begin();
		while (j != remover.end()) {
			indicesEmAux.erase( *j );
			j++;
		}
		indicesEmAux.push_front(i);

		if (aux2->at(i)==false) // insere
			aux2->insert_item(i, my_momkp->profit_.at(i), my_momkp->weight_.at(i));
		else  // remove
			aux2->remove_item(i, my_momkp->profit_.at(i), my_momkp->weight_.at(i));
		
		return true;
	}
	
	double getOHI(Knapsack *aux, std::list<int> &indicesEmAux, int i){
		double min = 1e9;
		double max = -1;
		bool mod_sup= false, mod_inf = false;
		double s_sup_obj1, s_sup_obj2, s_inf_obj1, s_inf_obj2 ;
		list<int>::iterator it1 = indicesEmAux.begin();

		*aux2 = *aux;

		if (aux2->at(i)==false) // insere
			aux2->insert_item(i, my_momkp->profit_.at(i), my_momkp->weight_.at(i));
		else  // remove
			aux2->remove_item(i, my_momkp->profit_.at(i), my_momkp->weight_.at(i));
		
		while (it1!=indicesEmAux.end()){

			if (aux->at(*it1)==false) // insere
				aux->insert_item(*it1, my_momkp->profit_.at(*it1), my_momkp->weight_.at(*it1));
			else  // remove
				aux->remove_item(*it1, my_momkp->profit_.at(*it1), my_momkp->weight_.at(*it1));
			

			//OHI
			if (aux->profit(1)>aux2->profit(1)){
				if (aux->profit(1)<min){
					min = aux->profit(1);
					s_sup_obj1 = aux->profit(0);
					s_sup_obj2 = aux->profit(1);
					mod_sup = true;
				}
			} else if (aux->profit(1) < aux2->profit(1)){
				if (aux->profit(1)>max){
					max = aux->profit(1);
					s_inf_obj1 = aux->profit(0);
					s_inf_obj2 = aux->profit(1);
					mod_inf = true;
				}
			}

			// volta ao que era antes
			if (aux->at(*it1)==false) // insere (antes era true)
				aux->insert_item(*it1, my_momkp->profit_.at(*it1), my_momkp->weight_.at(*it1));
			else  // remove (antes era false)
				aux->remove_item(*it1, my_momkp->profit_.at(*it1), my_momkp->weight_.at(*it1));

			it1++;
		}

		
		double obj1 = aux2->profit(0);
		double obj2 = aux2->profit(1);
		if (aux2->at(i)==false) // insere
			aux2->insert_item(i, my_momkp->profit_.at(i), my_momkp->weight_.at(i));
		else  // remove
			aux2->remove_item(i, my_momkp->profit_.at(i), my_momkp->weight_.at(i));
		

		if (mod_inf == false){
			return 2*ohvc(s_sup_obj1,s_sup_obj2, obj1, obj2);
		} else if (mod_sup == false){
			return 2*ohvc(obj1, obj2, s_inf_obj1, s_inf_obj2);
		} else {
			return ohvc(s_sup_obj1,s_sup_obj2, obj1, obj2)+ohvc(obj1, obj2, s_inf_obj1, s_inf_obj2);
		}

	}

	void PR(ParetoSet *retorno, Knapsack *s1, Knapsack *s2, int &contAval){
		
		// para todo elemento em s1 cujo bit é diferente do bit correspondnte em s2. 
		// Ou seja, para todo elemento em s1 que nao está em s2 e vice-versa
		*aux = *s1;
		bool continua = true;
		


		do{
			std::list<int> indicesEmAux; // guarda os indices de aux que precisam ser invertidos
			for (int i=0; i<my_momkp->size_; i++){
				if (aux->at(i)!=s2->at(i)){
					bool insere = true;

					if (aux->at(i)==false){ 
					// se item i é candidato a ser inserido em aux, 
					//entao, precisa-se verificar que sua inserçao nao inviabiliza a soluaçao
						for (int j=0; j<aux->weight().size(); j++){ // j dimenssao
							if (aux->weight(j)+my_momkp->weight_.at(i).at(j) > my_momkp->capacity_.at(j))
								insere = false;
						}
					} 

					if (insere==true) {
						insertPareto(aux, indicesEmAux, i);
					}
				}
			}

			if (indicesEmAux.size()>0){
				continua = true;
				list<int>::iterator it1 = indicesEmAux.begin();
				int item = *it1; // min
				double maxOHI = getOHI(aux, indicesEmAux, item);
				it1++;
				
				while (it1!=indicesEmAux.end()){
					int auxItem = *it1;

					double newOHI = getOHI(aux, indicesEmAux, auxItem);

					if (newOHI>maxOHI){ // MAXIMIZA OHI
						maxOHI = newOHI;
						item = auxItem;
					}
					it1++;
				}
					
				
				if (aux->at(item)==false){ // insere
					aux->insert_item(item, my_momkp->profit_.at(item), my_momkp->weight_.at(item));
				}else{  // remove
					aux->remove_item(item, my_momkp->profit_.at(item), my_momkp->weight_.at(item));
				}
				retorno->adicionarSol(aux);
				contAval++; // sempre que um item é inserido ou removido da mochila, precisa somar ou subtrarir algum valor da funçao objetivo. Por isso, contabiliza-se a quantidade de avaliaçoes aqui
					
			} else continua = false;
		} while (continua==true && contAval<NUM_AVALIACOES);
	}

public:
	Version8(){
		aux = new Knapsack(my_momkp->size_, my_momkp->dimension_);
		aux2 = new Knapsack(my_momkp->size_, my_momkp->dimension_);
	}

};


#endif
