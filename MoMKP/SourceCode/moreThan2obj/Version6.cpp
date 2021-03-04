#ifndef VERSION6_CPP
#define VERSION6_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2018)
	
This code implements the first PathRelinking version for Multiobjective Multidimensional Knapsack problem 

=====================================================================================
*/

/// Pareto-based



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


class Version6: public PathRelinking{

private:
	Knapsack *aux, *aux2;

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
	// void retiraDominadas(Knapsack *aux, std::list<int> &indicesEmAux ){
	// 	list<int>::iterator it1 = indicesEmAux.begin();
	// 	*aux2 = *aux;
	// 	while (it1!=indicesEmAux.end()){
	// 		if (aux->at(*it1)==false) // insere
	// 			aux->insert_item(*it1, my_momkp->profit_.at(*it1), my_momkp->weight_.at(*it1));
	// 		else  // remove
	// 			aux->remove_item(*it1, my_momkp->profit_.at(*it1), my_momkp->weight_.at(*it1));
			
	// 		list<int>::iterator it2 = indicesEmAux.begin();
	// 		bool domina = false;
	// 		while (it2!=indicesEmAux.end() && domina==false){
	// 			if (aux2->at(*it2)==false) // insere
	// 				aux2->insert_item(*it2, my_momkp->profit_.at(*it2), my_momkp->weight_.at(*it2));
	// 			else  // remove
	// 				aux2->remove_item(*it2, my_momkp->profit_.at(*it2), my_momkp->weight_.at(*it2));
				
	// 			if (*aux2 >> *aux) domina=true;


	// 			if (aux2->at(*it2)==false) // insere (antes era true)
	// 				aux2->insert_item(*it2, my_momkp->profit_.at(*it2), my_momkp->weight_.at(*it2));
	// 			else  // remove
	// 				aux2->remove_item(*it2, my_momkp->profit_.at(*it2), my_momkp->weight_.at(*it2));
				

	// 			it2++;
	// 		}

	// 		// volta ao que era antes
	// 		if (aux->at(*it1)==false) // insere (antes era true)
	// 			aux->insert_item(*it1, my_momkp->profit_.at(*it1), my_momkp->weight_.at(*it1));
	// 		else  // remove (antes era false)
	// 			aux->remove_item(*it1, my_momkp->profit_.at(*it1), my_momkp->weight_.at(*it1));
			
	// 		if (domina==true){
	// 			it1 = indicesEmAux.erase(it1);
	// 		} else it1++;
	// 	}
	// }

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
				// cout<<"indicesEmAux.size() = "<<indicesEmAux.size()<<endl;
				
				continua = true;
				
				int index1 = pick_a_number(0, indicesEmAux.size()-1);
				list<int>::iterator it = indicesEmAux.begin();
				for (int cont_i = 0; it!=indicesEmAux.end() && cont_i<index1; cont_i++) it++;
				
				int item = (*it);
				if (aux->at(item)==false){ // insere
					aux->insert_item(item, my_momkp->profit_.at(item), my_momkp->weight_.at(item));
				}else  // remove
					aux->remove_item(item, my_momkp->profit_.at(item), my_momkp->weight_.at(item));
				
				retorno->adicionarSol(aux);
				contAval++; // sempre que um item é inserido ou removido da mochila, precisa somar ou subtrarir algum valor da funçao objetivo. Por isso, contabiliza-se a quantidade de avaliaçoes aqui
				// aux->print_solution(stdout);
				// if (aux->overload(my_momkp->capacity_)) cout<<"\tERRRROOOO"<<endl; // habilitar esta linha para verificar que todas as soluçoes intermediarias sao de fato viaveis 
			} else continua = false;
		} while (continua==true && contAval<NUM_AVALIACOES);
	}

public:
	Version6(){
		aux = new Knapsack(my_momkp->size_, my_momkp->dimension_);
		aux2 = new Knapsack(my_momkp->size_, my_momkp->dimension_);
	}

};


#endif
