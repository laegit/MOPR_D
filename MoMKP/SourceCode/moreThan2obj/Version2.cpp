#ifndef VERSION2_CPP
#define VERSION2_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2018)
	
This code implements the second PathRelinking version for Multiobjective Multidimensional Knapsack problem 

=====================================================================================

///Todas as soluções da sequência PR(s_o, s_d) são obtidas otimizando um só objetivo

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


class Version2: public PathRelinking{

private:
	Knapsack *aux;
	int obj;

	void PR(ParetoSet *retorno, Knapsack *s1, Knapsack *s2, int &contAval){
		
		// para todo elemento em s1 cujo bit é diferente do bit correspondnte em s2. 
		// Ou seja, para todo elemento em s1 que nao está em s2 e vice-versa
		*aux = *s1;
		bool continua = true;
		
		do{
			std::vector<int> indicesEmAux; // guarda os indices de aux que precisam ser invertidos
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

					if (insere==true) indicesEmAux.push_back(i);
				}
			}

			if (indicesEmAux.size()>0){
				// cout<<"indicesEmAux.size() = "<<indicesEmAux.size()<<endl;
				continua = true;
				int item = indicesEmAux[0]; // min
				
				// simula adiçao o item 0
				int max = aux->profit(obj)+my_momkp->profit_.at(item).at(obj);
				if (aux->at(item)==true) // simula remoçao do item 0
					max = aux->profit(obj)-my_momkp->profit_.at(item).at(obj);

				for (int jjt = 0; jjt<indicesEmAux.size(); jjt++){
					int auxItem = indicesEmAux[jjt];
					// simula adiçao do item indicesEmAux[jjt]
					int jjtMin = aux->profit(obj)+my_momkp->profit_.at(auxItem).at(obj); 
					if (aux->at(auxItem)==true){
						jjtMin = aux->profit(obj)-my_momkp->profit_.at(auxItem).at(obj);
					}
					if (jjtMin>max){ // MAXIMIZA OBJETIVO 'obj'
						max = jjtMin;
						item = auxItem;
					}
				}
					
				
				if (aux->at(item)==false){ // insere
					aux->insert_item(item, my_momkp->profit_.at(item), my_momkp->weight_.at(item));
				}else{  // remove
					aux->remove_item(item, my_momkp->profit_.at(item), my_momkp->weight_.at(item));
				}
				retorno->adicionarSol(aux);
				contAval++; // sempre que um item é inserido ou removido da mochila, precisa somar ou subtrarir algum valor da funçao objetivo. Por isso, contabiliza-se a quantidade de avaliaçoes aqui
				// aux->print_solution(stdout);
				// if (aux->overload(my_momkp->capacity_)) cout<<"\tERRRROOOO"<<endl; // habilitar esta linha para verificar que todas as soluçoes intermediarias sao de fato viaveis 
			} else continua = false;
		} while (continua==true && contAval<NUM_AVALIACOES);
	}


public:
	Version2(int obj1){
		aux = new Knapsack(my_momkp->size_, my_momkp->dimension_);
		obj = obj1;
	}

};


#endif
