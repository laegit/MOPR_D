#ifndef VERSION4_CPP
#define VERSION4_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2018)
	
This code implements the third PathRelinking version for Multiobjective Multidimensional Knapsack problem 

=====================================================================================
*/

/*
	Version4 = Escolha de elementos que otimizam uma função de escalarização
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


class Version4: public PathRelinking{

private:
	Knapsack *aux;
	std::vector<double> lambda;

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
				int item = indicesEmAux[0]; 
				double max = 0.0;
				if (aux->at(item)==false){ // simula adiçao do item 0
					for (int pl = 0; pl<my_momkp->dimension_; pl++){
						max += lambda[pl]*((double) aux->profit(pl)+my_momkp->profit_.at(item).at(pl));
					}
				} else { // simula remoçao do item 0
					for (int pl = 0; pl<my_momkp->dimension_; pl++){
						max += lambda[pl]*((double) aux->profit(pl)-my_momkp->profit_.at(item).at(pl));
					}
				}

				for (int jjt = 0; jjt<indicesEmAux.size(); jjt++){
					int auxItem = indicesEmAux[jjt];
					double jjtMAX = 0.0;

					if (aux->at(auxItem)==false){ // simula adiçao do item indicesEmAux[jjt]
						for (int pl = 0; pl<my_momkp->dimension_; pl++){
							jjtMAX += lambda[pl]*((double) aux->profit(pl)+my_momkp->profit_.at(auxItem).at(pl));
						}
					} else { // simula remoçao do item indicesEmAux[jjt]
						for (int pl = 0; pl<my_momkp->dimension_; pl++){
							jjtMAX += lambda[pl]*((double) aux->profit(pl)-my_momkp->profit_.at(auxItem).at(pl));
						}
					}

					if (jjtMAX>max){ // MAXIMIZA OBJETIVO 'obj'
						max = jjtMAX;
						item = auxItem;
					}
				}
					
				// cout<<"max = "<<max<<endl;
				if (aux->at(item)==false){ // insere
					aux->insert_item(item, my_momkp->profit_.at(item), my_momkp->weight_.at(item));
				}else{  // remove
					aux->remove_item(item, my_momkp->profit_.at(item), my_momkp->weight_.at(item));
					// cout<<"remove"<<endl;
				}
				retorno->adicionarSol(aux);
				contAval++; // sempre que um item é inserido ou removido da mochila, precisa somar ou subtrarir algum valor da funçao objetivo. Por isso, contabiliza-se a quantidade de avaliaçoes aqui
				// aux->print_solution(stdout);
				// if (aux->overload(my_momkp->capacity_)) cout<<"\tERRRROOOO"<<endl; // habilitar esta linha para verificar que todas as soluçoes intermediarias sao de fato viaveis 
			

			} else continua = false;
		} while (continua==true && contAval<NUM_AVALIACOES);
	}


public:
	Version4(std::vector<double> lambda1){
		aux = new Knapsack(my_momkp->size_, my_momkp->dimension_);
		lambda = lambda1;

		// for (int i=0; i<lambda.size(); i++){
		// 	cout<<lambda[i]<<" ";
		// }
		// cout<<endl;
	}

};


#endif
