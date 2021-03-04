#ifndef NSGAII_CPP
#define NSGAII_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2019)
	
This code implements a NSGAII algorithm for Multiobjective Multidimensional Knapsack problem


=====================================================================================
*/

/// ESTA IMPLEMENTACAO NAO UTILIZA ARQUIVO EXTERNO

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

#include "BoundedParetoSet.cpp"
#include "Knapsack.cpp"
#include "param_geral.h"


using namespace std;



int objetivoOrdenacao; // esta variavel é utilizada para designar qual objetivo será utilizado para ordenar as soluçoes
		

bool compare1(Knapsack *s1, Knapsack *s2){
	return s1->profit(objetivoOrdenacao) < s2->profit(objetivoOrdenacao);
}

bool compare2(Knapsack *s1, Knapsack *s2){
	return s1->distance > s2->distance; // maior distância
}


class nsgaii {
	private:
		Knapsack *uniao[TAMANHOPOPULACAO_NSGAII*2];
		Knapsack *populacao[TAMANHOPOPULACAO_NSGAII];
		Knapsack *novaPop[TAMANHOPOPULACAO_NSGAII]; // populaçao de descendentes gerada ao fim de uma geraçao
		Knapsack * filho;
		int N[TAMANHOPOPULACAO_NSGAII*2]; // para cada p \in P, guarda a quantidade de soluçoes que dominam p
		int intContAval;	

		// deb et al (2002)
		// P[TAMANHOPOPULACAO_NSGAII*2] = populacao corrente unido com a populaçao criada
		// tem que ser EXATAMENTE TAMANHOPOPULACAO_NSGAII*2
		// F é o vetor de fronts. Cada front (f0, f1, f2, ... ) é uma list. A quantidade maxima de fronts é TAMANHOPOPULACAO_NSGAII*2
		void fast_non_dominanted_sort(Knapsack *P[TAMANHOPOPULACAO_NSGAII*2], list<Knapsack *> F[TAMANHOPOPULACAO_NSGAII*2], int &quantFronts){
			list< Knapsack* > S[TAMANHOPOPULACAO_NSGAII*2]; // para cada p \in P, guarda a lista de soluçoes dominadas por p
			
			quantFronts = 0;
			for (int p=0; p<TAMANHOPOPULACAO_NSGAII*2; p++){
				N[p] = 0;
			}
			for (int p=0; p<TAMANHOPOPULACAO_NSGAII*2; p++){
				P[p]->posicaoListaNSGAII = p;
				for (int q=0; q<TAMANHOPOPULACAO_NSGAII*2; q++){
					if (*P[p]>>*P[q]){
						S[p].push_back(P[q]);
					} else if (*P[q]>>*P[p]){
						N[p] = N[p] + 1;
					}
				}

				if (N[p] == 0){
					P[p]->rank = 0;
					F[0].push_back(P[p]); // primeiro rank. RANK COMEÇA DO 0
				}
			}
			int i=0; 
			while (F[i].size()>0){
				for (list< Knapsack* >::iterator p = F[i].begin(); p!=F[i].end(); p++){
					int pp = (*p)->posicaoListaNSGAII;
					for (list< Knapsack* >::iterator q = S[pp].begin(); q!=S[pp].end(); q++){
						int qq = (*q)->posicaoListaNSGAII;
						N[qq] = N[qq] - 1;
						if (N[qq] == 0){ //q vai pro próximo rank
							F[i+1].push_back((*q));
							(*q)->rank = i+1;
						}
					}
				}
				i++;
			}
			quantFronts = i;

		}

		void crownding_distance_assigment(list<Knapsack *> &I){
			#define INF 1e9
			int l = I.size();
			list<Knapsack *>::iterator it = I.begin();
			while (it!=I.end()){
				(*it)->distance = 0;
				it++;
			}
			for (int m=0; m<my_momkp->dimension_; m++){ // para cada objetivo m
				objetivoOrdenacao = m;
				I.sort(compare1);
				(*I.begin())->distance = INF;
				(I.back())->distance = INF; // É ISSO MESMO: back() NAO retorna um ponteiro interator!!! 
				it = I.begin();
				it++; // começa do 1 
				for (int i = 1; i<(l-1); i++){
					it++; // pos
					double objPos = (*it)->profit(m);
					it--; // volta para o original
					it--; //previous
					double objPrev = (*it)->profit(m);
					it++; // volta para o original
					(*it)->distance = (*it)->distance + (objPos - objPrev);
					it++; // avança
				}
			}
			#undef INF
		}

		// atualiza a populaçao (global)
		void atualizaPopulacaoNSGAII(Knapsack *novaPop[TAMANHOPOPULACAO_NSGAII]){
			for (int i=0; i<TAMANHOPOPULACAO_NSGAII*2; i++){
				if (i<TAMANHOPOPULACAO_NSGAII){
					*uniao[i] = *populacao[i];
				} else {
					*uniao[i] = *novaPop[i-TAMANHOPOPULACAO_NSGAII];
				}
			}
			int sizeFront = 0;
			list<Knapsack *> F[TAMANHOPOPULACAO_NSGAII*2];
			fast_non_dominanted_sort(uniao, F, sizeFront);
			int cont = 0;
			int i = 0;
			while (cont + F[i].size() < TAMANHOPOPULACAO_NSGAII && i<sizeFront){
				for (list< Knapsack* >::iterator p = F[i].begin(); p!=F[i].end(); p++){
					*populacao[cont++] = **p;
				}
				i++;
			}
			crownding_distance_assigment(F[i]);
			F[i].sort(compare2); // ordena por CD
			for (list< Knapsack* >::iterator p = F[i].begin(); p!=F[i].end() && cont<TAMANHOPOPULACAO_NSGAII; p++){
				*populacao[cont++] = **p;
			}
		}



	public:
		
		void editaPopulacao(){


			std::vector<Knapsack *> population_greedy;
			my_momkp->generate_population(population_greedy, TAMANHOPOPULACAO_NSGAII,0, 0.1, intContAval); // 1-Random 0-Greedy
			///TODO: melhorar isso!
			///TODO: VER POR QUE O GREEDY ESTÁ GERANDO MUITAS SOLUCOES REPETIDAS.
			int ontt=0;
			for (std::vector<Knapsack *>::iterator i = population_greedy.begin(); i != population_greedy.end() && ontt<TAMANHOPOPULACAO_NSGAII; ++i){
				*populacao[ontt++] = **i;
			}
			// cout<<"populacao initial "<<endl;
			// for (int i=0; i<TAMANHOPOPULACAO_NSGAII; i++){
			// 	cout<<populacao[i]->profit(0)<<" "<<populacao[i]->profit(1)<<endl;
			// }
			// cout<<endl;
		
			// filho = my_momkp->crossover(populacao[0], populacao[1], intContAval);
			// populacao[0]->print_solution(stdout);
			// populacao[1]->print_solution(stdout);
			// filho->print_solution(stdout);
			// cout<<filho->overload(my_momkp->capacity_)<<endl;

			// populacao[0]->print_solution(stdout);
			// my_momkp->mutation(populacao[0], intContAval);
			// populacao[0]->print_solution(stdout);
		}

		nsgaii(){ 
		
			intContAval = 0;
			// filho = new Knapsack(my_momkp->size_, my_momkp->dimension_); 
			for (int i=0; i<TAMANHOPOPULACAO_NSGAII; i++){
				populacao[i] = new Knapsack(my_momkp->size_, my_momkp->dimension_);
				novaPop[i] = new Knapsack(my_momkp->size_, my_momkp->dimension_);
			}
			for (int i=0; i<TAMANHOPOPULACAO_NSGAII*2; i++){
				uniao[i] = new Knapsack(my_momkp->size_, my_momkp->dimension_); // permanente
			}
			
		}

		Knapsack **executar(){ // retorna um array de ponteiros
			
			intContAval=0;
			editaPopulacao();
			int p1,p2,p3,p4;
			int pai, mae;	
			int geracao = 0;
			while (intContAval<NUM_AVALIACOES){
				geracao++;
				// cout<<"geracao = "<<geracao<<" intContAval = "<<intContAval<<endl;
			
				for (int j=0; j<TAMANHOPOPULACAO_NSGAII; j++){ // deve-se criar TAMANHOPOPULACAO_NSGAII novos individuos
					double p = pick_a_number(0.0,1.0); //genrand64_real3();;//rg.Random();
					if (p<TAXADECRUZAMENTO){
						/*SORTEIA 4 individuos*/
						/*Faz-se o torneio binario entre eles*/
						p1 = pick_a_number(0, TAMANHOPOPULACAO_NSGAII-1);
						p2 = pick_a_number(0, TAMANHOPOPULACAO_NSGAII-1);
						p3 = pick_a_number(0, TAMANHOPOPULACAO_NSGAII-1);
						p4 = pick_a_number(0, TAMANHOPOPULACAO_NSGAII-1);

						int objjjj_ram = pick_a_number(0, my_momkp->dimension_-1);
						if(populacao[p1]->profit(objjjj_ram)>populacao[p2]->profit(objjjj_ram)){ // compete com o primeiro objetivo
							pai = p1;;
						} else {
							pai = p2;
						}

						objjjj_ram = pick_a_number(0, my_momkp->dimension_-1);
						if (populacao[p3]->profit(objjjj_ram)>populacao[p4]->profit(objjjj_ram)){ // compete com o primeiro objetivo
							mae = p3;
						} else {
							mae = p4;
						}
						
						filho = my_momkp->crossover(populacao[pai], populacao[mae],intContAval);
						// filho foi definido; Agora aplica-se mutaçao
						p = pick_a_number(0.0,1.0);
						if (p<TAXADEMUTACAO){
							my_momkp->mutation(filho,intContAval);

						}
						*novaPop[j] = *filho;

					} else {
						my_momkp->generate_solution(filho,0,intContAval); // greedy
						// my_momkp->BuildSolution(filho, 0.1, intContAval);
						*novaPop[j] = *filho;
					}
				}
				atualizaPopulacaoNSGAII(novaPop);
			}

			return populacao;
		}
};



#endif