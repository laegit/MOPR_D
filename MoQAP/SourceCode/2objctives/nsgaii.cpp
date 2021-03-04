#ifndef NSGAII_CPP
#define NSGAII_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2019)
	
This code implements a NSGAII algorithm for Multi-objective Quadratic Assignment Problem


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

#include "Solution.cpp"
#include "param_geral.h"


using namespace std;



int objetivoOrdenacao; // esta variavel é utilizada para designar qual objetivo será utilizado para ordenar as soluçoes
		

bool compare1(Solution *s1, Solution *s2){
	return s1->getObj(objetivoOrdenacao) < s2->getObj(objetivoOrdenacao);
}

bool compare2(Solution *s1, Solution *s2){
	return s1->distance > s2->distance; // maior distância
}


class nsgaii {
	private:
		Solution *uniao[TAMANHOPOPULACAO_NSGAII*2];
		Solution *populacao[TAMANHOPOPULACAO_NSGAII];
		Solution *novaPop[TAMANHOPOPULACAO_NSGAII]; // populaçao de descendentes gerada ao fim de uma geraçao
		Solution * filho;
		int N[TAMANHOPOPULACAO_NSGAII*2]; // para cada p \in P, guarda a quantidade de soluçoes que dominam p
		int intContAval;	

		// deb et al (2002)
		// P[TAMANHOPOPULACAO_NSGAII*2] = populacao corrente unido com a populaçao criada
		// tem que ser EXATAMENTE TAMANHOPOPULACAO_NSGAII*2
		// F é o vetor de fronts. Cada front (f0, f1, f2, ... ) é uma list. A quantidade maxima de fronts é TAMANHOPOPULACAO_NSGAII*2
		void fast_non_dominanted_sort(Solution *P[TAMANHOPOPULACAO_NSGAII*2], list<Solution *> F[TAMANHOPOPULACAO_NSGAII*2], int &quantFronts){
			list< Solution* > S[TAMANHOPOPULACAO_NSGAII*2]; // para cada p \in P, guarda a lista de soluçoes dominadas por p
			
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
				for (list< Solution* >::iterator p = F[i].begin(); p!=F[i].end(); p++){
					int pp = (*p)->posicaoListaNSGAII;
					for (list< Solution* >::iterator q = S[pp].begin(); q!=S[pp].end(); q++){
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

		void crownding_distance_assigment(list<Solution *> &I){
			#define INF 1e9
			int l = I.size();
			list<Solution *>::iterator it = I.begin();
			while (it!=I.end()){
				(*it)->distance = 0;
				it++;
			}
			for (int m=0; m<n_objetivos; m++){ // para cada objetivo m
				objetivoOrdenacao = m;
				I.sort(compare1);
				(*I.begin())->distance = INF;
				(I.back())->distance = INF; // É ISSO MESMO: back() NAO retorna um ponteiro interator!!! 
				it = I.begin();
				it++; // começa do 1 
				for (int i = 1; i<(l-1); i++){
					it++; // pos
					double objPos = (*it)->getObj(m);
					it--; // volta para o original
					it--; //previous
					double objPrev = (*it)->getObj(m);
					it++; // volta para o original
					(*it)->distance = (*it)->distance + (objPos - objPrev);
					it++; // avança
				}
			}
			#undef INF
		}

		// atualiza a populaçao (global)
		void atualizaPopulacaoNSGAII(Solution *novaPop[TAMANHOPOPULACAO_NSGAII]){
			for (int i=0; i<TAMANHOPOPULACAO_NSGAII*2; i++){
				if (i<TAMANHOPOPULACAO_NSGAII){
					*uniao[i] = *populacao[i];
				} else {
					*uniao[i] = *novaPop[i-TAMANHOPOPULACAO_NSGAII];
				}
			}
			int sizeFront = 0;
			list<Solution *> F[TAMANHOPOPULACAO_NSGAII*2];
			fast_non_dominanted_sort(uniao, F, sizeFront);
			int cont = 0;
			int i = 0;
			while (cont + F[i].size() < TAMANHOPOPULACAO_NSGAII && i<sizeFront){
				for (list< Solution* >::iterator p = F[i].begin(); p!=F[i].end(); p++){
					*populacao[cont++] = **p;
				}
				i++;
			}
			crownding_distance_assigment(F[i]);
			F[i].sort(compare2); // ordena por CD
			for (list< Solution* >::iterator p = F[i].begin(); p!=F[i].end() && cont<TAMANHOPOPULACAO_NSGAII; p++){
				*populacao[cont++] = **p;
			}
		}



	public:
		
		void iniPopulacao(){

			for (int i=0; i<TAMANHOPOPULACAO_NSGAII; i++){
				populacao[i]->generateRandomSolution(intContAval);
			}
		}

		nsgaii(){ 
		
			intContAval = 0;
			filho = new Solution(n_locations, n_objetivos); 
			for (int i=0; i<TAMANHOPOPULACAO_NSGAII; i++){
				populacao[i] = new Solution(n_locations, n_objetivos);
				novaPop[i] = new Solution(n_locations, n_objetivos);
			}
			for (int i=0; i<TAMANHOPOPULACAO_NSGAII*2; i++){
				uniao[i] = new Solution(n_locations, n_objetivos); // permanente
			}
			
		}

		Solution **executar(){ // retorna um array de ponteiros
			
			intContAval=0;
			iniPopulacao();
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

						int objjjj_ram = pick_a_number(0, n_objetivos-1);
						if(populacao[p1]->getObj(objjjj_ram)<populacao[p2]->getObj(objjjj_ram)){ // compete com o primeiro objetivo
							pai = p1;;
						} else {
							pai = p2;
						}

						objjjj_ram = pick_a_number(0, n_objetivos-1);
						if (populacao[p3]->getObj(objjjj_ram)<populacao[p4]->getObj(objjjj_ram)){ // compete com o primeiro objetivo
							mae = p3;
						} else {
							mae = p4;
						}
						
						filho->crossover(populacao[pai], populacao[mae],intContAval);
						// filho foi definido; Agora aplica-se mutaçao
						p = pick_a_number(0.0,1.0);
						if (p<TAXADEMUTACAO){
							novaPop[j]->mutation(filho,intContAval);
						} else *novaPop[j] = *filho;

					} else {
						novaPop[j]->generateRandomSolution(intContAval); // random
					}
				}
				atualizaPopulacaoNSGAII(novaPop);
			}
			return populacao;
		}
};



#endif