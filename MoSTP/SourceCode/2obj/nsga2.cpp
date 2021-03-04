#ifndef NSGA2_CPP
#define NSGA2_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2018)
	
This code implements a NSGAII algorithm for Bi-objective spanning tree 
The data structure and some functions were kindly provided by Monteiro (2010)


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
#include "SolucaoEdgeSet.cpp"
#include "rmcPrim.cpp"
#include "popInicial.cpp"
#include "param_NSGAII.h"


using namespace std;


int objetivoOrdenacao; // esta variavel é utilizada para designar qual objetivo será utilizado para ordenar as soluçoes
		

bool compare1(SolucaoEdgeSet *s1, SolucaoEdgeSet *s2){
	return s1->getObj(objetivoOrdenacao) < s2->getObj(objetivoOrdenacao);
}

bool compare2(SolucaoEdgeSet *s1, SolucaoEdgeSet *s2){
	return s1->distance > s2->distance; // maior distância
}


class NSGA2 {
	private:
		SolucaoEdgeSet *uniao[TAMANHOPOPULACAO_NSGAII*2];
		SolucaoEdgeSet *populacao[TAMANHOPOPULACAO_NSGAII];
		SolucaoEdgeSet *novaPop[TAMANHOPOPULACAO_NSGAII]; // populaçao de descendentes gerada ao fim de uma geraçao
		SolucaoEdgeSet * filho;
		int N[TAMANHOPOPULACAO_NSGAII*2]; // para cada p \in P, guarda a quantidade de soluçoes que dominam p
		int intContAval;	

		// deb et al (2002)
		// P[TAMANHOPOPULACAO_NSGAII*2] = populacao corrente unido com a populaçao criada
		// tem que ser EXATAMENTE TAMANHOPOPULACAO_NSGAII*2
		// F é o vetor de fronts. Cada front (f0, f1, f2, ... ) é uma list. A quantidade maxima de fronts é TAMANHOPOPULACAO_NSGAII*2
		void fast_non_dominanted_sort(SolucaoEdgeSet *P[TAMANHOPOPULACAO_NSGAII*2], list<SolucaoEdgeSet *> F[TAMANHOPOPULACAO_NSGAII*2], int &quantFronts){
			list< SolucaoEdgeSet* > S[TAMANHOPOPULACAO_NSGAII*2]; // para cada p \in P, guarda a lista de soluçoes dominadas por p
			
			// int prank[TAMANHOPOPULACAO_NSGAII*2]; // para cada p \in P, guarda o rank de p
			quantFronts = 0;
			for (int p=0; p<TAMANHOPOPULACAO_NSGAII*2; p++){
				N[p] = 0;
				// prank[p] = 0;
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
					// prank[p] = 0; // primeiro rank. RANK COMEÇA DO 0 
					P[p]->rank = 0;
					F[0].push_back(P[p]); // primeiro rank. RANK COMEÇA DO 0
				}
			}
			int i=0; 
			while (F[i].size()>0){
				for (list< SolucaoEdgeSet* >::iterator p = F[i].begin(); p!=F[i].end(); p++){
					int pp = (*p)->posicaoListaNSGAII;
					for (list< SolucaoEdgeSet* >::iterator q = S[pp].begin(); q!=S[pp].end(); q++){
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

		// // verifica se a lista de Fronts está correta: apenas para fins de verificaçao da corretude
		// bool verificador(list<SolucaoEdgeSet *> F[TAMANHOPOPULACAO_NSGAII*2], int sizeFront){
		// 	for (int i=0; i<sizeFront; i++){
		// 		for (list< SolucaoEdgeSet* >::iterator p = F[i].begin(); p!=F[i].end(); p++){
		// 			for (list< SolucaoEdgeSet* >::iterator pp = F[i].begin(); pp!=F[i].end(); pp++){
		// 				if (p!=pp){
		// 					if (**p>>**pp || **pp>>**p) {
		// 						cout<<"RROOOOOOOOOOOOOOOOOOOOOOOOOO 1"<<endl;
		// 						return false;
		// 					}
		// 				}
		// 			}
		// 			bool encontrou = false;
		// 			for (int ii = i-1; ii>=0 && encontrou==false; ii--){
		// 				for (list< SolucaoEdgeSet* >::iterator pp = F[ii].begin(); pp!=F[ii].end() && encontrou==false; pp++){
		// 					if (**pp>>**p) encontrou = true;
		// 				}
		// 			}	
		// 			if (i!=0 && encontrou==false){
		// 				cout<<"RROOOOOOOOOOOOOOOOOOOOOOOOOO 2 "<<endl;//i = "<<i<<" ==> "<<(*p)->getObj(0)<<" "<<(*p)->getObj(1)<<endl;
		// 				return false;
		// 			} 

		// 		}
		// 	}
		// 	return true;
		// }

		void crownding_distance_assigment(list<SolucaoEdgeSet *> &I){
			#define INF 1e9
			int l = I.size();
			list<SolucaoEdgeSet *>::iterator it = I.begin();
			while (it!=I.end()){
				(*it)->distance = 0;
				it++;
			}
			for (int m=0; m<NUMOBJETIVOS; m++){ // para cada objetivo m
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
		void atualizaPopulacaoNSGAII(SolucaoEdgeSet *novaPop[TAMANHOPOPULACAO_NSGAII]){
			for (int i=0; i<TAMANHOPOPULACAO_NSGAII*2; i++){
				if (i<TAMANHOPOPULACAO_NSGAII){
					*uniao[i] = *populacao[i];
				} else {
					*uniao[i] = *novaPop[i-TAMANHOPOPULACAO_NSGAII];
				}
			}
			int sizeFront = 0;
			list<SolucaoEdgeSet *> F[TAMANHOPOPULACAO_NSGAII*2];
			fast_non_dominanted_sort(uniao, F, sizeFront);
			int cont = 0;
			int i = 0;
			while (cont + F[i].size() < TAMANHOPOPULACAO_NSGAII && i<sizeFront){
				for (list< SolucaoEdgeSet* >::iterator p = F[i].begin(); p!=F[i].end(); p++){
					*populacao[cont++] = **p;
				}
				i++;
			}
			crownding_distance_assigment(F[i]);
			F[i].sort(compare2); // ordena por CD
			for (list< SolucaoEdgeSet* >::iterator p = F[i].begin(); p!=F[i].end() && cont<TAMANHOPOPULACAO_NSGAII; p++){
				*populacao[cont++] = **p;
			}
		}

	public:
		
		void editaPopulacao(){
			int init = 0;
			
			gerarPopulacao1(populacao, TAMANHOPOPULACAO_NSGAII, init, intContAval);
			
		}

		NSGA2(){  // deve ser chamado em loadDecisionStrategies() (Silvia et al 2017)
			
		
			filho = new SolucaoEdgeSet(NUMEROVERTICES-1); 
			for (int i=0; i<TAMANHOPOPULACAO_NSGAII; i++){
				populacao[i] = new SolucaoEdgeSet(NUMEROVERTICES-1);
				novaPop[i] = new SolucaoEdgeSet(NUMEROVERTICES-1);
			}
			for (int i=0; i<TAMANHOPOPULACAO_NSGAII*2; i++){
				uniao[i] = new SolucaoEdgeSet(NUMEROVERTICES-1); // permanente
			}
			
		}

		SolucaoEdgeSet **executar(){ // retorna um array de ponteiros
			
			intContAval=0;
			editaPopulacao();
			int p1,p2,p3,p4;
			int pai, mae;	
		
			int geracao = 0;
			while (intContAval<NUM_AVALIACOES){
				// if(geracao%100 == 0) {
				// 	cout<<"geracao NSGA2 = "<<geracao<<" avaliacoes = "<<intContAval<<endl; //DEVE SAIR
				// }
				geracao++;
				for (int j=0; j<TAMANHOPOPULACAO_NSGAII; j++){ // deve-se criar TAMANHOPOPULACAO_NSGAII novos individuos
					double p = genrand64_real3();;//rg.Random();
					if (p<TAXADECRUZAMENTO){
						/*SORTEIA 4 individuos*/
						/*Faz-se o torneio binario entre eles*/
						p1 = IRandom(0, TAMANHOPOPULACAO_NSGAII-1);
						p2 = IRandom(0, TAMANHOPOPULACAO_NSGAII-1);
						p3 = IRandom(0, TAMANHOPOPULACAO_NSGAII-1);
						p4 = IRandom(0, TAMANHOPOPULACAO_NSGAII-1);

						if(populacao[p1]->getObj(0)<populacao[p2]->getObj(0)){ // compete com o primeiro objetivo
							pai = p1;;
						} else {
							pai = p2;
						}
						if (populacao[p3]->getObj(1)<populacao[p4]->getObj(1)){ // compete com o primeiro objetivo
							mae = p3;
						} else {
							mae = p4;
						}
						
						filho->crossover(populacao[pai], populacao[mae],intContAval);
						// filho foi definido; Agora aplica-se mutaçao
						p = genrand64_real3();;//rg.Random();
						if (p<TAXADEMUTACAO){
							novaPop[j]->mutacao(*filho,intContAval);
						} else{
							*novaPop[j] = *filho;
						}

					} else {
						filho->doRandomWalk(intContAval);
						*novaPop[j] = *filho;
					}
				}
				atualizaPopulacaoNSGAII(novaPop);
			}

			return populacao;
		}
};



#endif