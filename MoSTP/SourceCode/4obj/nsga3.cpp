#ifndef NSGA3_CPP
#define NSGA3_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2019)
	
This code implements a NSGA-III algorithm for Multi-objective Spanning Tree Problem

NSGA-III functions were adapted from https://github.com/adanjoga/nsga3

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
#include "param_NSGAIII.h"
#include "PontoReferencia.cpp"

#include "rmcPrim.cpp"
#include "popInicial.cpp"


using namespace std;



int objetivoOrdenacao; // esta variavel é utilizada para designar qual objetivo será utilizado para ordenar as soluçoes
		

bool compare1(SolucaoEdgeSet *s1, SolucaoEdgeSet *s2){
	return s1->getObj(objetivoOrdenacao) < s2->getObj(objetivoOrdenacao);
}

// obsolete
bool compare2(SolucaoEdgeSet *s1, SolucaoEdgeSet *s2){
	return s1->distance > s2->distance; // maior distância
}


class nsga3 {
	private:
		SolucaoEdgeSet *uniao[TAMANHOPOPULACAO_NSGAII*2];
		SolucaoEdgeSet *populacao[TAMANHOPOPULACAO_NSGAII];
		SolucaoEdgeSet *novaPop[TAMANHOPOPULACAO_NSGAII]; // populaçao de descendentes gerada ao fim de uma geraçao
		SolucaoEdgeSet * filho;
		double *idealpoint;
		double *nadirpoint; 
		int N[TAMANHOPOPULACAO_NSGAII*2]; // para cada p \in P, guarda a quantidade de soluçoes que dominam p
		int intContAval;	
		vector<PontoReferencia> refepoints;	// reference points
		int p_boundary;
		// deb et al (2002)
		// P[TAMANHOPOPULACAO_NSGAII*2] = populacao corrente unido com a populaçao criada
		// tem que ser EXATAMENTE TAMANHOPOPULACAO_NSGAII*2
		// F é o vetor de fronts. Cada front (f0, f1, f2, ... ) é uma list. A quantidade maxima de fronts é TAMANHOPOPULACAO_NSGAII*2
		void fast_non_dominanted_sort(SolucaoEdgeSet *P[TAMANHOPOPULACAO_NSGAII*2], list<SolucaoEdgeSet *> F[TAMANHOPOPULACAO_NSGAII*2], int &quantFronts){
			list< SolucaoEdgeSet* > S[TAMANHOPOPULACAO_NSGAII*2]; // para cada p \in P, guarda a lista de soluçoes dominadas por p
			
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

		// void crownding_distance_assigment(list<SolucaoEdgeSet *> &I){
		// 	#define INF 1e9
		// 	int l = I.size();
		// 	list<SolucaoEdgeSet *>::iterator it = I.begin();
		// 	while (it!=I.end()){
		// 		(*it)->distance = 0;
		// 		it++;
		// 	}
		// 	for (int m=0; m<NUMOBJETIVOS; m++){ // para cada objetivo m
		// 		objetivoOrdenacao = m;
		// 		I.sort(compare1);
		// 		(*I.begin())->distance = INF;
		// 		(I.back())->distance = INF; // É ISSO MESMO: back() NAO retorna um ponteiro interator!!! 
		// 		it = I.begin();
		// 		it++; // começa do 1 
		// 		for (int i = 1; i<(l-1); i++){
		// 			it++; // pos
		// 			double objPos = (*it)->getObj(m);
		// 			it--; // volta para o original
		// 			it--; //previous
		// 			double objPrev = (*it)->getObj(m);
		// 			it++; // volta para o original
		// 			(*it)->distance = (*it)->distance + (objPos - objPrev);
		// 			it++; // avança
		// 		}
		// 	}
		// 	#undef INF
		// }


		/*
			=============================================================================
										Funcoes para o NSGA-III
			=============================================================================
		*/
		double distnciaPerpendicular(vector<double> &ponto_referencia, std::vector<double> &solution){
			double a = 0, b = 0;
			for (int i = 0; i < ponto_referencia.size(); i++) {
				a += ponto_referencia[i] * solution[i];
				b += ponto_referencia[i]*ponto_referencia[i];
			} 
			double a__b = (double) a/b;
			double distancia = 0;
			for (int i = 0; i < ponto_referencia.size(); i++) {
				distancia += (a__b * ponto_referencia[i] - solution[i])*(a__b * ponto_referencia[i] - solution[i]);
			}
			return sqrt(distancia);
		}

		size_t find_ponto_referencia(vector<PontoReferencia> &pontos_referencia, vector<bool> rp_isable){
			size_t min_size =  1e9;;
			for (size_t i = 0; i < pontos_referencia.size(); ++i) {
				if (rp_isable[i]){ // Se o ponto de referência tem membros que estao no ultimo front considerado
					if (pontos_referencia[i].getQuantMembros()<min_size)
						min_size = pontos_referencia[i].getQuantMembros();
				}
			}

			vector<size_t> min_pontos_referencia;
			for (size_t i = 0; i < pontos_referencia.size(); ++i)
				if (rp_isable[i]){ // Se o ponto de referência tem membros que estao no ultimo front considerado
					if(pontos_referencia[i].getQuantMembros() == min_size)
						min_pontos_referencia.push_back(i);
				}
			return min_pontos_referencia[IRandom(0, min_pontos_referencia.size()-1)];
		}

		SolucaoEdgeSet * seleciona_membro(PontoReferencia &ponto_referencia){
			if (ponto_referencia.tem_membros()) {
				if (ponto_referencia.getQuantMembros() == 0)
					return ponto_referencia.get_membro_mais_proximo(); // se getQuantMembros() == 0, como pode selecionar algum elemento aqui?
				else
					return ponto_referencia.get_membro_random();
			}
			return NULL;
		}



		void normalize(list<SolucaoEdgeSet *> Fronts[TAMANHOPOPULACAO_NSGAII*2], int quantFronts) { 
			// calcula o ponto ideal (basta consultar apenas o primeiro front, pois todos os outros fronts sao dominados pelo primeiro)
			for (std::list<SolucaoEdgeSet *>::iterator i = Fronts[0].begin(); i != Fronts[0].end(); ++i){
				for (size_t obj = 0; obj < NUMOBJETIVOS; ++obj) {
					if((*i)->getObj(obj)<idealpoint[obj]){ // O PROBLEMA EH DE MINIMIZACAO
						idealpoint[obj] = (*i)->getObj(obj);
					}
				}
			}


			//calcula o ponto Nadir
			for (size_t i = 0; i < quantFronts; ++i) {
				for (std::list<SolucaoEdgeSet *>::iterator itt = Fronts[i].begin(); itt != Fronts[i].end(); ++itt){
					for (size_t obj = 0; obj < NUMOBJETIVOS; ++obj) {
						if ((*itt)->getObj(obj)>nadirpoint[obj]){ // O PROBLEMA EH DE MINIMIZACAO
							nadirpoint[obj] = (*itt)->getObj(obj);
						}
					}
				}
			}

			
			// Normalize objectives
			for (size_t i = 0; i < quantFronts; ++i){
				for (std::list<SolucaoEdgeSet *>::iterator itt = Fronts[i].begin(); itt != Fronts[i].end(); ++itt){
					for (int obj = 0; obj < NUMOBJETIVOS; ++obj){
						if (fabs(nadirpoint[obj] - idealpoint[obj]) > 10e-10) //  nadirpoint[obj] - idealpoint[obj] diferente de 0
							(*itt)->getObj_NORMALIZADO[obj] = abs((*itt)->getObj(obj) - idealpoint[obj]) / abs(nadirpoint[obj] - idealpoint[obj]);
						else
							(*itt)->getObj_NORMALIZADO[obj] = abs((*itt)->getObj(obj) - idealpoint[obj]) / 10e-10;
					}
				}
			}
		}


		void associate(list<SolucaoEdgeSet *> Fronts[TAMANHOPOPULACAO_NSGAII*2], int quantFronts, int id_lastfront){
			for (size_t i = 0; i < refepoints.size(); ++i)
				refepoints[i].clear_ponto_Referencia();

			for (size_t i = 0; i < quantFronts; ++i) {
				for (std::list<SolucaoEdgeSet *>::iterator itt = Fronts[i].begin(); itt != Fronts[i].end(); ++itt){  // determina onde associar este ponto

					double mindist =  1e9;;
					size_t idx;
					for (size_t r = 0; r < refepoints.size(); ++r) {
						double pdist = distnciaPerpendicular(refepoints[r].ponto_Referencia,(*itt)->getObj_NORMALIZADO);
						if (pdist < mindist) {
							mindist = pdist;
							idx = r;
						}
					}
					if (i == id_lastfront)
						refepoints[idx].add_member((*itt), mindist);
					else
						refepoints[idx].sum_member();
				}
			}
		}

		void niching(int &conttt){ // prenche populacao
			vector<bool> rp_isable (refepoints.size(), true);
			while (conttt < TAMANHOPOPULACAO_NSGAII) {
				size_t rp_idx = find_ponto_referencia(refepoints, rp_isable);
				SolucaoEdgeSet *solution = seleciona_membro(refepoints[rp_idx]);
				
				if (solution==NULL) rp_isable[rp_idx] = false;
				else {
					refepoints[rp_idx].sum_member();
					refepoints[rp_idx].remove_member(solution);
					*populacao[conttt++] = *solution;
				}
			}
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
			int id_lastfront = 0; // id_lastfront
			// cout<<"sizeFront = "<<sizeFront<<endl;
			// cout<<"F[0].size() = "<<F[id_lastfront].size()<<endl;
			while (cont + F[id_lastfront].size() < TAMANHOPOPULACAO_NSGAII && id_lastfront<sizeFront){
				for (list< SolucaoEdgeSet* >::iterator p = F[id_lastfront].begin(); p!=F[id_lastfront].end(); p++){
					*populacao[cont++] = **p;
				}
				id_lastfront++;
			}

			if (cont == TAMANHOPOPULACAO_NSGAII) return;
			normalize(F, sizeFront);
			associate(F, sizeFront, id_lastfront);
			
			// for (int i = 0; i < refepoints.size(); ++i)
			// {
			// 	for (int j = 0; j < NUMOBJETIVOS; ++j)
			// 	{
			// 		cout<<refepoints[i].ponto_Referencia[j]<<" ";
			// 	}
			// 	cout<<endl;

			// 	for (int j = 0; j < refepoints[i].membros.size(); ++j)
			// 	{
			// 		cout<<"\t "<<refepoints[i].membros[j].first->getObj(0)<<" "<<refepoints[i].membros[j].first->getObj(1)<<" "<<refepoints[i].membros[j].first->getObj(2)<<endl;
			// 	}

			// }
			// exit(-1);
			// cout<<"Antes cont = "<<cont<<endl;
			niching(cont);
			// cout<<"Depois cont = "<<cont<<endl;


			// crownding_distance_assigment(F[i]);
			// F[i].sort(compare2); // ordena por CD
			// for (list< SolucaoEdgeSet* >::iterator p = F[i].begin(); p!=F[i].end() && cont<TAMANHOPOPULACAO_NSGAII; p++){
			// 	*populacao[cont++] = **p;
			// }
		}









		void gera_pontos_boundary(vector<PontoReferencia> *pontos, PontoReferencia *pt, size_t left, size_t count){
			if (count == NUMOBJETIVOS-1) {
				pt->ponto_Referencia[count] = 1.0 * left / p_boundary;
				pontos->push_back(*pt);
			}
			else {
				for (size_t i = 0; i <= left; ++i) {
					pt->ponto_Referencia[count] = 1.0 * i / p_boundary;
					gera_pontos_boundary(pontos, pt, left-i, count+1);
				}
			}
		}








	public:
		
		void editaPopulacao(){
			int init = 0;
			gerarPopulacao1(populacao, TAMANHOPOPULACAO_NSGAII, init, intContAval);
			
		}

		nsga3(){  // deve ser chamado em loadDecisionStrategies() (Silvia et al 2017)
			
			if (NUMOBJETIVOS==3){
				p_boundary = 23;
			} else if (NUMOBJETIVOS==4){
				p_boundary = 20;//10//11;
			}else if (NUMOBJETIVOS==5){
				p_boundary = 6;
			}

			PontoReferencia rpoint(NUMOBJETIVOS);
			gera_pontos_boundary(&refepoints, &rpoint, p_boundary, 0);
		

			// for (int i = 0; i < refepoints.size(); ++i)
			// {
			// 	for (int j = 0; j < NUMOBJETIVOS; ++j)
			// 	{
			// 		cout<<refepoints[i].ponto_Referencia[j]<<" ";
			// 	}
			// 	cout<<endl;
			// }



			idealpoint = new double[NUMOBJETIVOS];
			nadirpoint = new double[NUMOBJETIVOS];
			for (int i = 0; i < NUMOBJETIVOS; i++) {
				idealpoint[i] =  1e9;; // PROBLEMA DE MINIMIZACAO: o ponto que domina todos os pontos da fronteira
				nadirpoint[i] = INT_MIN; // PROBLEMA DE MINIMIZACAO: o ponto que eh dominado por todos os pontos da fronteira
			}


			intContAval=0;
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
			// for (int i=0; i<TAMANHOPOPULACAO_NSGAII; i++){
			// 	populacao[i]->printPonto(stdout);
			// }
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

						int ookn = IRandom(0, NUMOBJETIVOS-1);
						if(populacao[p1]->getObj(ookn)<populacao[p2]->getObj(ookn)){ // compete com o primeiro objetivo
							pai = p1;;
						} else {
							pai = p2;
						}

						ookn = IRandom(0, NUMOBJETIVOS-1);
						if (populacao[p3]->getObj(ookn)<populacao[p4]->getObj(ookn)){ // compete com o primeiro objetivo
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
			// for (int i=0; i<TAMANHOPOPULACAO_NSGAII; i++){
			// 	populacao[i]->printPonto(stdout);
			// }
			// cout<<endl;
			// cout<<endl;
			return populacao;
		}

};



#endif