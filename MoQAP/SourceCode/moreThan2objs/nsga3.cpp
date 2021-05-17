#ifndef NSGA3_CPP
#define NSGA3_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2019)
	
This code implements a NSGA-III algorithm for Multi-objective Quadratic Assignment Problem


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
#include "Solution.cpp"
#include "param_geral.h"
#include "PontoReferencia.cpp"


using namespace std;



int objetivoOrdenacao; // esta variavel é utilizada para designar qual objetivo será utilizado para ordenar as soluçoes
		

bool compare1(Solution *s1, Solution *s2){
	return s1->getObj(objetivoOrdenacao) < s2->getObj(objetivoOrdenacao);
}

bool compare2(Solution *s1, Solution *s2){
	return s1->distance > s2->distance; // maior distância
}


class nsga3 {
	private:
		Solution *uniao[TAMANHOPOPULACAO_NSGAII*2];
		Solution *populacao[TAMANHOPOPULACAO_NSGAII];
		Solution *novaPop[TAMANHOPOPULACAO_NSGAII]; // populaçao de descendentes gerada ao fim de uma geraçao
		Solution * filho;
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

		// void crownding_distance_assigment(list<Solution *> &I){
		// 	#define INF 1e9
		// 	int l = I.size();
		// 	list<Solution *>::iterator it = I.begin();
		// 	while (it!=I.end()){
		// 		(*it)->distance = 0;
		// 		it++;
		// 	}
		// 	for (int m=0; m<n_objetivos; m++){ // para cada objetivo m
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
			size_t min_size = INT_MAX;
			for (size_t i = 0; i < pontos_referencia.size(); ++i) {
				if (rp_isable[i]){ // Se o ponto de referência tem membros que estao no ultimo front considerado
					if (pontos_referencia[i].getQuantMembros()<min_size)
						min_size = pontos_referencia[i].getQuantMembros();
				}
			}

			vector<size_t> min_pontos_referencia;
			for (size_t i = 0; i < pontos_referencia.size(); ++i)
				if(pontos_referencia[i].getQuantMembros() == min_size)
					min_pontos_referencia.push_back(i);

			return min_pontos_referencia[pick_a_number(0, min_pontos_referencia.size()-1)];
		}

		Solution * seleciona_membro(PontoReferencia &ponto_referencia){
			if (ponto_referencia.tem_membros()) {
				if (ponto_referencia.getQuantMembros() == 0)
					return ponto_referencia.get_membro_mais_proximo(); // se getQuantMembros() == 0, como pode selecionar algum elemento aqui?
				else
					return ponto_referencia.get_membro_random();
			}
			return NULL;
		}



		void normalize(list<Solution *> Fronts[TAMANHOPOPULACAO_NSGAII*2], int quantFronts) { 
			// calcula o ponto ideal (basta consultar apenas o primeiro front, pois todos os outros fronts sao dominados pelo primeiro)
			for (std::list<Solution *>::iterator i = Fronts[0].begin(); i != Fronts[0].end(); ++i){
				for (size_t obj = 0; obj < n_objetivos; ++obj) {
					if((*i)->getObj(obj)<idealpoint[obj]){ // O PROBLEMA EH DE MINIMIZACAO
						idealpoint[obj] = (*i)->getObj(obj);
					}
				}
			}


			//calcula o ponto Nadir
			for (size_t i = 0; i < quantFronts; ++i) {
				for (std::list<Solution *>::iterator itt = Fronts[i].begin(); itt != Fronts[i].end(); ++itt){
					for (size_t obj = 0; obj < n_objetivos; ++obj) {
						if ((*itt)->getObj(obj)>nadirpoint[obj]){ // O PROBLEMA EH DE MINIMIZACAO
							nadirpoint[obj] = (*itt)->getObj(obj);
						}
					}
				}
			}

			
			// Normalize objectives
			for (size_t i = 0; i < quantFronts; ++i){
				for (std::list<Solution *>::iterator itt = Fronts[i].begin(); itt != Fronts[i].end(); ++itt){
					for (int obj = 0; obj < n_objetivos; ++obj){
						if (fabs(nadirpoint[obj] - idealpoint[obj]) > 10e-10) //  nadirpoint[obj] - idealpoint[obj] diferente de 0
							(*itt)->getObj_NORMALIZADO[obj] = ((*itt)->getObj(obj) - idealpoint[obj]) / (nadirpoint[obj] - idealpoint[obj]);
						else
							(*itt)->getObj_NORMALIZADO[obj] = ((*itt)->getObj(obj) - idealpoint[obj]) / 10e-10;
					}
				}
			}
		}


		void associate(list<Solution *> Fronts[TAMANHOPOPULACAO_NSGAII*2], int quantFronts, int id_lastfront){
			for (size_t i = 0; i < refepoints.size(); ++i)
				refepoints[i].clear_ponto_Referencia();

			for (size_t i = 0; i < quantFronts; ++i) {
				for (std::list<Solution *>::iterator itt = Fronts[i].begin(); itt != Fronts[i].end(); ++itt){  // determina onde associar este ponto

					double mindist = INT_MAX;
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

		void niching(int conttt){ // prenche populacao
			vector<bool> rp_isable (refepoints.size(), true);
			while (conttt < TAMANHOPOPULACAO_NSGAII) {
				size_t rp_idx = find_ponto_referencia(refepoints, rp_isable);
				Solution *solution = seleciona_membro(refepoints[rp_idx]);
				
				if (solution==NULL) rp_isable[rp_idx] = false;
				else {
					refepoints[rp_idx].sum_member();
					refepoints[rp_idx].remove_member(solution);
					*populacao[conttt++] = *solution;
				}
			}
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
			int id_lastfront = 0; // id_lastfront
			while (cont + F[id_lastfront].size() < TAMANHOPOPULACAO_NSGAII && id_lastfront<sizeFront){
				for (list< Solution* >::iterator p = F[id_lastfront].begin(); p!=F[id_lastfront].end(); p++){
					*populacao[cont++] = **p;
				}
				id_lastfront++;
			}

			if (cont == TAMANHOPOPULACAO_NSGAII) return;
			normalize(F, sizeFront);
			associate(F, sizeFront, id_lastfront);
			
			// for (int i = 0; i < refepoints.size(); ++i)
			// {
			// 	for (int j = 0; j < n_objetivos; ++j)
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

			niching(cont);
			


			// crownding_distance_assigment(F[i]);
			// F[i].sort(compare2); // ordena por CD
			// for (list< Solution* >::iterator p = F[i].begin(); p!=F[i].end() && cont<TAMANHOPOPULACAO_NSGAII; p++){
			// 	*populacao[cont++] = **p;
			// }
		}









		void gera_pontos_boundary(vector<PontoReferencia> *pontos, PontoReferencia *pt, size_t left, size_t count){
			if (count == n_objetivos-1) {
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
		
		void iniPopulacao(){

			for (int i=0; i<TAMANHOPOPULACAO_NSGAII; i++){
				populacao[i]->generateRandomSolution(intContAval);
			}
		}

		nsga3(){ 

			if (n_objetivos==3){
				p_boundary = 23;
			} else if (n_objetivos==4){
				p_boundary = 10;
			}else if (n_objetivos==5){
				p_boundary = 6;
			}

			PontoReferencia rpoint(n_objetivos);
			gera_pontos_boundary(&refepoints, &rpoint, p_boundary, 0);


			// for (int i = 0; i < refepoints.size(); ++i)
			// {
			// 	for (int j = 0; j < n_objetivos; ++j)
			// 	{
			// 		cout<<refepoints[i].ponto_Referencia[j]<<" ";
			// 	}
			// 	cout<<endl;
			// }


			idealpoint = new double[n_objetivos];
			nadirpoint = new double[n_objetivos];
			for (int i = 0; i < n_objetivos; i++) {
				idealpoint[i] = INT_MAX; // PROBLEMA DE MINIMIZACAO: o ponto que domina todos os pontos da fronteira
				nadirpoint[i] = INT_MIN; // PROBLEMA DE MINIMIZACAO: o ponto que eh dominado por todos os pontos da fronteira
			}

		
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