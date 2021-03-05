#ifndef NSGA3_CPP
#define NSGA3_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2019)
	
This code implements a NSGA-III algorithm for Multiobjective Multidimensional Knapsack problem

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
#include "Knapsack.cpp"
#include "param_geral.h"
#include "PontoReferencia.cpp"


using namespace std;



int objetivoOrdenacao; // esta variavel é utilizada para designar qual objetivo será utilizado para ordenar as soluçoes
		

bool compare1(Knapsack *s1, Knapsack *s2){
	return s1->profit(objetivoOrdenacao) < s2->profit(objetivoOrdenacao);
}

bool compare2(Knapsack *s1, Knapsack *s2){
	return s1->distance > s2->distance; // maior distância
}


class nsga3 {
	private:
		Knapsack *uniao[TAMANHOPOPULACAO_NSGAII*2];
		Knapsack *populacao[TAMANHOPOPULACAO_NSGAII];
		Knapsack *novaPop[TAMANHOPOPULACAO_NSGAII]; // populaçao de descendentes gerada ao fim de uma geraçao
		Knapsack * filho;
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

		// void crownding_distance_assigment(list<Knapsack *> &I){
		// 	#define INF 1e9
		// 	int l = I.size();
		// 	list<Knapsack *>::iterator it = I.begin();
		// 	while (it!=I.end()){
		// 		(*it)->distance = 0;
		// 		it++;
		// 	}
		// 	for (int m=0; m<my_momkp->dimension_; m++){ // para cada objetivo m
		// 		objetivoOrdenacao = m;
		// 		I.sort(compare1);
		// 		(*I.begin())->distance = INF;
		// 		(I.back())->distance = INF; // É ISSO MESMO: back() NAO retorna um ponteiro interator!!! 
		// 		it = I.begin();
		// 		it++; // começa do 1 
		// 		for (int i = 1; i<(l-1); i++){
		// 			it++; // pos
		// 			double objPos = (*it)->profit(m);
		// 			it--; // volta para o original
		// 			it--; //previous
		// 			double objPrev = (*it)->profit(m);
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

		Knapsack * seleciona_membro(PontoReferencia &ponto_referencia){
			if (ponto_referencia.tem_membros()) {
				if (ponto_referencia.getQuantMembros() == 0)
					return ponto_referencia.get_membro_mais_proximo(); // se getQuantMembros() == 0, como pode selecionar algum elemento aqui?
				else
					return ponto_referencia.get_membro_random();
			}
			return NULL;
		}



		void normalize(list<Knapsack *> Fronts[TAMANHOPOPULACAO_NSGAII*2], int quantFronts) { 
			// calcula o ponto ideal (basta consultar apenas o primeiro front, pois todos os outros fronts sao dominados pelo primeiro)
			for (std::list<Knapsack *>::iterator i = Fronts[0].begin(); i != Fronts[0].end(); ++i){
				for (size_t obj = 0; obj < my_momkp->dimension_; ++obj) {
					if((*i)->profit(obj)>idealpoint[obj]){ // O PROBLEMA EH DE MAXIMIZACAO
						idealpoint[obj] = (*i)->profit(obj);
					}
				}
			}


			//calcula o ponto Nadir
			for (size_t i = 0; i < quantFronts; ++i) {
				for (std::list<Knapsack *>::iterator itt = Fronts[i].begin(); itt != Fronts[i].end(); ++itt){
					for (size_t obj = 0; obj < my_momkp->dimension_; ++obj) {
						if ((*itt)->profit(obj)<nadirpoint[obj]){ // O PROBLEMA EH DE MAXIMIZACAO
							nadirpoint[obj] = (*itt)->profit(obj);
						}
					}
				}
			}

			
			// Normalize objectives
			for (size_t i = 0; i < quantFronts; ++i){
				for (std::list<Knapsack *>::iterator itt = Fronts[i].begin(); itt != Fronts[i].end(); ++itt){
					for (int obj = 0; obj < my_momkp->dimension_; ++obj){
						if (fabs(nadirpoint[obj] - idealpoint[obj]) > 10e-10) //  nadirpoint[obj] - idealpoint[obj] diferente de 0
							(*itt)->profit_NORMALIZADO[obj] = ((*itt)->profit(obj) - idealpoint[obj]) / (nadirpoint[obj] - idealpoint[obj]);
						else
							(*itt)->profit_NORMALIZADO[obj] = ((*itt)->profit(obj) - idealpoint[obj]) / 10e-10;
					}
				}
			}
		}


		void associate(list<Knapsack *> Fronts[TAMANHOPOPULACAO_NSGAII*2], int quantFronts, int id_lastfront){
			for (size_t i = 0; i < refepoints.size(); ++i)
				refepoints[i].clear_ponto_Referencia();

			for (size_t i = 0; i < quantFronts; ++i) {
				for (std::list<Knapsack *>::iterator itt = Fronts[i].begin(); itt != Fronts[i].end(); ++itt){  // determina onde associar este ponto

					double mindist = INT_MAX; 
					size_t idx;
					for (size_t r = 0; r < refepoints.size(); ++r) {
						double pdist = distnciaPerpendicular(refepoints[r].ponto_Referencia,(*itt)->profit_NORMALIZADO);
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
				Knapsack *solution = seleciona_membro(refepoints[rp_idx]);
				
				if (solution==NULL) rp_isable[rp_idx] = false;
				else {
					refepoints[rp_idx].sum_member();
					refepoints[rp_idx].remove_member(solution);
					*populacao[conttt++] = *solution;
				}
			}
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
			int id_lastfront = 0; // id_lastfront
			while (cont + F[id_lastfront].size() < TAMANHOPOPULACAO_NSGAII && id_lastfront<sizeFront){
				for (list< Knapsack* >::iterator p = F[id_lastfront].begin(); p!=F[id_lastfront].end(); p++){
					*populacao[cont++] = **p;
				}
				id_lastfront++;
			}

			if (cont == TAMANHOPOPULACAO_NSGAII) return;
			normalize(F, sizeFront);
			associate(F, sizeFront, id_lastfront);
			
			// for (int i = 0; i < refepoints.size(); ++i)
			// {
			// 	for (int j = 0; j < my_momkp->dimension_; ++j)
			// 	{
			// 		cout<<refepoints[i].ponto_Referencia[j]<<" ";
			// 	}
			// 	cout<<endl;

			// 	for (int j = 0; j < refepoints[i].membros.size(); ++j)
			// 	{
			// 		cout<<"\t "<<refepoints[i].membros[j].first->profit(0)<<" "<<refepoints[i].membros[j].first->profit(1)<<" "<<refepoints[i].membros[j].first->profit(2)<<endl;
			// 	}

			// }
			// exit(-1);

			niching(cont);
			


			// crownding_distance_assigment(F[i]);
			// F[i].sort(compare2); // ordena por CD
			// for (list< Knapsack* >::iterator p = F[i].begin(); p!=F[i].end() && cont<TAMANHOPOPULACAO_NSGAII; p++){
			// 	*populacao[cont++] = **p;
			// }
		}









		void gera_pontos_boundary(vector<PontoReferencia> *pontos, PontoReferencia *pt, size_t left, size_t count){
			if (count == my_momkp->dimension_-1) {
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


			std::vector<Knapsack *> population_greedy;
			my_momkp->generate_population(population_greedy, TAMANHOPOPULACAO_NSGAII,0, 0.1, intContAval); // 1-Random 0-Greedy
			int ontt=0;
			for (std::vector<Knapsack *>::iterator i = population_greedy.begin(); i != population_greedy.end() && ontt<TAMANHOPOPULACAO_NSGAII; ++i){
				*populacao[ontt++] = **i;
			}
		}

		nsga3(){ 

			if (my_momkp->dimension_==3){
				p_boundary = 23;
			} else if (my_momkp->dimension_==4){
				p_boundary = 10;
			}else if (my_momkp->dimension_==5){
				p_boundary = 6;
			}

			PontoReferencia rpoint(my_momkp->dimension_);
			gera_pontos_boundary(&refepoints, &rpoint, p_boundary, 0);


			// for (int i = 0; i < refepoints.size(); ++i)
			// {
			// 	for (int j = 0; j < my_momkp->dimension_; ++j)
			// 	{
			// 		cout<<refepoints[i].ponto_Referencia[j]<<" ";
			// 	}
			// 	cout<<endl;
			// }


			idealpoint = new double[my_momkp->dimension_];
			nadirpoint = new double[my_momkp->dimension_];
			for (int i = 0; i < my_momkp->dimension_; i++) {
				idealpoint[i] = INT_MIN; // PROBLEMA DE MAXIMIZACAO: o ponto que domina todos os pontos da fronteira
				nadirpoint[i] = INT_MAX; // PROBLEMA DE MAXIMIZACAO: o ponto que eh dominado por todos os pontos da fronteira
			}

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
				//cout<<"geracao = "<<geracao<<" intContAval = "<<intContAval<<endl;
			
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