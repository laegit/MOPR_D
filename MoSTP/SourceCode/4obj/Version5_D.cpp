#ifndef VERSION5_D_CPP
#define VERSION5_D_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2019)

Baseado nos papers 
	"A New Evolutionary Algorithm Based on Decomposition for Multi-objective Optimization Problems" (Dai e Lei (2016)).
	"A new multi-objective particle swarm optimization algorithm based on decomposition"

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

typedef struct vviz
{
	int subproblema;
	double distance;
} Vizinho;



class Version5_D: public PathRelinking{
	private:
		int amostralSai[NUMEROVERTICES-1]; // guarda indexes de arestas (de s1) que devem sair de s1 (vide abaixo)
		int amostralEnt[NUMEROVERTICES-1]; // guarda indexes de arestas (de s2) que entrar em s1 (vide abaixo)

		SolucaoEdgeSet *aux, *s2;
		double lambda[NUMSUMPROBLEMAS][NUMOBJETIVOS]; // o lambda para cada subproblema
		list<SolucaoEdgeSet * > subespacos[NUMSUMPROBLEMAS];
		double ideal[NUMSUMPROBLEMAS][NUMOBJETIVOS]; // um ponto ideal para cada subproblema
		double Z_min[NUMOBJETIVOS];
		double Z_max[NUMOBJETIVOS];
		double Z_min_sub[NUMSUMPROBLEMAS][NUMOBJETIVOS];
		Vizinho B_vizinhos[NUMSUMPROBLEMAS][NUM_VIZINHOS]; // para cada subproblema (lambda), temos os seus vizinhos
		BoundedParetoSet *retorno;
		int subproblemaATUAL;

		double calcula_t(BoundedParetoSet *approx_set){
			double min_t = INT_MAX;
			for ( list<SolucaoEdgeSet *>::iterator it_1=approx_set->sol.begin(); it_1!=approx_set->sol.end(); it_1++){
				// calcula-se a norma de it_1 - Z_min (distância euclidiana)
				double norma_t = 0.0;
				for (int i=0; i<NUMOBJETIVOS; i++){
					norma_t += (Z_min[i]-(*it_1)->getObj(i))*(Z_min[i]-(*it_1)->getObj(i));
				}
				double de_t = sqrt(norma_t);
				if (de_t<min_t) min_t = de_t;
			}
			return min_t;
		}

		void calculaZ_por_subespaco(int sub){
			Z_min_sub[sub][0] = INT_MAX;
			Z_min_sub[sub][1] = INT_MAX;
			Z_min_sub[sub][2] = INT_MAX;
			Z_min_sub[sub][3] = INT_MAX;
			for (list<SolucaoEdgeSet *>::iterator i = subespacos[sub].begin(); i!=subespacos[sub].end(); i++){
				if ((*i)->getObj(0)<Z_min_sub[sub][0]) Z_min_sub[sub][0] = (*i)->getObj(0);
				if ((*i)->getObj(1)<Z_min_sub[sub][1]) Z_min_sub[sub][1] = (*i)->getObj(1); 
				if ((*i)->getObj(2)<Z_min_sub[sub][2]) Z_min_sub[sub][2] = (*i)->getObj(2); 
				if ((*i)->getObj(3)<Z_min_sub[sub][3]) Z_min_sub[sub][3] = (*i)->getObj(3); 
			
			}
		}

		void calculaZ_inicial(BoundedParetoSet *approx_set){
			Z_min[0] = INT_MAX;
			Z_min[1] = INT_MAX;
			Z_min[2] = INT_MAX;
			Z_min[3] = INT_MAX;

			Z_max[0] = INT_MIN;
			Z_max[1] = INT_MIN;
			Z_max[2] = INT_MIN;
			Z_max[3] = INT_MIN;
			for ( list<SolucaoEdgeSet *>::iterator it_1=approx_set->sol.begin(); it_1!=approx_set->sol.end(); it_1++){
				if ((*it_1)->getObj(0)<Z_min[0]) Z_min[0] = (*it_1)->getObj(0);
				if ((*it_1)->getObj(1)<Z_min[1]) Z_min[1] = (*it_1)->getObj(1); 
				if ((*it_1)->getObj(2)<Z_min[2]) Z_min[2] = (*it_1)->getObj(2); 
				if ((*it_1)->getObj(3)<Z_min[3]) Z_min[3] = (*it_1)->getObj(3); 


				if ((*it_1)->getObj(0)>Z_max[0]) Z_max[0] = (*it_1)->getObj(0);
				if ((*it_1)->getObj(1)>Z_max[1]) Z_max[1] = (*it_1)->getObj(1);
				if ((*it_1)->getObj(2)>Z_max[2]) Z_max[2] = (*it_1)->getObj(2);
				if ((*it_1)->getObj(3)>Z_max[3]) Z_max[3] = (*it_1)->getObj(3);
				
			}
		}	



		// pre-condiçao : 0 <= subpro < NUMSUMPROBLEMAS
	 	double delta(SolucaoEdgeSet *fx, int subpro){
			double norma_fx = 0.0;
			double norma_lambda = 0.0;
			double produto_escalar = 0.0;

			for (int i=0; i<NUMOBJETIVOS; i++){
				norma_lambda += lambda[subpro][i]*lambda[subpro][i];
				norma_fx += ( fx->getObj(i)-Z_min[i])*( fx->getObj(i)-Z_min[i]); 
				produto_escalar += (lambda[subpro][i]*(fx->getObj(i)-Z_min[i]));// *(lambda[subpro][i]*(objj_normalizado-0.0));
			}

			double norma_lambda_aux = sqrt(norma_lambda);
			double norma_fx_aux = sqrt(norma_fx);
			return (double) produto_escalar/(norma_fx_aux*norma_lambda_aux);
		}



		// insere nova_solu se ela nao for dominada. As dominadas saem. Se atingir o tamanho maximo da subpopulaçao, retira-se a de maior distância para o ponto ideal
		bool add_solution_subespaco(SolucaoEdgeSet *nova_solu, int subb){
			
			/* percorre o vetor de solucoes e de valores e, caso exista solucao dominada, retira e retorna true. caso contrario, retorna false */
			list<SolucaoEdgeSet *>::iterator i = subespacos[subb].begin();
			list< list<SolucaoEdgeSet *>::iterator > remover;
			list<SolucaoEdgeSet *>::iterator maisPopuloso;
			double max_distance_i = INT_MIN;
			while (i != subespacos[subb].end()) {
				// se a solucao que vai entrar domina a solucao i 
				if (*nova_solu >> **i) {
					remover.push_back(i); 
				}
				
				double distance_it = sqrt(((*i)->getObj(0)-ideal[subb][0])*((*i)->getObj(0)-ideal[subb][0]) + ((*i)->getObj(1)-ideal[subb][1])*((*i)->getObj(1)-ideal[subb][1])  +  ((*i)->getObj(2)-ideal[subb][2])*((*i)->getObj(2)-ideal[subb][2])   +  ((*i)->getObj(3)-ideal[subb][3])*((*i)->getObj(3)-ideal[subb][3]) );
				if (distance_it>max_distance_i){ // get a soluçao de maior distância para o ideal[subb]
					max_distance_i = distance_it;
					maisPopuloso = i;
				}
				
				if (**i >> *nova_solu || **i == *nova_solu){ // caso algum "i" domine "s", para e retorna false
					return false;
				}
				i++;
			}

			// se a solucao que vai entrar nao domina nenhuma e o tamanho da populaçao ja estah no maximo
			// (se nenhuma solucao vai sair do conjunto), remove a de maior distância para o ponto ideal
			if (remover.size() == 0 && subespacos[subb].size()+1 > MAX_SUB_POPULATION) {
				
				double distance_nova_solu = sqrt((nova_solu->getObj(0)-ideal[subb][0])*(nova_solu->getObj(0)-ideal[subb][0]) + (nova_solu->getObj(1)-ideal[subb][1])*(nova_solu->getObj(1)-ideal[subb][1]) +  (nova_solu->getObj(2)-ideal[subb][2])*(nova_solu->getObj(2)-ideal[subb][2])  +   (nova_solu->getObj(3)-ideal[subb][3])*(nova_solu->getObj(3)-ideal[subb][3]));
				if (distance_nova_solu>max_distance_i) return false;

				remover.push_back(maisPopuloso);
			}

			list< list<SolucaoEdgeSet *>::iterator >::iterator j = remover.begin();
			while (j != remover.end()) {
				delete( **j ); 
				// remove do conjunto pareto
				subespacos[subb].erase( *j );
				j++;
			}

			SolucaoEdgeSet *t = new SolucaoEdgeSet(NUMEROVERTICES-1);
			*t = *nova_solu; /// CRIA COPIA
			// adiciona ao conjunto pareto
			subespacos[subb].push_front( t );
			if (subespacos[subb].size() > MAX_SUB_POPULATION) {
				fprintf(stderr,"ERRO in MAX_SUB_POPULATION!\n");
				exit(-1);
			}
			return true;		
		}



		void initialize(BoundedParetoSet *approx_set){
			calculaZ_inicial(approx_set);
			// divide os subespaços
			for ( list<SolucaoEdgeSet *>::iterator it_1=approx_set->sol.begin(); it_1!=approx_set->sol.end(); it_1++){
				retorno->adicionarSol(*it_1); // copia, previamente, o conjunto inicial.
			
				//calcula os sub-espaços objetivo para cada subproblema
				double cosine_max = INT_MIN;
				int subpro = -1;
				for (int i=0; i<NUMSUMPROBLEMAS; i++){
					
					double cosine = delta((*it_1), i);
					if (cosine>cosine_max){
						cosine_max = cosine;
						subpro = i;
					}
				}
				if (subpro!=-1){
					// subespacos[subpro].push_back(copia);
					add_solution_subespaco(*it_1, subpro);
				}
			}


			// se algum sub-espaço estiver vazio, entao inserimos nele a soluçao de menor ângulo 
			// QUESTAO: E se fizéssimos assim para todos os sub-espaços?????
			for (int i=0; i<NUMSUMPROBLEMAS; i++){

				if (subespacos[i].size()==0){
					double cosine_max = INT_MIN;
					list<SolucaoEdgeSet *>::iterator sol_mais_prox;
					for ( list<SolucaoEdgeSet *>::iterator it_1=approx_set->sol.begin(); it_1!=approx_set->sol.end(); it_1++){
						double cosine = delta((*it_1), i);
						if (cosine>cosine_max){
							cosine_max = cosine;
							sol_mais_prox = it_1;
						}
					}
					add_solution_subespaco(*sol_mais_prox, i);
				}
			}
			

			//inicializa o ponto ideal (onde Z_min é calculado dentro de cada sub-espaço)
			for (int i = 0; i<NUMSUMPROBLEMAS; i++){
				calculaZ_por_subespaco(i);
				for (int obj = 0; obj<NUMOBJETIVOS; obj++){
					ideal[i][obj] = Z_min_sub[i][obj]; 
				}
			}

			

			// calcula os NUM_VIZINHOS mais proximos de cada vetor i
			for (int i=0; i<NUMSUMPROBLEMAS; i++){
				int contVizinhos = 0;
				for (int j=0; j<NUMSUMPROBLEMAS; j++){
					if (j!=i){
						double sumDistance =  0;
						for (int obj=0; obj<NUMOBJETIVOS; obj++){
							sumDistance+=(lambda[i][obj]-lambda[j][obj])*(lambda[i][obj]-lambda[j][obj]);
						}	
						double distanceEuclidiana = sqrt(sumDistance);
						
						int inserted = 0;
						if (contVizinhos==NUM_VIZINHOS){
							if (distanceEuclidiana<B_vizinhos[i][contVizinhos-1].distance){
								B_vizinhos[i][contVizinhos-1].distance = distanceEuclidiana;
								B_vizinhos[i][contVizinhos-1].subproblema = j;
								inserted = contVizinhos-1;
							}
						} else {
							B_vizinhos[i][contVizinhos].distance = distanceEuclidiana;
							B_vizinhos[i][contVizinhos].subproblema = j;
							inserted = contVizinhos;
							contVizinhos++;
						}
						while (inserted>0 && B_vizinhos[i][inserted].distance<B_vizinhos[i][inserted-1].distance){
							double aux = B_vizinhos[i][inserted].distance;
							B_vizinhos[i][inserted].distance = B_vizinhos[i][inserted-1].distance;
							B_vizinhos[i][inserted-1].distance = aux;
								
							int aux2 = B_vizinhos[i][inserted].subproblema;
							B_vizinhos[i][inserted].subproblema = B_vizinhos[i][inserted-1].subproblema;
							B_vizinhos[i][inserted-1].subproblema = aux2;

							inserted--;
						}
					}
				}
			}
		}


		void PR(BoundedParetoSet *retorno, SolucaoEdgeSet *s11, SolucaoEdgeSet *s222, int &contAval){
		
			int contttttt = 0;
			*aux = *s11; // O(n)
			*s2 = *s222; // O(n)

			// O(n^2)
			int totalSai = inS1_NotInS2(aux, s2, amostralSai); // arestas de s1 que nao estao em s2. Tais arestas devem sair de s1
			int totalEnt = inS1_NotInS2(s2, aux, amostralEnt); // arestas de s2 que nao estao em s1. Tais arestas devem entrar em s1

			// retorno->adicionarSol(aux);

			for (int cont = 0; cont<totalSai && contAval<NUM_AVALIACOES; cont++){
				
				int sai = -1;
				int entra = -1;
				// int subbespaco = -1;
				double min_objetivo = INT_MAX; // minimiza a maior distância
				
				// para cada par possivel (sai, entra)		
				for (int iii=0; iii<totalSai; iii++){ // para cada aresta i que deve sair de s1 O(n)

					if (amostralSai[iii]!=-1){

						int i = amostralSai[iii]; // aresta i que sai
						int origem = aux->edges[i][0];
						int destino= aux->edges[i][1];

						// controi o UnionFind O(nlogn)
						aux->uf.clear();
						for (int j=0;j<NUMEROVERTICES-1;j++){ // para toda aresta j diferente de i (contabilizar o conjunto uf) 
							if (i!=j) 
								aux->uf.unionClass(aux->edges[j][0], aux->edges[j][1]); 
						}			

						for (int jd=0; jd<totalEnt; jd++){
							int j = amostralEnt[jd];
							int a = s2->edges[j][0];
							int b = s2->edges[j][1];
							if (aux->uf.sameClass(a,b)==false){
								 // (origem,destino) sai
								 // (a,b) entra
								double new_objj_0 = aux->getObj(0)-custos[0][origem][destino]+custos[0][a][b];
								double new_objj_1 = aux->getObj(1)-custos[1][origem][destino]+custos[1][a][b];
								double new_objj_2 = aux->getObj(2)-custos[2][origem][destino]+custos[2][a][b];
								double new_objj_3 = aux->getObj(3)-custos[3][origem][destino]+custos[3][a][b];
								
								//atualizamos o Z - OK!
								if (new_objj_0>Z_max[0]) Z_max[0]=new_objj_0;
								if (new_objj_0<Z_min[0]) Z_min[0]=new_objj_0;
								if (new_objj_1>Z_max[1]) Z_max[1]=new_objj_1;
								if (new_objj_1<Z_min[1]) Z_min[1]=new_objj_1;
								if (new_objj_2>Z_max[2]) Z_max[2]=new_objj_2;
								if (new_objj_2<Z_min[2]) Z_min[2]=new_objj_2;
								if (new_objj_3>Z_max[3]) Z_max[3]=new_objj_3;
								if (new_objj_3<Z_min[3]) Z_min[3]=new_objj_3;
								

								////////

								
								// distância de aux para o ponto ideal
								double distance = sqrt((new_objj_0-ideal[subproblemaATUAL][0])*(new_objj_0-ideal[subproblemaATUAL][0]) + (new_objj_1-ideal[subproblemaATUAL][1])*(new_objj_1-ideal[subproblemaATUAL][1])   +   (new_objj_2-ideal[subproblemaATUAL][2])*(new_objj_2-ideal[subproblemaATUAL][2])   +    (new_objj_3-ideal[subproblemaATUAL][3])*(new_objj_3-ideal[subproblemaATUAL][3]));
								
								if (distance<min_objetivo){
									min_objetivo = distance;
									entra = jd;
							 		sai = iii;
								}

								//NOVO: atualiza o ponto ideal
								if (new_objj_0<ideal[subproblemaATUAL][0]) ideal[subproblemaATUAL][0] = new_objj_0;
								if (new_objj_1<ideal[subproblemaATUAL][1]) ideal[subproblemaATUAL][1] = new_objj_1;
								if (new_objj_2<ideal[subproblemaATUAL][2]) ideal[subproblemaATUAL][2] = new_objj_2;
								if (new_objj_3<ideal[subproblemaATUAL][3]) ideal[subproblemaATUAL][3] = new_objj_3;
								

							}
						}
					}
				}

				if (sai != -1 && entra != -1){
					int i_sai = amostralSai[sai]; // aresta i_sai que sai
					int o_sai = aux->edges[i_sai][0];
					int d_sai = aux->edges[i_sai][1];
					amostralSai[sai] = -1;

					int j_entra = amostralEnt[entra]; // aresta j_entra que entra
					int o_entra = s2->edges[j_entra][0];
					int d_entra = s2->edges[j_entra][1];

					aux->edges[i_sai][0] = o_entra;
					aux->edges[i_sai][1] = d_entra;
					aux->setObj(0, aux->getObj(0)-custos[0][o_sai][d_sai]+custos[0][o_entra][d_entra]);
					aux->setObj(1, aux->getObj(1)-custos[1][o_sai][d_sai]+custos[1][o_entra][d_entra]);
					aux->setObj(2, aux->getObj(2)-custos[2][o_sai][d_sai]+custos[2][o_entra][d_entra]);
					aux->setObj(3, aux->getObj(3)-custos[3][o_sai][d_sai]+custos[3][o_entra][d_entra]);
									
					/// insere a nova soluçao (aux) no subproblma correspondente
					// remove a soluçao de maior distância para o ponto ideal
					retorno->adicionarSol(aux); 

					double cosine_max = INT_MIN;
					int subprol = -1;
					for (int subk=0; subk<NUMSUMPROBLEMAS; subk++){

						double norma_lambda = lambda[subk][0]*lambda[subk][0] + lambda[subk][1]*lambda[subk][1]  + lambda[subk][2]*lambda[subk][2]   +  lambda[subk][3]*lambda[subk][3];
						double norma_aux = (aux->getObj(0)-Z_min[0])*(aux->getObj(0)-Z_min[0]) + (aux->getObj(1)-Z_min[1])*(aux->getObj(1)-Z_min[1])  +  (aux->getObj(2)-Z_min[2])*(aux->getObj(2)-Z_min[2])    + (aux->getObj(3)-Z_min[3])*(aux->getObj(3)-Z_min[3]);; 
						double produto_escalar = (lambda[subk][0]*(aux->getObj(0) - Z_min[0]))   +  (lambda[subk][1]*(aux->getObj(1) - Z_min[1]))  +  (lambda[subk][2]*(aux->getObj(2) - Z_min[2])) + (lambda[subk][3]*(aux->getObj(3) - Z_min[3]));
						double cosino_novo = produto_escalar/(sqrt(norma_aux)*sqrt(norma_lambda));

						if (cosino_novo>cosine_max){
							cosine_max = cosino_novo;
							subprol = subk;
						}
					}
					add_solution_subespaco(aux, subprol); // criar copia de aux
					
					//NOVO: atualiza o ponto ideal
								if (aux->getObj(0)<ideal[subprol][0]) ideal[subprol][0] = aux->getObj(0);
								if (aux->getObj(1)<ideal[subprol][1]) ideal[subprol][1] = aux->getObj(1);
								if (aux->getObj(2)<ideal[subprol][2]) ideal[subprol][2] = aux->getObj(2);
								if (aux->getObj(3)<ideal[subprol][3]) ideal[subprol][3] = aux->getObj(3);
								
					contAval++;
					contttttt++;
					// contabiliza avaliaçoes da FO aqui	

				} else {
					cout<<"ERROOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO -1"<<endl;
					exit(-1);

				}
			}

		}
public:




	BoundedParetoSet *MOPR_D(int &contPares){
		int contAval = 0;
		SolucaoEdgeSet *origm_p = new SolucaoEdgeSet(NUMEROVERTICES-1);
		SolucaoEdgeSet *dest_p = new SolucaoEdgeSet(NUMEROVERTICES-1);
		do{
			for (int i=0; i<NUMSUMPROBLEMAS  && contAval<NUM_AVALIACOES; i++){  
				subproblemaATUAL = i;
				int sub_a = B_vizinhos[i][IRandom(0,NUM_VIZINHOS-1)].subproblema;
				int sub_b = B_vizinhos[i][IRandom(0,NUM_VIZINHOS-1)].subproblema;
				
				int a1 = IRandom(0,subespacos[sub_a].size()-1);
				int b1 = IRandom(0,subespacos[sub_b].size()-1);
				int a2 = IRandom(0,subespacos[sub_a].size()-1);
				int b2 = IRandom(0,subespacos[sub_b].size()-1);
				
				list<SolucaoEdgeSet * >::iterator it_a1 = subespacos[sub_a].begin();
				for ( int ont = 0; ont<a1; ont++) it_a1++;
				
				list<SolucaoEdgeSet * >::iterator it_b1 = subespacos[sub_b].begin();
				for ( int ont = 0; ont<b1; ont++) it_b1++;

				list<SolucaoEdgeSet * >::iterator it_a2 = subespacos[sub_a].begin();
				for ( int ont = 0; ont<a2; ont++) it_a2++;
				
				list<SolucaoEdgeSet * >::iterator it_b2 = subespacos[sub_b].begin();
				for ( int ont = 0; ont<b2; ont++) it_b2++;
				
				double dis_a1 = sqrt(((*it_a1)->getObj(0)-ideal[sub_a][0])*((*it_a1)->getObj(0)-ideal[sub_a][0]) + ((*it_a1)->getObj(1)-ideal[sub_a][1])*((*it_a1)->getObj(1)-ideal[sub_a][1])  + ((*it_a1)->getObj(2)-ideal[sub_a][2])*((*it_a1)->getObj(2)-ideal[sub_a][2]) + ((*it_a1)->getObj(3)-ideal[sub_a][3])*((*it_a1)->getObj(3)-ideal[sub_a][3]));
				double dis_b1 = sqrt(((*it_b1)->getObj(0)-ideal[sub_b][0])*((*it_b1)->getObj(0)-ideal[sub_b][0]) + ((*it_b1)->getObj(1)-ideal[sub_b][1])*((*it_b1)->getObj(1)-ideal[sub_b][1])  + ((*it_b1)->getObj(2)-ideal[sub_b][2])*((*it_b1)->getObj(2)-ideal[sub_b][2]) + ((*it_b1)->getObj(3)-ideal[sub_b][3])*((*it_b1)->getObj(3)-ideal[sub_b][3]));
				double dis_a2 = sqrt(((*it_a2)->getObj(0)-ideal[sub_a][0])*((*it_a2)->getObj(0)-ideal[sub_a][0]) + ((*it_a2)->getObj(1)-ideal[sub_a][1])*((*it_a2)->getObj(1)-ideal[sub_a][1])  + ((*it_a2)->getObj(2)-ideal[sub_a][2])*((*it_a2)->getObj(2)-ideal[sub_a][2]) + ((*it_a2)->getObj(3)-ideal[sub_a][3])*((*it_a2)->getObj(3)-ideal[sub_a][3]));
				double dis_b2 = sqrt(((*it_b2)->getObj(0)-ideal[sub_b][0])*((*it_b2)->getObj(0)-ideal[sub_b][0]) + ((*it_b2)->getObj(1)-ideal[sub_b][1])*((*it_b2)->getObj(1)-ideal[sub_b][1])  + ((*it_b2)->getObj(2)-ideal[sub_b][2])*((*it_b2)->getObj(2)-ideal[sub_b][2]) + ((*it_b2)->getObj(3)-ideal[sub_b][3])*((*it_b2)->getObj(3)-ideal[sub_b][3]));
								

				SolucaoEdgeSet *origm, *dest;
				if (dis_a1 < dis_a2)  origm = *it_a1;
				else origm = *it_a2;

				if (dis_b1 < dis_b2)  dest = *it_b1;
				else dest = *it_b2; 

				*origm_p = *origm;
				*dest_p = *dest;
				PR(retorno, origm, dest, contAval);
				contPares++;
				PR(retorno, dest_p, origm_p, contAval);
				contPares++;
			}
			// cout<<"contAval = "<<contAval<<endl;
		}while(contAval<NUM_AVALIACOES);

		
	 	return retorno;

	}




	Version5_D(BoundedParetoSet *approx_set, double lambda2[NUMSUMPROBLEMAS][NUMOBJETIVOS]){
		for (int i=0; i<NUMSUMPROBLEMAS; i++){
			for (int j=0; j<NUMOBJETIVOS; j++){ 
				lambda[i][j] = lambda2[i][j];
				ideal[i][j] = 0.0;
			}
		}
		aux= new SolucaoEdgeSet(NUMEROVERTICES-1);
		s2 = new SolucaoEdgeSet(NUMEROVERTICES-1);
		retorno = new BoundedParetoSet();
		initialize(approx_set);
		

	}

};


#endif