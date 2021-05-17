#ifndef VERSION5_D_CPP
#define VERSION5_D_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2019)
	
This code implements the first PathRelinking version for Multi-objective Quadratic Assignment Problem

=====================================================================================
*/

/// Version5_D-msc


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

typedef struct vviz
{
	int subproblema;
	double distance;
} Vizinho;

class Version5_D: public PathRelinking{

private:
	Solution *aux,*s2;

	std::vector<double> lambda[NUMSUMPROBLEMAS]; // o lambda para cada subproblema
	list<Solution * > subespacos[NUMSUMPROBLEMAS];
	std::vector<double> ideal[NUMSUMPROBLEMAS]; // um ponto ideal para cada subproblema
	std::vector<double>  Z_min;
	std::vector<double>  Z_max;
	std::vector<double>  Z_min_sub[NUMSUMPROBLEMAS]; // um z_max para cada subproblema
	Vizinho B_vizinhos[NUMSUMPROBLEMAS][NUM_VIZINHOS]; // para cada subproblema (lambdas), temos os seus vizinhos
	BoundedParetoSet *retorno;
	int subproblemaATUAL;

	void calculaZ_por_subespaco(int sub){
		for (int obj = 0; obj < n_objetivos; ++obj){
			Z_min_sub[sub][obj] = INT_MAX;
		}
		for (list<Solution *>::iterator i = subespacos[sub].begin(); i!=subespacos[sub].end(); i++){
			for (int obj = 0; obj < n_objetivos; ++obj){
				if ((*i)->getObj(obj)<Z_min_sub[sub][obj]) Z_min_sub[sub][obj] = (*i)->getObj(obj);
			}
		}
	}

	void calculaZ_inicial(BoundedParetoSet *approx_set){
		for (int obj = 0; obj < n_objetivos; ++obj){
			Z_min[obj] = INT_MAX;
			Z_max[obj] = INT_MIN;
		}
		for ( list<Solution *>::iterator it_1=approx_set->sol.begin(); it_1!=approx_set->sol.end(); it_1++){
			for (int obj = 0; obj < n_objetivos; ++obj){
				if ((*it_1)->getObj(obj)<Z_min[obj]) Z_min[obj] = (*it_1)->getObj(obj);
				if ((*it_1)->getObj(obj)>Z_max[obj]) Z_max[obj] = (*it_1)->getObj(obj);
			}
		}
	}	



	// pre-condiçao : 0 <= subpro < NUMSUMPROBLEMAS
 	double delta_1(Solution *fx, int subpro){
		double norma_fx = 0.0;
		double norma_lambda = 0.0;
		double produto_escalar = 0.0;

		for (int i=0; i<n_objetivos; i++){
			norma_lambda += lambda[subpro][i]*lambda[subpro][i];
			norma_fx += ( fx->getObj(i)-Z_min[i])*( fx->getObj(i)-Z_min[i]); 
			produto_escalar += (lambda[subpro][i]*abs(fx->getObj(i)-Z_min[i]));// *(lambda[subpro][i]*(objj_normalizado-0.0));
		}

		double norma_lambda_aux = sqrt(norma_lambda);
		double norma_fx_aux = sqrt(norma_fx);
		return (double) produto_escalar/(norma_fx_aux*norma_lambda_aux);
	}




		// insere nova_solu se ela nao for dominada. As dominadas saem. Se atingir o tamanho maximo da subpopulaçao, retira-se a de maior distância para o ponto ideal
		bool add_solution_subespaco(Solution *nova_solu, int subb){
			
			/* percorre o vetor de solucoes e de valores e, caso exista solucao dominada, retira e retorna true. caso contrario, retorna false */
			list<Solution *>::iterator i = subespacos[subb].begin();
			list< list<Solution *>::iterator > remover;
			list<Solution *>::iterator maisPopuloso;
			double max_distance_i = INT_MIN;
			while (i != subespacos[subb].end()) {
				// se a solucao que vai entrar domina a solucao i 
				if (*nova_solu >> **i) {
					remover.push_back(i); 
				}
				
				double distance_it = 0.0;
				for (int obj = 0; obj < n_objetivos; ++obj)  distance_it += ((*i)->getObj(obj)-ideal[subb][obj])*((*i)->getObj(obj)-ideal[subb][obj]);
				
				distance_it = sqrt(distance_it);
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
				
				double distance_nova_solu = 0.0;
				for (int obj = 0; obj < n_objetivos; ++obj) distance_nova_solu+= (nova_solu->getObj(obj)-ideal[subb][obj])*(nova_solu->getObj(obj)-ideal[subb][obj]);
				distance_nova_solu = sqrt(distance_nova_solu);
				if (distance_nova_solu>max_distance_i) return false;

				remover.push_back(maisPopuloso);
			}

			list< list<Solution *>::iterator >::iterator j = remover.begin();
			while (j != remover.end()) {
				delete( **j ); 
				// remove do conjunto pareto
				subespacos[subb].erase( *j );
				j++;
			}

			Solution *t = new Solution(n_locations,n_objetivos);
			*t = *nova_solu; /// CRIA COPIA
			// adiciona ao conjunto pareto
			subespacos[subb].push_front( t );
			// if (subespacos[subb].size() > MAX_SUB_POPULATION) {
			// 	fprintf(stderr,"ERRO in MAX_SUB_POPULATION!\n");
			// 	exit(-1);
			// }
			return true;		
		}



	void initialize(BoundedParetoSet *approx_set){
		calculaZ_inicial(approx_set);
		// divide os subespaços
		for ( list<Solution *>::iterator it_1=approx_set->sol.begin(); it_1!=approx_set->sol.end(); it_1++){
			retorno->adicionarSol(*it_1); // copia, previamente, o conjunto inicial.
		
			//calcula os sub-espaços objetivo para cada subproblema
			double cosine_max = INT_MIN;
			int subpro = -1;
			for (int i=0; i<NUMSUMPROBLEMAS; i++){
				
				double cosine = delta_1((*it_1), i);
				if (cosine>cosine_max){
					cosine_max = cosine;
					subpro = i;
				}
			}
			if (subpro!=-1){
				add_solution_subespaco(*it_1, subpro);
			}
		}


		// se algum sub-espaço estiver vazio, entao inserimos nele a soluçao de menor ângulo 
		// QUESTAO: E se fizéssimos assim para todos os sub-espaços?????
		for (int i=0; i<NUMSUMPROBLEMAS; i++){

			if (subespacos[i].size()==0){
				double cosine_max = INT_MIN;
				list<Solution *>::iterator sol_mais_prox;
				for ( list<Solution *>::iterator it_1=approx_set->sol.begin(); it_1!=approx_set->sol.end(); it_1++){
					double cosine = delta_1((*it_1), i);
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
			for (int obj = 0; obj<n_objetivos; obj++){
				ideal[i][obj] = Z_min_sub[i][obj]; 
			}
		}

		

		// calcula os NUM_VIZINHOS mais proximos de cada vetor i
		for (int i=0; i<NUMSUMPROBLEMAS; i++){
			int contVizinhos = 0;
			for (int j=0; j<NUMSUMPROBLEMAS; j++){
				if (j!=i){
					double sumDistance =  0;
					for (int obj=0; obj<n_objetivos; obj++){
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




	void PR(BoundedParetoSet *retorno, Solution *s11, Solution *s222, int &contAval){
		*aux = *s11; // O(n)
		*s2 = *s222; // O(n)
		std::vector<pair<int, int> > paresTroca; // troca o first com o second


		do{
			paresTroca.clear();

			// ============================= Computa pares O(n^2) ============================= //
			for (int i=0; i<n_locations; i++){ // varre a origem (aux)
				// para cada facility que está na posicao i em aux e na posicao j em s2, tal que i!=j, entao cria-se o par (i,j) em paresTroca
				if (aux->getAssignment(i)!=s2->getAssignment(i)) { // isso assegura que a posicao de aux->getAssignment(i) em s2 não é i
					// procura-se a posicao j de aux->getAssignment(i) em s2
					for (int j = 0; j < n_locations; ++j){
						if (aux->getAssignment(i)==s2->getAssignment(j)){
							// devido o if (aux->getAssignment(i)!=s2->getAssignment(i)), pode-se assegurar que i!=j
							paresTroca.push_back(make_pair(i,j)); // O(n) insercoes
						}
					}
				}
			}

			// ============================= Critério de escolha ============================= //
			if (paresTroca.size()>0){

				double minDist = INT_MAX;
				int i = -1;
				int j = -1;

				for (int opcao = 0; opcao < paresTroca.size(); ++opcao){ // Essa linha é O(n) (no maximo, n valores em aux estarao fora do lugar). Total do laço O(qn^2), pois o delta eh O(n) 
					// cout<<"troca facilidades = "<<auxiliar->getAssignment(i)<<" "<<auxiliar->getAssignment(j)<<endl;
					
					double distance = 0.0;

					for (int obj = 0; obj < n_objetivos; ++obj){
						double obj_novo = aux->getObj(obj)+delta(aux->getAssignment_(),obj,paresTroca[opcao].first,paresTroca[opcao].second); // i j ou j i
						
						if (obj_novo > Z_max[obj]) Z_max[obj]=obj_novo;
						if (obj_novo < Z_min[obj]) Z_min[obj]=obj_novo;

						// distância de aux para o ponto ideal
						distance+=(obj_novo-ideal[subproblemaATUAL][obj])*(obj_novo-ideal[subproblemaATUAL][obj]);

						//atualiza o ponto ideal
						if (obj_novo < ideal[subproblemaATUAL][obj]) ideal[subproblemaATUAL][obj] = obj_novo;
								
					}


					distance = sqrt(distance);
					if (distance < minDist){
						minDist = distance;
						i = paresTroca[opcao].first;
						j = paresTroca[opcao].second;
					}
				
				}

				int a_i =  aux->getAssignment(i);
				aux->setAssignment_(i,aux->getAssignment(j));
				aux->setAssignment_(j,a_i);
				aux->calculaObjetivos(contAval); // contAval++
				retorno->adicionarSol(aux);





				////////
				double cosine_max = INT_MIN;
				int subprol = -1;

				for (int subk=0; subk<NUMSUMPROBLEMAS; subk++){
					
					double norma_aux = 0.0;
					double norma_lambda = 0.0;
					double produto_escalar = 0.0;

					for (int obj=0; obj<n_objetivos; obj++){
						norma_lambda += lambda[subk][obj]*lambda[subk][obj];
						norma_aux += (aux->getObj(obj)-Z_min[obj])*( aux->getObj(obj)-Z_min[obj]); 
						produto_escalar += (lambda[subk][obj]*abs(aux->getObj(obj)-Z_min[obj]));
					}

					
					double cosino_novo = (double) produto_escalar/(sqrt(norma_aux)*sqrt(norma_lambda));

					if (cosino_novo>cosine_max){
						cosine_max = cosino_novo;
						subprol = subk;
					}
				}
				add_solution_subespaco(aux, subprol); // criar copia de aux


			}

		} while (paresTroca.size()>0);
		
	}








public:





	BoundedParetoSet *MOPR_D(int &contPares){
		int contAval = 0;
		Solution *origm_p = new Solution(n_locations,n_objetivos);
		Solution  *dest_p = new Solution(n_locations,n_objetivos);
		do{
			for (int i=0; i<NUMSUMPROBLEMAS  && contAval<NUM_AVALIACOES; i++){  
				subproblemaATUAL = i;
				int sub_a = B_vizinhos[i][pick_a_number(0,NUM_VIZINHOS-1)].subproblema;
				int sub_b = B_vizinhos[i][pick_a_number(0,NUM_VIZINHOS-1)].subproblema;
				
				int a1 = pick_a_number(0,subespacos[sub_a].size()-1);
				int b1 = pick_a_number(0,subespacos[sub_b].size()-1);
				int a2 = pick_a_number(0,subespacos[sub_a].size()-1);
				int b2 = pick_a_number(0,subespacos[sub_b].size()-1);
				
				list<Solution * >::iterator it_a1 = subespacos[sub_a].begin();
				for ( int ont = 0; ont<a1; ont++) it_a1++;
				
				list<Solution * >::iterator it_b1 = subespacos[sub_b].begin();
				for ( int ont = 0; ont<b1; ont++) it_b1++;

				list<Solution * >::iterator it_a2 = subespacos[sub_a].begin();
				for ( int ont = 0; ont<a2; ont++) it_a2++;
				
				list<Solution * >::iterator it_b2 = subespacos[sub_b].begin();
				for ( int ont = 0; ont<b2; ont++) it_b2++;
				
				double dis_a1 = 0.0;
				double dis_b1 = 0.0;
				double dis_a2 = 0.0;
				double dis_b2 = 0.0;
				for (int obj = 0; obj < n_objetivos; ++obj){
					dis_a1 += ((*it_a1)->getObj(obj)-ideal[sub_a][obj])*((*it_a1)->getObj(obj)-ideal[sub_a][obj]);
					dis_b1 += ((*it_b1)->getObj(obj)-ideal[sub_b][obj])*((*it_b1)->getObj(obj)-ideal[sub_b][obj]);
					dis_a2 += ((*it_a2)->getObj(obj)-ideal[sub_a][obj])*((*it_a2)->getObj(obj)-ideal[sub_a][obj]);
					dis_b2 += ((*it_b2)->getObj(obj)-ideal[sub_b][obj])*((*it_b2)->getObj(obj)-ideal[sub_b][obj]);
				}
				dis_a1 = sqrt(dis_a1);
				dis_b1 = sqrt(dis_b1);
				dis_a2 = sqrt(dis_a2);
				dis_b2 = sqrt(dis_b2);
				// double dis_a1 = sqrt(((*it_a1)->getObj(0)-ideal[sub_a][0])*((*it_a1)->getObj(0)-ideal[sub_a][0]) + ((*it_a1)->getObj(1)-ideal[sub_a][1])*((*it_a1)->getObj(1)-ideal[sub_a][1])  + ((*it_a1)->getObj(2)-ideal[sub_a][2])*((*it_a1)->getObj(2)-ideal[sub_a][2]) );
				// double dis_b1 = sqrt(((*it_b1)->getObj(0)-ideal[sub_b][0])*((*it_b1)->getObj(0)-ideal[sub_b][0]) + ((*it_b1)->getObj(1)-ideal[sub_b][1])*((*it_b1)->getObj(1)-ideal[sub_b][1])  + ((*it_b1)->getObj(2)-ideal[sub_b][2])*((*it_b1)->getObj(2)-ideal[sub_b][2]) );
				// double dis_a2 = sqrt(((*it_a2)->getObj(0)-ideal[sub_a][0])*((*it_a2)->getObj(0)-ideal[sub_a][0]) + ((*it_a2)->getObj(1)-ideal[sub_a][1])*((*it_a2)->getObj(1)-ideal[sub_a][1])  + ((*it_a2)->getObj(2)-ideal[sub_a][2])*((*it_a2)->getObj(2)-ideal[sub_a][2]) );
				// double dis_b2 = sqrt(((*it_b2)->getObj(0)-ideal[sub_b][0])*((*it_b2)->getObj(0)-ideal[sub_b][0]) + ((*it_b2)->getObj(1)-ideal[sub_b][1])*((*it_b2)->getObj(1)-ideal[sub_b][1])  + ((*it_b2)->getObj(2)-ideal[sub_b][2])*((*it_b2)->getObj(2)-ideal[sub_b][2]) );
								

				Solution *origm, *dest;
				if (dis_a1 < dis_a2)  origm = *it_a1;
				else origm = *it_a2;

				if (dis_b1 < dis_b2)  dest = *it_b1;
				else dest = *it_b2; 

				*origm_p = *origm;
				*dest_p = *dest;
				PR(retorno, origm, dest, contAval);
				// return retorno; //**************** REMOVE ****************
				contPares++;
				PR(retorno, dest_p, origm_p, contAval);
				contPares++;
			}
			// cout<<"contAval = "<<contAval<<endl;
		}while(contAval<NUM_AVALIACOES);

		
	 	return retorno;

	}



	
	Version5_D(BoundedParetoSet *approx_set, std::vector<double> lambda1[NUMSUMPROBLEMAS]){
		aux = new  Solution(n_locations,n_objetivos);
		s2 = new  Solution(n_locations,n_objetivos);
		for (int i=0; i<NUMSUMPROBLEMAS; i++){
			ideal[i].resize(n_objetivos);
			Z_min_sub[i].resize(n_objetivos);
			lambda[i] = lambda1[i];
		}
		retorno = new BoundedParetoSet();
		Z_max.resize(n_objetivos); 
		Z_min.resize(n_objetivos);
		initialize(approx_set);


		// for (int i=0; i<NUMSUMPROBLEMAS; i++){

		// 	cout<<"Subproblema "<<i<<" ("<< lambda[i][0] <<" "<<lambda[i][1]<<" "<<lambda[i][2] <<")"<<endl; //( "<< bestSolution[i]->getObj(0)<<" "<< bestSolution[i]->getObj(1)<<"): "<<endl;
			
		// 	for (list<Solution *>::iterator it_1=subespacos[i].begin(); it_1!=subespacos[i].end(); it_1++){
		// 		// cout<<"\t"<<(*it_1)->getObj(0)<<" "<<(*it_1)->getObj(1)<<" "<<(*it_1)->getObj(2)<<endl;
		// 		fprintf(stdout,"\t %f %f %f\n", (*it_1)->getObj(0), (*it_1)->getObj(1), (*it_1)->getObj(2)); 
		// 	}
		// 	cout<<"ideal = "<<ideal[i][0]<<" "<<ideal[i][1]<<" "<<ideal[i][2]<<endl;
		// 	for (int viz=0; viz<NUM_VIZINHOS; viz++){
		// 		cout<<B_vizinhos[i][viz].subproblema<<" ";
		// 	}
		// 	cout<<endl;
		// }




	}

};


#endif
