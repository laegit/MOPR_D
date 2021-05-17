#ifndef VERSION5_D_CPP
#define VERSION5_D_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2018)
	
This code implements the 5th PathRelinking version for Multiobjective Multidimensional Knapsack problem 

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
#include "BoundedParetoSet.cpp"

typedef struct vviz
{
	int subproblema;
	double distance;
} Vizinho;

class Version5_D{//: public PathRelinking{

private:
	Knapsack *aux, *s2;
	std::vector<double> lambdas[NUMSUMPROBLEMAS]; // o lambda para cada subproblema
	list<Knapsack * > subespacos[NUMSUMPROBLEMAS];
	std::vector<double> ideal[NUMSUMPROBLEMAS]; // um ponto ideal para cada subproblema
	std::vector<double>  Z_min;
	std::vector<double>  Z_max;
	std::vector<double>  Z_max_sub[NUMSUMPROBLEMAS]; // um z_max para cada subproblema
	Vizinho B_vizinhos[NUMSUMPROBLEMAS][NUM_VIZINHOS]; // para cada subproblema (lambdas), temos os seus vizinhos
	BoundedParetoSet *retorno;
	int subproblemaATUAL;




	// um ponto ideal para cada sub-espaço
	void calculaZ_por_subespaco(int sub){
		for (int obj = 0; obj < my_momkp->dimension_; ++obj){
			Z_max_sub[sub][obj] = INT_MIN;
		}
		for (list<Knapsack *>::iterator i = subespacos[sub].begin(); i!=subespacos[sub].end(); i++){
			for (int obj = 0; obj < my_momkp->dimension_; ++obj){
				if ((*i)->profit(obj) > Z_max_sub[sub][obj]) Z_max_sub[sub][obj] = (*i)->profit(obj);
			}
		}
	}


	void calculaZ_inicial(BoundedParetoSet *approx_set){
		for (int obj = 0; obj < my_momkp->dimension_; ++obj){
			Z_min[obj] = INT_MAX;
			Z_max[obj] = INT_MIN;
		}
		for ( list<Knapsack *>::iterator it_1=approx_set->sol.begin(); it_1!=approx_set->sol.end(); it_1++){
			for (int obj = 0; obj < my_momkp->dimension_; ++obj){
				if ((*it_1)->profit(obj)<Z_min[obj]) Z_min[obj] = (*it_1)->profit(obj);
				if ((*it_1)->profit(obj)>Z_max[obj]) Z_max[obj] = (*it_1)->profit(obj);
			}
		}
	}	


	// pre-condiçao : 0 <= subpro < NUMSUMPROBLEMAS
 	double delta(Knapsack *fx, int subpro){
		double norma_fx = 0.0;
		double norma_lambda = 0.0;
		double produto_escalar = 0.0;

		for (int i=0; i<my_momkp->dimension_; i++){
			norma_lambda += lambdas[subpro][i]*lambdas[subpro][i];
			norma_fx += abs( fx->profit(i)-Z_max[i])*abs( fx->profit(i)-Z_max[i]); 
			produto_escalar += (lambdas[subpro][i]*abs(fx->profit(i)-Z_max[i]));// 
		}

		double norma_lambda_aux = sqrt(norma_lambda);
		double norma_fx_aux = sqrt(norma_fx);
		return (double) produto_escalar/(norma_fx_aux*norma_lambda_aux);
	}



	// insere nova_solu se ela nao for dominada. As dominadas saem. Se atingir o tamanho maximo da subpopulaçao, retira-se a de maior distância para o ponto ideal
	bool add_solution_subespaco(Knapsack *nova_solu, int subb){
		
		/* percorre o vetor de solucoes e de valores e, caso exista solucao dominada, retira e retorna true. caso contrario, retorna false */
		list<Knapsack *>::iterator i = subespacos[subb].begin();
		list< list<Knapsack *>::iterator > remover;
		list<Knapsack *>::iterator maisPopuloso;
		double max_distance_i = INT_MIN;
		while (i != subespacos[subb].end()) {
			// se a solucao que vai entrar domina a solucao i 
			if (*nova_solu >> **i) {
				remover.push_back(i); 
			}
			
			double distance_it = 0.0;
			for (int obj = 0; obj < my_momkp->dimension_; ++obj){
				distance_it += ((*i)->profit(obj)-ideal[subb][obj])*((*i)->profit(obj)-ideal[subb][obj]);
			}
			distance_it = sqrt(distance_it);
			// = sqrt(((*i)->profit(0)-ideal[subb][0])*((*i)->profit(0)-ideal[subb][0]) + ((*i)->getObj(1)-ideal[subb][1])*((*i)->getObj(1)-ideal[subb][1])  +  ((*i)->getObj(2)-ideal[subb][2])*((*i)->getObj(2)-ideal[subb][2])   +  ((*i)->getObj(3)-ideal[subb][3])*((*i)->getObj(3)-ideal[subb][3]) );
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
			// sqrt((nova_solu->getObj(0)-ideal[subb][0])*(nova_solu->getObj(0)-ideal[subb][0]) + (nova_solu->getObj(1)-ideal[subb][1])*(nova_solu->getObj(1)-ideal[subb][1]) +  (nova_solu->getObj(2)-ideal[subb][2])*(nova_solu->getObj(2)-ideal[subb][2])  +   (nova_solu->getObj(3)-ideal[subb][3])*(nova_solu->getObj(3)-ideal[subb][3]));
			for (int obj = 0; obj < my_momkp->dimension_; ++obj){
				distance_nova_solu += (nova_solu->profit(obj)-ideal[subb][obj])*(nova_solu->profit(obj)-ideal[subb][obj]);
			}
			distance_nova_solu = sqrt(distance_nova_solu);
			if (distance_nova_solu>max_distance_i) return false;

			remover.push_back(maisPopuloso);
		}

		list< list<Knapsack *>::iterator >::iterator j = remover.begin();
		while (j != remover.end()) {
			delete( **j ); 
			// remove do conjunto pareto
			subespacos[subb].erase( *j );
			j++;
		}

		Knapsack *t = new Knapsack(my_momkp->size_, my_momkp->dimension_);
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
		for ( list<Knapsack *>::iterator it_1=approx_set->sol.begin(); it_1!=approx_set->sol.end(); it_1++){
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
				add_solution_subespaco(*it_1, subpro);
			}
		}


		// se algum sub-espaço estiver vazio, entao inserimos nele a soluçao de menor ângulo 
		// QUESTAO: E se fizéssimos assim para todos os sub-espaços?????
		for (int i=0; i<NUMSUMPROBLEMAS; i++){

			if (subespacos[i].size()==0){
				double cosine_max = INT_MIN;
				list<Knapsack *>::iterator sol_mais_prox;
				for ( list<Knapsack *>::iterator it_1=approx_set->sol.begin(); it_1!=approx_set->sol.end(); it_1++){
					double cosine = delta((*it_1), i);
					if (cosine>cosine_max){
						cosine_max = cosine;
						sol_mais_prox = it_1;
					}
				}
				add_solution_subespaco(*sol_mais_prox, i);
			}
		}
		

		//inicializa o ponto ideal para cada subespaço (onde Z_max_sub é calculado dentro de cada sub-espaço)
		for (int i = 0; i<NUMSUMPROBLEMAS; i++){
			calculaZ_por_subespaco(i);
			for (int obj = 0; obj<my_momkp->dimension_; obj++){
				ideal[i][obj] = Z_max_sub[i][obj]; 
			}
		}

		

		// calcula os NUM_VIZINHOS mais proximos de cada vetor i
		for (int i=0; i<NUMSUMPROBLEMAS; i++){
			int contVizinhos = 0;
			for (int j=0; j<NUMSUMPROBLEMAS; j++){
				if (j!=i){
					double sumDistance =  0;
					for (int obj=0; obj<my_momkp->dimension_; obj++){
						sumDistance+=(lambdas[i][obj]-lambdas[j][obj])*(lambdas[i][obj]-lambdas[j][obj]);
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



	void PR(ParetoSet *retorno, Knapsack *s111, Knapsack *s222, int &contAval){
		// para todo elemento em s1 cujo bit é diferente do bit correspondnte em s2. 
		// Ou seja, para todo elemento em s1 que nao está em s2 e vice-versa
		*aux = *s111;
		*s2 = *s222;
		bool continua = true;

		// retorno->adicionarSol(aux);	
		
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

			if (indicesEmAux.size()>1){
				// cout<<"indicesEmAux.size() = "<<indicesEmAux.size()<<endl;
				continua = true;
				int item = -1; 
				double min_objetivo = INT_MAX; // minimiza a distância ao ponto ideal
				
				for (int jjt = 0; jjt<indicesEmAux.size(); jjt++){ // para cada elemento
					int auxItem = indicesEmAux[jjt];
					// std::vector<double> novos_objetivos(my_momkp->dimension_);
					int multiplicador = 1; // // simula adiçao do item indicesEmAux[jjt]
					double distance = 0.0;
					if (aux->at(auxItem)==true){ // simula remoçao do item indicesEmAux[jjt]
						multiplicador = -1;
					}
					for (int pl = 0; pl<my_momkp->dimension_; pl++){
						double novo_obj_pl = (double) aux->profit(pl)+(my_momkp->profit_.at(auxItem).at(pl)*multiplicador);
						// novos_objetivos.at(novo_obj_pl);
						
						//atualizamos o Z - OK!
						if (novo_obj_pl>Z_max[pl]) Z_max[pl]=novo_obj_pl;
						if (novo_obj_pl<Z_min[pl]) Z_min[pl]=novo_obj_pl;


						// distância de aux para o ponto ideal
						distance += (novo_obj_pl-ideal[subproblemaATUAL][pl])*(novo_obj_pl-ideal[subproblemaATUAL][pl]);
					
						//NOVO: atualiza o ponto ideal
						if (novo_obj_pl>ideal[subproblemaATUAL][pl]) ideal[subproblemaATUAL][pl] = novo_obj_pl;
						
					}
					distance = sqrt(distance);
					if (distance<min_objetivo){
						min_objetivo = distance;
						item = auxItem;
					}
				}
				

				if (item!=-1){
					if (aux->at(item)==false){ // insere
						aux->insert_item(item, my_momkp->profit_.at(item), my_momkp->weight_.at(item));
					}else{  // remove
						aux->remove_item(item, my_momkp->profit_.at(item), my_momkp->weight_.at(item));
					}
					// aux->print_solution(stdout);
					retorno->adicionarSol(aux);
					contAval++; // sempre que um item é inserido ou removido da mochila, precisa somar ou subtrarir algum valor da funçao objetivo. Por isso, contabiliza-se a quantidade de avaliaçoes aqui
					

					////////
					double cosine_max = INT_MIN;
					int subprol = -1;

					for (int subk=0; subk<NUMSUMPROBLEMAS; subk++){
						
						double norma_aux = 0.0;
						double norma_lambda = 0.0;
						double produto_escalar = 0.0;

						for (int i=0; i<my_momkp->dimension_; i++){
							norma_lambda += lambdas[subk][i]*lambdas[subk][i];
							norma_aux += abs( aux->profit(i)-Z_max[i])*abs( aux->profit(i)-Z_max[i]); 
							produto_escalar += (lambdas[subk][i]*abs(aux->profit(i)-Z_max[i]));// 
						}

						
						double cosino_novo = (double) produto_escalar/(sqrt(norma_aux)*sqrt(norma_lambda));

						if (cosino_novo>cosine_max){
							cosine_max = cosino_novo;
							subprol = subk;
						}
					}
					add_solution_subespaco(aux, subprol); // criar copia de aux
					


				} else {
					cout<<"ERROOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO -1"<<endl;
					exit(-1);
				}

			} else continua = false;
		} while (continua==true);
	}


public:




	BoundedParetoSet *MOPR_D(int &contPares){
		int contAval = 0;
		Knapsack *origm_p = new Knapsack(my_momkp->size_, my_momkp->dimension_);
		Knapsack *dest_p = new Knapsack(my_momkp->size_, my_momkp->dimension_);
		do{
			for (int i=0; i<NUMSUMPROBLEMAS && contAval<NUM_AVALIACOES; i++){
				subproblemaATUAL = i;
				int sub_a = B_vizinhos[i][pick_a_number(0,NUM_VIZINHOS-1)].subproblema;
				int sub_b = B_vizinhos[i][pick_a_number(0,NUM_VIZINHOS-1)].subproblema;
				
				int a1 = pick_a_number(0,subespacos[sub_a].size()-1);
				int b1 = pick_a_number(0,subespacos[sub_b].size()-1);
				int a2 = pick_a_number(0,subespacos[sub_a].size()-1);
				int b2 = pick_a_number(0,subespacos[sub_b].size()-1);
				
				list<Knapsack * >::iterator it_a1 = subespacos[sub_a].begin();
				for ( int ont = 0; ont<a1; ont++) it_a1++;
				
				list<Knapsack * >::iterator it_b1 = subespacos[sub_b].begin();
				for ( int ont = 0; ont<b1; ont++) it_b1++;

				list<Knapsack * >::iterator it_a2 = subespacos[sub_a].begin();
				for ( int ont = 0; ont<a2; ont++) it_a2++;
				
				list<Knapsack * >::iterator it_b2 = subespacos[sub_b].begin();
				for ( int ont = 0; ont<b2; ont++) it_b2++;
				
				double dis_a1 = 0.0;
				double dis_b1 = 0.0;
				double dis_a2 = 0.0;
				double dis_b2 = 0.0;
				for (int objj = 0; objj < my_momkp->dimension_; ++objj){
					dis_a1 += ((*it_a1)->profit(objj)-ideal[sub_a][objj])*((*it_a1)->profit(objj)-ideal[sub_a][objj]);
					dis_b1 += ((*it_b1)->profit(objj)-ideal[sub_b][objj])*((*it_b1)->profit(objj)-ideal[sub_b][objj]);
					dis_a2 += ((*it_a2)->profit(objj)-ideal[sub_a][objj])*((*it_a2)->profit(objj)-ideal[sub_a][objj]);
					dis_b2 += ((*it_b2)->profit(objj)-ideal[sub_b][objj])*((*it_b2)->profit(objj)-ideal[sub_b][objj]);
				}
				dis_a1 = sqrt(dis_a1);
				dis_b1 = sqrt(dis_b1);
				dis_a2 = sqrt(dis_a2);
				dis_b2 = sqrt(dis_b2);

				// double dis_a1 = sqrt(((*it_a1)->getObj(0)-ideal[sub_a][0])*((*it_a1)->getObj(0)-ideal[sub_a][0]) + ((*it_a1)->getObj(1)-ideal[sub_a][1])*((*it_a1)->getObj(1)-ideal[sub_a][1]) );
				// double dis_b1 = sqrt(((*it_b1)->getObj(0)-ideal[sub_b][0])*((*it_b1)->getObj(0)-ideal[sub_b][0]) + ((*it_b1)->getObj(1)-ideal[sub_b][1])*((*it_b1)->getObj(1)-ideal[sub_b][1]) );
				// double dis_a2 = sqrt(((*it_a2)->getObj(0)-ideal[sub_a][0])*((*it_a2)->getObj(0)-ideal[sub_a][0]) + ((*it_a2)->getObj(1)-ideal[sub_a][1])*((*it_a2)->getObj(1)-ideal[sub_a][1]) );
				// double dis_b2 = sqrt(((*it_b2)->getObj(0)-ideal[sub_b][0])*((*it_b2)->getObj(0)-ideal[sub_b][0]) + ((*it_b2)->getObj(1)-ideal[sub_b][1])*((*it_b2)->getObj(1)-ideal[sub_b][1]) );
								

				Knapsack *origm, *dest;
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



	Version5_D(BoundedParetoSet *approx_set, std::vector<double> lambda1[NUMSUMPROBLEMAS]){
		aux = new Knapsack(my_momkp->size_, my_momkp->dimension_);
		s2 = new Knapsack(my_momkp->size_, my_momkp->dimension_);
		for (int i=0; i<NUMSUMPROBLEMAS; i++){
			ideal[i].resize(my_momkp->dimension_);
			Z_max_sub[i].resize(my_momkp->dimension_);
			lambdas[i] = lambda1[i];
		}
		retorno = new BoundedParetoSet();
		Z_max.resize(my_momkp->dimension_); 
		Z_min.resize(my_momkp->dimension_);
		initialize(approx_set);
		// for (int i=0; i<NUMSUMPROBLEMAS; i++){

		// 	cout<<"Subproblema "<<i<<" ("<< lambdas[i][0] <<" "<<lambdas[i][1]<<" "<<lambdas[i][2] <<")"<<endl; //( "<< bestSolution[i]->getObj(0)<<" "<< bestSolution[i]->getObj(1)<<"): "<<endl;
			
		// 	for (list<Knapsack *>::iterator it_1=subespacos[i].begin(); it_1!=subespacos[i].end(); it_1++){
		// 		cout<<"\t"<<(*it_1)->profit(0)<<" "<<(*it_1)->profit(1)<<" "<<(*it_1)->profit(2)<<endl;
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
