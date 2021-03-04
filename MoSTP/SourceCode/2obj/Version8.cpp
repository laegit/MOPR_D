#ifndef VERSION8_CPP
#define VERSION8_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2018)
	
This code implements the second PathRelinking version for Bi-objective spanning tree 
The data structure and some functions were kindly provided by Monteiro (2010)


hypervolume-based selection 

baseado no trabalho do anyTimePLS

apenas para 2 objetivos

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



class Version8: public PathRelinking{
private:
	int amostralSai[NUMEROVERTICES-1]; // guarda indexes de arestas (de s1) que devem sair de s1 (vide abaixo)
	int amostralEnt[NUMEROVERTICES-1]; // guarda indexes de arestas (de s2) que entrar em s1 (vide abaixo)

	SolucaoEdgeSet *aux;
	std::list<auxSoluacao > nonDominadedSolutions;



	bool A_domina_B(double a_obj1, double a_obj2, double b_obj1, double b_obj2){
		return (a_obj1<=b_obj1 && a_obj2<=b_obj2 && (a_obj1<b_obj1 || a_obj2<b_obj2)); //minimizaçao
	}

	double ohvc(double obj1_s1, double obj2_s1, double obj1_s2, double obj2_s2){
		return (obj1_s1 - obj1_s2)*(obj2_s2 - obj2_s1);
	}

	double addNonDominandedSolution(double obj1, double obj2, int sai,  int entra){
		list<auxSoluacao>::iterator i = nonDominadedSolutions.begin();
		list< list<auxSoluacao>::iterator > remover;
		auxSoluacao s_sup,  s_inf;
		bool mod_sup= false, mod_inf = false;
		double min = 1e9;
		double max = -1;
		while (i != nonDominadedSolutions.end()) {
			// if (*nova >> *((*i).sol))
			if (A_domina_B(obj1,obj2, (*i).obj1, (*i).obj2))
				remover.push_back(i);
			// if (*((*i).sol) >> *nova)
			if (A_domina_B((*i).obj1, (*i).obj2, obj1,obj2) || ((*i).obj1==obj1 && (*i).obj2==obj2))
				return -1e9;
			

			//OHI
			if ((*i).obj2>obj2){
				if ((*i).obj2<min){
					min = (*i).obj2;
					s_sup = *i;
					mod_sup = true;
				}
			} else if ((*i).obj2 < obj2){
				if ((*i).obj2>max){
					max = (*i).obj2;
					s_inf = *i;
					mod_inf = true;
				}
			}


			i++;
		}

		list< list<auxSoluacao>::iterator >::iterator j = remover.begin();
		while (j != remover.end()) {
			nonDominadedSolutions.erase( *j );
			j++;
		}
		nonDominadedSolutions.push_front(auxSoluacao {obj1, obj2, sai, entra});
		

		if (mod_inf == false){
			return 2*ohvc(s_sup.obj1,s_sup.obj2, obj1, obj2);
		} else if (mod_sup == false){
			return 2*ohvc(obj1, obj2, s_inf.obj1, s_inf.obj2);
		} else {
			return ohvc(s_sup.obj1,s_sup.obj2, obj1, obj2)+ohvc(obj1, obj2, s_inf.obj1, s_inf.obj2);
		}
	}

	


	// double OHI(double obj1, double obj2, BoundedParetoSet *retorno){

	// 	list<SolucaoEdgeSet*> solucoes = retorno->getElementos();
	// 	SolucaoEdgeSet * s_sup = NULL,  *s_inf = NULL;


	// 	//get s_sup and s_inf
	// 	double min = 1e9;
	// 	double max = -1;
	// 	list<SolucaoEdgeSet*>::iterator it = solucoes.begin();
	// 	while (it!=solucoes.end()){
	// 		if ((*it)->getObj(1)>obj2){
	// 			if ((*it)->getObj(1)<min){
	// 				min = (*it)->getObj(1);
	// 				s_sup = *it;
	// 			}
	// 		} else if ((*it)->getObj(1) < obj2){
	// 			if ((*it)->getObj(1)>max){
	// 				max = (*it)->getObj(1);
	// 				s_inf = *it;
	// 			}
	// 		}
	// 		it++;
	// 	}


	// 	if (s_inf == NULL){
	// 		return 2*ohvc(s_sup->getObj(0),s_sup->getObj(1), obj1, obj2);
	// 	} else if (s_sup == NULL){
	// 		return 2*ohvc(obj1, obj2, s_inf->getObj(0), s_inf->getObj(1));
	// 	} else {
	// 		return ohvc(s_sup->getObj(0),s_sup->getObj(1), obj1, obj2)+ohvc(obj1, obj2, s_inf->getObj(0), s_inf->getObj(1));
	// 	}

	// }




	void PR(BoundedParetoSet *retorno, SolucaoEdgeSet *s1, SolucaoEdgeSet *s2, int &contAval){
	
		int contttttt = 0;
		// O(n^2)
		int totalSai = inS1_NotInS2(s1, s2, amostralSai); // arestas de s1 que nao estao em s2. Tais arestas devem sair de s1
		int totalEnt = inS1_NotInS2(s2, s1, amostralEnt); // arestas de s2 que nao estao em s1. Tais arestas devem entrar em s1

		*aux = *s1; // O(n)

		for (int cont = 0; cont<totalSai && contAval<NUM_AVALIACOES; cont++){

			double maxOHI = -1e9;
			int sai = -1;
			int entra = -1;
			
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

							 double newObj_0 = aux->getObj(0)-custos[0][origem][destino]+custos[0][a][b];
							 double newObj_1 = aux->getObj(1)-custos[1][origem][destino]+custos[1][a][b];
							 double newOHI = addNonDominandedSolution(newObj_0,newObj_1, iii, jd);
						// cout<<"newOHI = "<<newOHI<<" maxOHI= "<<maxOHI<<endl;
						 	 if (newOHI>maxOHI){
						 	 	maxOHI = newOHI;
							 	entra = jd;
							 	sai = iii;

							 
						 	 } 
						 	 // else if (nonDominadedSolutions.size()==0){
						 	 // 	maxOHI = -1e9;
						 	 // 	entra = jd;
							 	// sai = iii;
						 	 // }
						 	
						}
					}
				}
			}

			if (sai != -1 && entra != -1){
				//cout<<"entrou"<<endl;
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
				// aux->isTree();  // TODO: retirar
				retorno->adicionarSol(aux);

				//esvazia nonDominadedSolutions
				nonDominadedSolutions.clear();
				
				contAval++;
				 
				contttttt++;
				// contabiliza avaliaçoes da FO aqui	

			} else {
				cout<<"ERROOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO -1"<<endl;
				exit(-1);

			}
		}
		// cout<<"Soucoes vizitadas = "<<contttttt<<endl;
		// cout<<"Destino = "<<s2->getObj(0)<<" "<<s2->getObj(1)<<endl;
	

	}

public:
	Version8(){
		aux= new SolucaoEdgeSet(NUMEROVERTICES-1);

	}

};


#endif
