#ifndef VERSION6_CPP
#define VERSION6_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2018)
	
This code implements the 6th PathRelinking version for Bi-objective spanning tree 
The data structure and some functions were kindly provided by Monteiro (2010)

Pareto-based

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


class Version6: public PathRelinking{
private:
	int amostralSai[NUMEROVERTICES-1]; // guarda indexes de arestas (de s1) que devem sair de s1 (vide abaixo)
	int amostralEnt[NUMEROVERTICES-1]; // guarda indexes de arestas (de s2) que entrar em s1 (vide abaixo)

	SolucaoEdgeSet *aux;
	std::list<auxSoluacao > nonDominadedSolutions;



	bool A_domina_B(double a_obj1, double a_obj2, double b_obj1, double b_obj2){
		return (a_obj1<=b_obj1 && a_obj2<=b_obj2 && (a_obj1<b_obj1 || a_obj2<b_obj2)); //minimizaçao
	}

	bool addNonDominandedSolution(double obj1, double obj2, int sai,  int entra){
		list<auxSoluacao>::iterator i = nonDominadedSolutions.begin();
		list< list<auxSoluacao>::iterator > remover;
		while (i != nonDominadedSolutions.end()) {
			// if (*nova >> *((*i).sol))
			if (A_domina_B(obj1,obj2, (*i).obj1, (*i).obj2))
				remover.push_back(i);
			// if (*((*i).sol) >> *nova)
			if (A_domina_B((*i).obj1, (*i).obj2, obj1,obj2) || ((*i).obj1==obj1 && (*i).obj2==obj2))
				return false;
			i++;
		}

		list< list<auxSoluacao>::iterator >::iterator j = remover.begin();
		while (j != remover.end()) {
			nonDominadedSolutions.erase( *j );
			j++;
		}
		nonDominadedSolutions.push_front(auxSoluacao {obj1, obj2, sai, entra});
		return true;
	}

	void PR(BoundedParetoSet *retorno, SolucaoEdgeSet *s1, SolucaoEdgeSet *s2, int &contAval){
		
		int contttttt = 0;
		// O(n^2)
		int totalSai = inS1_NotInS2(s1, s2, amostralSai); // arestas de s1 que nao estao em s2. Tais arestas devem sair de s1
		int totalEnt = inS1_NotInS2(s2, s1, amostralEnt); // arestas de s2 que nao estao em s1. Tais arestas devem entrar em s1

		*aux = *s1; // O(n)

		nonDominadedSolutions.clear();
		
		// cout<<"Origem = "<<s1->getObj(0)<<" "<<s1->getObj(1)<<endl;

		// Total: O(n^3logn)
		for (int cont = 0; cont<totalSai && contAval<NUM_AVALIACOES; cont++){
			// int sai = -1;
			// int entra = -1;
			
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
							
						 	 addNonDominandedSolution(newObj_0,newObj_1, iii, jd);
						 	 // add ao nonDominadedSolutions

						 	
						}
					}
				}
			}

			if (nonDominadedSolutions.size()>0){
				
				// cout<<"nonDominadedSolutions.size() = "<<nonDominadedSolutions.size()<<endl;
				int randdd = IRandom(0,nonDominadedSolutions.size()-1);
				list<auxSoluacao>::iterator it = nonDominadedSolutions.begin();
				for (int i = 0; i < randdd; ++i) it++;


				int i_sai = amostralSai[(*it).sai]; // aresta i_sai que sai
				int o_sai = aux->edges[i_sai][0];
				int d_sai = aux->edges[i_sai][1];
				amostralSai[(*it).sai] = -1;


				int j_entra = amostralEnt[(*it).entra]; // aresta j_entra que entra
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
				// cout<<aux->getObj(0)<<" "<<aux->getObj(1)<<endl;
				contttttt++;
				// contabiliza avaliaçoes da FO aqui	

			} else {
				cout<<"ERROOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO -1"<<endl;
				exit(-1);

			}
		}
	}

public:
	Version6(){
		
		aux= new SolucaoEdgeSet(NUMEROVERTICES-1);

	}

};


#endif
