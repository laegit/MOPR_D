#ifndef VERSION3_CPP
#define VERSION3_CPP

/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2018)
	
This code implements the third PathRelinking version for Bi-objective spanning tree 
The data structure and some functions were kindly provided by Monteiro (2010)

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



class Version3: public PathRelinking{
private:
	int amostralSai[NUMEROVERTICES-1]; // guarda indexes de arestas (de s1) que devem sair de s1 (vide abaixo)
	int amostralEnt[NUMEROVERTICES-1]; // guarda indexes de arestas (de s2) que entrar em s1 (vide abaixo)

	SolucaoEdgeSet *aux;

	void PR(BoundedParetoSet *retorno, SolucaoEdgeSet *s1, SolucaoEdgeSet *s2, int &contAval){
		
		int obj = 0;
		int contttttt = 0;
		// O(n^2)
		int totalSai = inS1_NotInS2(s1, s2, amostralSai); // arestas de s1 que nao estao em s2. Tais arestas devem sair de s1
		int totalEnt = inS1_NotInS2(s2, s1, amostralEnt); // arestas de s2 que nao estao em s1. Tais arestas devem entrar em s1

		*aux = *s1; // O(n)
		// retorno->adicionarSol(aux); // depreciado (cópia prévia)

		// Total: O(n^3logn)
		for (int cont = 0; cont<totalSai && contAval<NUM_AVALIACOES; cont++){
			// cout<<"OBJ = "<<obj<<endl;
			double mimObj = 1e9;
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
							 double newObj = aux->getObj(obj)-custos[obj][origem][destino]+custos[obj][a][b];
							 if (newObj<mimObj){
							 	mimObj = newObj;
							 	entra = jd;
							 	sai = iii;
							 }
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
				// aux->isTree();  // TODO: retirar
				retorno->adicionarSol(aux);
				// cout<<"retorno->getSize()_v3 = "<<retorno->getSize()<<endl;
				
				contAval++;

				contttttt++;
				// contabiliza avaliaçoes da FO aqui	

			} else {
				cout<<"ERROOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO -1"<<endl;
				exit(-1);

			}
		
			obj++;
			if (obj == NUMOBJETIVOS) obj = 0;
		}

		// cout<<"Soucoes vizitadas = "<<contttttt<<endl;
	}

public:
	Version3(){
		aux= new SolucaoEdgeSet(NUMEROVERTICES-1);

	}

};


#endif
