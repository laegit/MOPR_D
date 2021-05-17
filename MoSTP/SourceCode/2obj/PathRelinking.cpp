
#ifndef PATHRELINKING_CPP
#define PATHRELINKING_CPP

#include "SolucaoEdgeSet.cpp"
#include "BoundedParetoSet.cpp"
#include "param_NSGAII.h"


typedef struct aa {

	double obj1;
	double obj2;
	int sai; // index em amostralSai
	int entra; // index em amostralEnt

} auxSoluacao;

class PathRelinking{
	public:
		virtual void PR(BoundedParetoSet *retorno, SolucaoEdgeSet *s1, SolucaoEdgeSet *s2, int &contAval) = 0; 

		// guarda arestas de s1 que nao estao em s2. 
		int inS1_NotInS2(SolucaoEdgeSet *s1, SolucaoEdgeSet *s2, int amostral[NUMEROVERTICES-1]){ // O(n^2)
			int cont=0;
			for (int i=0; i<NUMEROVERTICES-1; i++){
				bool ha = false;
				for (int j=0; j<NUMEROVERTICES-1 && ha==false; j++){
					if ((s2->edges[j][0]==s1->edges[i][0] && s2->edges[j][1]==s1->edges[i][1]) || (s2->edges[j][0]==s1->edges[i][1] && s2->edges[j][1]==s1->edges[i][0])){
						ha = true;
					}
				}	
				if (ha==false) {
					amostral[cont] = i; // aresta i de s1 deve sair
					cont++;
				}
			}
			return cont;
		}
};

#endif
