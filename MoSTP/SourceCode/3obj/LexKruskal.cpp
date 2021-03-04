#ifndef LEXKRUSKAL_CPP
#define LEXKRUSKAL_CPP


#include "param_NSGAIII.h"
#include "SolucaoEdgeSet.cpp"
#include "auxEdgeStruct.h"
#include "UnionFind.cpp"
#include <iostream>
using namespace std;

class LexKruskal {
private:
	UnionFind uf;

public:
	LexKruskal() {}

	void executar(SolucaoEdgeSet &s, std::vector<auxEdgeStruct> arestas, int &contAval) {
		contAval++;
		for(int i = 0; i < NUMOBJETIVOS; i++) {
			s.setObj(i,0.0);
		}

		uf.clear();
		
		// coloca NUMEROVERTICES-1 arestas do grafo sem formar ciclo
		int cont = 0, edge = 0;
		while (cont < NUMEROVERTICES-1) {

			// anda ate a proxima aresta que pode ser inserida
			while (uf.sameClass(arestas[edge].a,arestas[edge].b)) edge++;

			// coloca a aresta na solucao
			s.edges[cont][0] = arestas[edge].a;
			s.edges[cont][1] = arestas[edge].b;
			s.setObj(0,s.getObj(0)+f(0,s.edges[cont][0],s.edges[cont][1]));
			s.setObj(1,s.getObj(1)+f(1,s.edges[cont][0],s.edges[cont][1]));
			s.setObj(2,s.getObj(2)+f(2,s.edges[cont][0],s.edges[cont][1]));
			uf.unionClass( arestas[edge].a, arestas[edge].b );
			cont++;
		}
		// assert (s.confereArestas());
	}
};

#endif

