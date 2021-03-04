#ifndef RMCKRUSKAL_CPP
#define RMCKRUSKAL_CPP

#include "param_NSGAII.h"
#include "UnionFind.cpp"
#include "SolucaoEdgeSet.cpp"
#include <vector>
#include <algorithm>

#define INF 1e9
#define for(i,n) for (int i=0;i<n;i++)

typedef struct {
	int a, b;
	double c;
} edge;

inline bool compEdge(edge a, edge b) {
	if (a.c < b.c) return true;
	return false;
}


class rmcKruskal {
private:
	edge arestas[(NUMEROVERTICES*(NUMEROVERTICES-1))/2];
	UnionFind uf;
	int nArestas;

public:
	rmcKruskal() {
		nArestas = 0;
		for(i,NUMEROVERTICES) {
			for (j,i) {
				arestas[nArestas].a = i;
				arestas[nArestas].b = j;
				nArestas++;
			}
		}
	}

	void executar(SolucaoEdgeSet &s, double lambda[NUMOBJETIVOS], double alfa, int &contAval) {
		for(i,NUMOBJETIVOS) {
			s.setObj(i,0.0);
		}
		contAval++;
		uf.clear();
		for(i,nArestas)
			arestas[i].c = lambda[0]*f(0,arestas[i].a,arestas[i].b) + lambda[1]*f(1,arestas[i].a,arestas[i].b);

		// ordena as arestas do grafo
		stable_sort(arestas,arestas+((NUMEROVERTICES*(NUMEROVERTICES-1))/2),compEdge);
		
		// coloca NUMEROVERTICES-1 arestas do grafo sem formar ciclo
		int cont = 0, edge = 0;
		while (cont < NUMEROVERTICES-1) {

			// anda ate a proxima aresta que pode ser inserida
			while (uf.sameClass(arestas[edge].a,arestas[edge].b)) edge++;

			// encontra as arestas que tem custo ate alfa% maior que o da aresta "edge"
			int edgeRCL = edge+1;
			double maxCost = arestas[edge].c*(1.0+alfa);
			while (edgeRCL < nArestas && arestas[edgeRCL].c < maxCost)
				edgeRCL++;

			// escolhe aleatoriamente uma aresta da lista ate que encontre uma que nao forma ciclo
			int esc = IRandom(edge,edgeRCL-1);
			while (uf.sameClass(arestas[esc].a,arestas[esc].b))
				esc = IRandom(edge,edgeRCL-1);
	
			//printf("%d: aresta %d (%d,%d) --> (%lf)\n",cont,esc,arestas[esc].a,arestas[esc].b,custoKruskal[arestas[esc].a][arestas[esc].b]);

			// coloca a aresta escolhida da RCL na solucao
			s.edges[cont][0] = arestas[esc].a;
			s.edges[cont][1] = arestas[esc].b;
			s.setObj(0,s.getObj(0)+f(0,s.edges[cont][0],s.edges[cont][1]));
			s.setObj(1,s.getObj(1)+f(1,s.edges[cont][0],s.edges[cont][1]));
			uf.unionClass( arestas[esc].a, arestas[esc].b );
			cont++;
		}
	}
};
#undef INF
#undef for

#endif
