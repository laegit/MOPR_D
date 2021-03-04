#ifndef RMC_PRIM
#define RMC_PRIM

#include "param_NSGAIII.h"
#include "SolucaoEdgeSet.cpp"

#include <iostream>
using namespace std;

/*This code file was kindly provided by Monteiro */


#define INF 1e9
#define PER_ACCEPTED 0.03 
#define for(i,n) for (int i=0;i<n;i++)


/*Procedimento rmc-prim normal ==> que gera uma árvore total 
populacao inicial do NSGAII */
void rmcPrim (SolucaoEdgeSet &s, double lambda[NUMOBJETIVOS], int &contAval) {
	//cout<<"============== rm1 ============== "<<endl;
	bool arvore[NUMEROVERTICES]; // elementos ja colocados na arvore
	double mindist[NUMEROVERTICES]; // menor distancia ate cada no
	int cont = 0; // numero de elementos ja colocados
	int viz, ant[NUMEROVERTICES];
	double maisprox;
	double maisdist;
	int lrc [NUMEROVERTICES];
	int lrctam = 0;
	int inicio = IRandom(0,NUMEROVERTICES-1);

	// inicialmente, nenhum elemento esta na arvore e
	// a menor distancia eh infinito
	for(i,NUMEROVERTICES) {
		arvore[i] = false;
		mindist[i] = INF;
	}

	// menor distancia ate o inicio eh zero
	mindist[inicio] = 0.0;
	viz = inicio;

	while (cont < NUMEROVERTICES-1) {
		// coloca o no na arvore
		arvore[viz] = true;
		// ajusta a distancia dos vizinhos do no
		for(i,NUMEROVERTICES) {
			double p = 0.0;
			for(j,NUMOBJETIVOS) p += (lambda[j]*custos[j][viz][i]);
			if (!arvore[i] && i != viz && mindist[i] > p) {
				mindist[i] = p;
				ant[i] = viz;
			}
		}

		viz = -1;
		maisprox = INF;
		maisdist = 0.0;
		lrctam = 0;

		// coloca os que nao tem distancia inf na lrc. o tamanho da lrc sera regulado dps de acordo com o custo
		for(i,NUMEROVERTICES) {
			if (!arvore[i] && mindist[i] < INF) {
				lrc [lrctam] = i;
				lrctam++;
				if (mindist[i] < maisprox)
					maisprox = mindist[i];
				if (mindist[i] > maisdist)
					maisdist = mindist[i];
			}
		}

		double maxcustoaceito = maisprox + (maisdist-maisprox)*PER_ACCEPTED;

		int auxlrc [NUMEROVERTICES];
		int auxlrctam = 0;
		for(i,lrctam) {
			if (mindist[lrc[i]] > maxcustoaceito) continue;
			else {
				auxlrc [auxlrctam] = lrc[i];
				auxlrctam++;
			}
		}

		viz = auxlrc [IRandom (0, auxlrctam-1)];

		// adiciona a aresta entre ant[viz] e viz
		if ( ant[viz]< viz){
			s.edges[cont][0] = ant[viz];
			s.edges[cont][1] = viz;
		} else {
			s.edges[cont][0] = viz;
			s.edges[cont][1] = ant[viz];
		}		

		cont++;
	}
	// atualiza os custos da solucao
	s.calcularFitness(contAval);
}

void rmcPrim (SolucaoEdgeSet &s, double lambda[NUMOBJETIVOS], int tamanho, int &contAval) {
	//cout<<"============== rm2 ============== "<<endl;
	
	bool arvore[NUMEROVERTICES]; // elementos ja colocados na arvore
	double mindist[NUMEROVERTICES]; // menor distancia ate cada no
	int cont = 0; // numero de elementos ja colocados
	int viz, ant[NUMEROVERTICES];
	double maisprox;
	double maisdist;
	int lrc [NUMEROVERTICES];
	int lrctam = 0;

	int inicio = IRandom(0,NUMEROVERTICES-1); // vértice inicial

	// inicialmente, nenhum elemento esta na arvore e
	// a menor distancia eh infinito
	for(i,NUMEROVERTICES) {
		arvore[i] = false;
		mindist[i] = INF;
	}

	// menor distancia ate o inicio eh zero
	mindist[inicio] = 0.0;
	viz = inicio;
	while (cont < tamanho) {
		// coloca o no na arvore
		arvore[viz] = true;
		// ajusta a distancia dos vizinhos do no
		for(i,NUMEROVERTICES) {
			double p = 0.0;
			for(j,NUMOBJETIVOS) p += (lambda[j]*custos[j][viz][i]);
			if (!arvore[i] && i != viz && mindist[i] > p) {
				mindist[i] = p;
				ant[i] = viz;
			}
		}

		viz = -1;
		maisprox = INF;
		maisdist = 0.0;
		lrctam = 0;
		
		// coloca os que nao tem distancia inf na lrc. o tamanho da lrc sera regulado dps de acordo com o custo
		for(i,NUMEROVERTICES) {
			if (!arvore[i] && mindist[i] < INF) {
				lrc [lrctam] = i;
				lrctam++;
				if (mindist[i] < maisprox)
					maisprox = mindist[i];
				if (mindist[i] > maisdist)
					maisdist = mindist[i];
			}
		}
		double maxcustoaceito = maisprox + (maisdist-maisprox)*PER_ACCEPTED;

		int auxlrc [NUMEROVERTICES];
		int auxlrctam = 0;
		for(i,lrctam) {
			if (mindist[lrc[i]] > maxcustoaceito) continue;
			else {
				auxlrc [auxlrctam] = lrc[i];
				auxlrctam++;
			}
		}

		viz = auxlrc [IRandom (0, auxlrctam-1)];

		// adiciona a aresta entre ant[viz] e viz
		
		if ( ant[viz]< viz){
			s.edges[cont][0] = ant[viz];
			s.edges[cont][1] = viz;
		} else {
			s.edges[cont][0] = viz;
			s.edges[cont][1] = ant[viz];
		}
		cont++;
	}
	// atualiza os custos da solucao
	s.calcularFitness(contAval);
}

// usado num transponson => gera uma arvore "s" que possui as arestas de "ref".
// @param reftam é a quantidade de arestas de ref
void mcprim(SolucaoEdgeSet &ref, int reftam, SolucaoEdgeSet &s, double lambda[NUMOBJETIVOS], int &contAval) {

	bool arvore[NUMEROVERTICES]; // elementos ja colocados na arvore
	double mindist[NUMEROVERTICES]; // menor distancia ate cada no
	int cont = 0; // numero de elementos ja colocados
	int viz, ant[NUMEROVERTICES];
	double maisprox;

	// inicialmente, nenhum elemento esta na arvore e
	// a menor distancia eh infinito
	for(i,NUMEROVERTICES) {
		arvore[i] = false;
		mindist[i] = INF;
	}

	// a seguir, todos os elementos de ref sao inseridos na arvore
	
	for (i, 2) {
		arvore [ref.edges[0][i]] = true;
		//atualizacao das distancias minimas dos nos adjacentes
		for(k,NUMEROVERTICES) {
			double p = 0.0;
			for(l,NUMOBJETIVOS) p += (lambda[l]*custos[l][ ref.edges[0][i] ][k]);
			if (!arvore[k] && k != ref.edges[0][i] && mindist[k] > p) {
				mindist[k] = p;
				ant[k] = ref.edges[0][i];
			}
		}
	}
	s.edges[cont][0] = ref.edges[0][0];
	s.edges[cont][1] = ref.edges[0][1];
	cont++;
	

	int ctrIn = 0;
	bool newNode;

	do {
		for (i, reftam) {
			newNode = false;
			ctrIn = 0;

			if (!arvore [ ref.edges[i][1] ] && arvore [ ref.edges [i][0] ] ) {
				arvore [ref.edges[i][1]] = true;
				//atualizacao das distancias minimas dos nos adjacentes
				for(k,NUMEROVERTICES) {
					double p = 0.0;
					for(l,NUMOBJETIVOS) p += (lambda[l]*custos[l][ ref.edges[i][1] ][k]);
					if (!arvore[k] && k != ref.edges[i][1] && mindist[k] > p) {
						mindist[k] = p;
						ant[k] = ref.edges[i][1];
					}
				}
				s.edges[cont][0] = ref.edges[i][0];
				s.edges[cont][1] = ref.edges[i][1];
				cont++;
				ctrIn++;
			}

			else if (arvore [ ref.edges[i][1] ] && !arvore [ ref.edges [i][0] ] ) {
				arvore [ref.edges[i][0]] = true;
				//atualizacao das distancias minimas dos nos adjacentes
				for(k,NUMEROVERTICES) {
					double p = 0.0;
					for(l,NUMOBJETIVOS) p += (lambda[l]*custos[l][ ref.edges[i][0] ][k]);
					if (!arvore[k] && k != ref.edges[i][0] && mindist[k] > p) {
						mindist[k] = p;
						ant[k] = ref.edges[i][0];
					}
				}
				s.edges[cont][0] = ref.edges[i][0];
				s.edges[cont][1] = ref.edges[i][1];
				cont++;
				ctrIn++;
			}
		}
	}while (ctrIn > 0);
	
	viz = -1;
	maisprox = INF;

	// procura o no que tem a menor distancia
	// em relacao a arvore formada pelas arestas em ref
	// e que ainda nao foi colocado
	for(i,NUMEROVERTICES) {
		if (!arvore[i] && mindist[i] < maisprox) {
			viz = i;
			maisprox = mindist[i];
		}
	}


	// adiciona a aresta entre ant[viz] e viz
	s.edges[cont][0] = ant[viz];
	s.edges[cont][1] = viz;
	cont++;

	while (cont < NUMEROVERTICES-1) {
		// coloca o no na arvore
		arvore[viz] = true;
		// ajusta a distancia dos vizinhos do no
		for(i,NUMEROVERTICES) {
			double p = 0.0;
			for(j,NUMOBJETIVOS) p += (lambda[j]*custos[j][viz][i]);
			if (!arvore[i] && i != viz && mindist[i] > p) {
				mindist[i] = p;
				ant[i] = viz;
			}
		}
		viz = -1;
		maisprox = INF;

		// procura o no que tem a menor distancia nesta iteracao
		// e que ainda nao foi colocado
		for(i,NUMEROVERTICES) {
			if (!arvore[i] && mindist[i] < maisprox) {
				viz = i;
				maisprox = mindist[i];
			}
		}
		if (viz == -1) {
			fprintf (stderr, "ERRO VIZ!!!");
			//sleep (1);
		}

		// adiciona a aresta entre ant[viz] e viz
		s.edges[cont][0] = ant[viz];
		s.edges[cont][1] = viz;
		cont++;
	}

	// atualiza os custos da solucao
	s.calcularFitness(contAval);

}

#undef INF
#undef for

#endif
