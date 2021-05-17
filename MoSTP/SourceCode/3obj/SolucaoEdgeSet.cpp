#ifndef SOLUCAOZG_CPP
#define SOLUCAOZG_CPP

/*This code file was kindly provided by Monteiro */


#include <cstdio>
#include <queue>
#include "Solucao.cpp"
#include "UnionFind.cpp"
#include "param_NSGAIII.h"
#include "auxEdgeStruct.h"
#include <algorithm>
#include <iostream>

using namespace std;

extern double custos[NUMOBJETIVOS][NUMEROVERTICES][NUMEROVERTICES];

typedef struct {
	int listaadj[NUMEROVERTICES][NUMEROVERTICES], graus[NUMEROVERTICES];
	bool completo;
} grafo;



class SolucaoEdgeSet : public Solucao {
	public:
	int edges[NUMEROVERTICES-1][2];
	int nEdges;
	UnionFind uf;
	grafo *g;
	double distance; // utilizada no NSGA-II (crownd distance) e na reciclagem do Hudson
	double fitness; // usado pelo SPEA2
	int strength; // S(i) do SPEA2
	int raw_fitness; // R(i) do SPEA2
	double antigof[NUMOBJETIVOS];
	int posicaoListaNSGAII; // guarda o index onde a soluçao é guardada na popupacao NUMPOPULACAO*2 do NSGA-II
	vector<double> distances; // utilizado para o SPEA2. Uma vez ordenado, pode-se escolher o kth vizinho mais proximo
	bool unionGraph[NUMEROVERTICES][NUMEROVERTICES];
	int vertlist [NUMEROVERTICES];
	bool visited;
	int rank;

	std::vector<double> getObj_NORMALIZADO;

	// // Copy constructor
	//   explicit SolucaoEdgeSet(const SolucaoEdgeSet & rhs)
	//   {
	//   	nEdges = rhs.nEdges;
	// 	f[0] = rhs.f[0];
	// 	f[1] = rhs.f[1];
	// 	f[2] = rhs.f[2];
	// 	memcpy(edges,rhs.edges,sizeof(edges));
	// 	antigof[0] = rhs.antigof[0];
	// 	antigof[1] = rhs.antigof[1];
	// 	antigof[2] = rhs.antigof[2];
	// 	distance = rhs.distance;
	// 	fitness = rhs.fitness;
	// 	strength = rhs.strength;
	// 	raw_fitness = rhs.raw_fitness;
	// 	distances = rhs.distances;
	// 	visited = rhs.visited;
	// 	rank = rhs.rank;
	// 	getObj_NORMALIZADO = rhs.getObj_NORMALIZADO;
	//   } // ctor

	SolucaoEdgeSet(int n):getObj_NORMALIZADO(NUMOBJETIVOS) {
		distance=0.0;
		fitness = 0.0;
		strength = 0;
		raw_fitness = 0;
		nEdges = n;
		f[0] = f[1] = f[2] = 0.0;
		g = NULL;
		visited = false;
		rank = -1;
	}
	~SolucaoEdgeSet() {

	}

    // copia somente os valores de fitness
    // util para guardar copias de solucoes somente para comparacao
    void shallowCopy(SolucaoEdgeSet &s) {
        s.f[0] = f[0];
        s.f[1] = f[1];
        s.f[2] = f[2];
    }

	/* Calcula o fitness atual da solucao
	 * Complexidade O(N) */
	void calcularFitness(int &contAval) {
		f[0] = f[1] = f[2] = 0.0;
		// contAvaliacoesObjetivoo++;
		for (int i=0;i<NUMEROVERTICES-1;i++)
			for (int j=0;j<NUMOBJETIVOS;j++)
				f[j] += f(j,edges[i][0],edges[i][1]);

		// cout<<f[0]<<" "<<f[1]<<" "<<f[2]<<endl;
		contAval++;
	}

	/* Calcula o fitness de uma soluçao parcial com tam arestas
	 * Complexidade O(N) */
	void calcularFitness(int tam, int &contAval) {
		// contAvaliacoesObjetivoo++;
		f[0] = f[1] = f[2] = 0.0;
		for (int i=0;i<tam;i++)
			for (int j=0;j<NUMOBJETIVOS;j++)
				f[j] += f(j,edges[i][0],edges[i][1]);

		// cout<<f[0]<<" "<<f[1]<<" "<<f[2]<<endl;
		contAval++;
	}



	/* gera uma arvore aleatoria
	 * Complexidade n */
    void doRandomWalk(int &contAval) {

       
        for (int i = 0; i < NUMEROVERTICES; i++)
            vertlist[i]=i;
        int vertctr = NUMEROVERTICES-1; //tamanho valido da lista;

        int v_ind = IRandom(0,vertctr);
        int v_esc = vertlist[v_ind];
        vertlist [v_ind] = vertlist[vertctr];
        vertlist [vertctr] = v_esc;
        vertctr--;
        int viz1, viz2, viz2_ind;

        int cont = 0;
        while (cont < NUMEROVERTICES-1) {

            viz1 = vertlist [IRandom(vertctr+1, NUMEROVERTICES-1)];
            viz2_ind = IRandom(0,vertctr);
            viz2 = vertlist[viz2_ind];

            vertlist [viz2_ind] = vertlist[vertctr];
            vertlist [vertctr] = viz2;
            vertctr--;

            edges[cont][0] = viz1;
            edges[cont][1] = viz2;
            cont++;
        }
        calcularFitness(contAval);
       // assert (confereArestas());
    }

    // para grafos completos e incompletos 
    void RandomWalk(){
    	uf.clear();
    	std::vector<pair<int, int> > amostral;
    	for (int i=0; i<NUMEROVERTICES; i++){
    		for (int j=i+1; j<NUMEROVERTICES; j++){
    			if(unionGraph[i][j]==true) {
    				amostral.push_back(make_pair(i,j));
    			}
    		}
    	}
    	int cont=0;
    	while (cont<NUMEROVERTICES-1 && amostral.size()>0){
			int are = IRandom(0,amostral.size()-1);

			if (uf.sameClass(amostral[are].first,amostral[are].second)==false){
				uf.unionClass(amostral[are].first,amostral[are].second);
				if (amostral[are].first<amostral[are].second){
					edges[cont][0] = amostral[are].first;
					edges[cont][1] = amostral[are].second;
				} else {
					edges[cont][0] = amostral[are].second;
					edges[cont][1] = amostral[are].first;
				}
				cont++;
			}
			amostral.erase(amostral.begin()+are);
		}
		if (cont<NUMEROVERTICES-1) cout<<"ERROR RandomWalk cont = "<<cont<<" amostral.size() = "<<amostral.size()<<endl;
    }

    /* Faz o crossover entre dois individuos.
    BAseado no crossover sugerido por Raidl and Julstrom (2003)*/
	void crossover(SolucaoEdgeSet *pai, SolucaoEdgeSet *mae, int &contAval) {
		
		//memset(unionGraph,false,sizeof(unionGraph));
		

		for (int i=0; i<NUMEROVERTICES; i++){
			for (int j=0; j<NUMEROVERTICES; j++){
				unionGraph[i][j] = false;
			}
		}
		for (int i=0;i<NUMEROVERTICES-1;i++) {
			if (!unionGraph[pai->edges[i][0]][pai->edges[i][1]]) {
				unionGraph[pai->edges[i][0]][pai->edges[i][1]] = true;
				unionGraph[pai->edges[i][1]][pai->edges[i][0]] = true;
			}
		}
		for (int i=0;i<NUMEROVERTICES-1;i++) {
			if (!unionGraph[mae->edges[i][0]][mae->edges[i][1]]) {
				unionGraph[mae->edges[i][0]][mae->edges[i][1]] = true;
				unionGraph[mae->edges[i][1]][mae->edges[i][0]] = true;
			}
		}
		RandomWalk();
		calcularFitness(contAval);
		// exit(1);
	}


	/* Gera um individuo aleatorio
	 * Complexidade O(n lg n) */
	void geraIndividuoAleatorio(int &contAval) {
	    g = new grafo;
		g->completo = true;

		for (int i=0;i<NUMEROVERTICES;i++)
				g->graus[i] = NUMEROVERTICES;
		doRandomWalk(contAval);
		delete g;
	}

	

	/* Faz a troca das arestas ai e aj, religando no formato 2-OPT
	 * Complexidade O(1) */
	void trocaArestas(int ai, int aj, char tipo, SolucaoEdgeSet &soloriginal) {
		int a = soloriginal.edges[ai][0];
		int b = soloriginal.edges[ai][1];
		int c = soloriginal.edges[aj][0];
		int d = soloriginal.edges[aj][1];

		//inicialmente temos a-b e c-d

		int novaA[2]={0,0}, novaB[2]={0,0};
		novaA[0] = a;
		novaB[0] = b;

		// se todos ficarem na mesma componente, junta AC/BD, se nao junta AD/BC
		if (tipo == 0) {
			// faz a uniao AC e BD
			novaA[1] = c;
			novaB[1] = d;
		} else if (tipo == 1){
			// faz a uniao AD e BC
			novaA[1] = d;
			novaB[1] = c;
		}
		if (novaA[0]<novaA[1]){
			edges[ai][0] = novaA[0];
			edges[ai][1] = novaA[1];
		} else {
			edges[ai][0] = novaA[1];
			edges[ai][1] = novaA[0];
		}
		if (novaB[0]<novaB[1]){
			edges[aj][0] = novaB[0];
			edges[aj][1] = novaB[1];
		} else {
			edges[aj][0] = novaB[1];
			edges[aj][1] = novaB[0];
		}
	}

	char calcularTrocaArestas(int ai, int aj, SolucaoEdgeSet &soloriginal, int &contAval) {
		// Calcula a possivel troca das arestas "ai" com "aj"
		// Complexidade O(n)
		// uniao busca
		contAval++;
		uf.clear();
		for (int k=0;k<NUMOBJETIVOS;k++) f[k] = 0.0; // (re)inicializa os objetivos
		
		for (int i=0;i<NUMEROVERTICES-1;i++)
			if (i != ai && i != aj) {
				edges[i][0] = soloriginal.edges[i][0];
				edges[i][1] = soloriginal.edges[i][1];
				uf.unionClass(edges[i][0],edges[i][1]);
				for (int k=0;k<NUMOBJETIVOS;k++)
					f[k] +=  f(k,edges[i][0],edges[i][1]);
			}
		int a = soloriginal.edges[ai][0];
		int b = soloriginal.edges[ai][1];
		int c = soloriginal.edges[aj][0];
		int d = soloriginal.edges[aj][1];

		char tipoTroca;
		// junta a-c e b-d
		uf.unionClass(a,c);
		uf.unionClass(b,d);
		// se todos ficarem na mesma componente, junta AC/BD, se nao junta AD/BC
		if (uf.sameClass(a,b) && uf.sameClass(a,c) && uf.sameClass(a,d)) {
			// faz a uniao AC e BD
			tipoTroca = 0;
			for (int k=0;k<NUMOBJETIVOS;k++)
				f[k] = f[k] + ( f(k,a,c) + f(k,b,d) );
		} else {
			// faz a uniao AD e BC
			tipoTroca = 1;
			for (int k=0;k<NUMOBJETIVOS;k++)
				f[k] = f[k] + ( f(k,a,d) + f(k,b,c) );
		}
		return tipoTroca;
	}

	//efetua remocao de aresta ai, adicionando a aresta de menor custo que reconecta a arvore
	void substAresta (int ai, int index, fitVecNode * fitVec[], int &contAval) {

		contAval++;
		uf.clear();

		antigof[0] = getObj(0);
		antigof[1] = getObj(1);
		antigof[2] = getObj(2);
		f[0] = f[1] = f[2] = 0.0;

		int a = edges[ai][0];
		int b = edges[ai][1];
		for (int i=0;i<NUMEROVERTICES-1;i++) {
			if (i != ai) {
				uf.unionClass(edges[i][0],edges[i][1]);
				for (int k=0;k<NUMOBJETIVOS;k++)
					f[k] +=  f(k,edges[i][0],edges[i][1]);
			}
		}

		int edge = 0;
		// anda ate a proxima aresta que pode ser inserida

		while (uf.sameClass(fitVec[index][edge].a,fitVec[index][edge].b)) edge++;

		if ( (fitVec[index][edge].a == a && fitVec[index][edge].b == b) ||
				(fitVec[index][edge].a == b && fitVec[index][edge].b == a) ) {
			edge++;
			while (uf.sameClass(fitVec[index][edge].a,fitVec[index][edge].b)) edge++;
		}
		// coloca a aresta na solucao
		edges[ai][0] = fitVec[index][edge].a;
		edges[ai][1] = fitVec[index][edge].b;
		setObj(0,getObj(0)+f(0,edges[ai][0],edges[ai][1]));
		setObj(1,getObj(1)+f(1,edges[ai][0],edges[ai][1]));
		setObj(2,getObj(2)+f(2,edges[ai][0],edges[ai][1]));
		uf.unionClass( fitVec[index][edge].a, fitVec[index][edge].b );

	//	assert (confereArestas());
	}

	void desfazerCalculoTroca() {
		f[0] = antigof[0];
		f[1] = antigof[1];
		f[2] = antigof[2];
	}

	// TEM DE SER CHAMADO ANTES DA TROCA!
	int getArestaEntrou(int i, int j, char tipo, int arest, int no) {
		if (tipo == 0) {
			int nos[2][2] = { {edges[i][0],edges[j][0]}, {edges[i][1],edges[j][1]} };
			return nos[ arest ][ no ];
		} else {
			int nos[2][2] = { {edges[i][0],edges[j][1]}, {edges[i][1],edges[j][0]} };
			return nos[ arest ][ no ];
		}
	}
	int getArestaSaiu(int i, int j, char tipo, int arest, int no) {
		int nos[2][2] = { {edges[i][0],edges[i][1]}, {edges[j][0],edges[j][1]} };
		return nos[ arest ][ no ];
	}


	void mutacao(SolucaoEdgeSet &sol, int &contAval){
		// contMutacoes++; // estatistica
		int a1 = IRandom(0,NUMEROVERTICES-1-1), a2;
		while ((a2 = IRandom(0,NUMEROVERTICES-1-1)) == a1);
		trocaArestas(a1,a2,calcularTrocaArestas(a1,a2,sol,contAval),sol);
	}

	void operator=(SolucaoEdgeSet &s) {
		nEdges = s.nEdges;
		f[0] = s.f[0];
		f[1] = s.f[1];
		f[2] = s.f[2];
		memcpy(edges,s.edges,sizeof(edges));
		antigof[0] = s.antigof[0];
		antigof[1] = s.antigof[1];
		antigof[2] = s.antigof[2];
		distance = s.distance;
		fitness = s.fitness;
		strength = s.strength;
		raw_fitness = s.raw_fitness;
		distances = s.distances;
		visited = s.visited;
		rank = s.rank;
		getObj_NORMALIZADO = s.getObj_NORMALIZADO;
	}

	void printSolucao(FILE *f) {
		fprintf(f,"Custos = (%.3lf,%.3lf,%.3lf)\n",this->f[0],this->f[1], this->f[2]);
		fprintf(f,"Arestas = \n");
		for (int i=0;i<NUMEROVERTICES-1;i++)
			fprintf(f,"(%d - %d)\n",edges[i][0],edges[i][1]);
		fprintf(f,"\n");
	}

	void printPonto(FILE *f) {
		fprintf(f,"%.6lf %.6lf %.6lf\n",getObj(0),getObj(1), getObj(2));
	}

	bool confereArestas() {
		uf.clear();
		bool usados[NUMEROVERTICES] = {false};
		int numUsados = 0;
		for (int i=0;i<NUMEROVERTICES-1;i++) {
			if (uf.sameClass(edges[i][0],edges[i][1])) return false;
			uf.unionClass(edges[i][0],edges[i][1]);
			for (int j=0;j<2;j++) {
		                if (!usados[edges[i][j]]) {
	               			usados[edges[i][j]] = true;
        			        numUsados++;
                		}
			}
		}
		if (numUsados != NUMEROVERTICES) return false;
		return true;
	}

	bool confereObjetivos() {
		double v[3] = {0.0,0.0,0.0};
		for (int i=0;i<NUMEROVERTICES-1;i++) {
			v[0] += f(0,edges[i][0],edges[i][1]);
			v[1] += f(1,edges[i][0],edges[i][1]);
			v[2] += f(2,edges[i][0],edges[i][1]);
		}
		if (fabs(v[0])-EPS > f[0] || fabs(v[1])-EPS > f[1] || fabs(v[2])-EPS > f[2]) return false;
		return true;
	}




	bool isTree(){ // verificador
		UnionFind uf;
		bool vetkktk[NUMEROVERTICES];
		for(int i=0; i<NUMEROVERTICES; i++)vetkktk[i] = false;
		for (int i=0; i<NUMEROVERTICES-1; i++){
			if (uf.sameClass(this->edges[i][0],this->edges[i][1])==false){
				vetkktk[this->edges[i][0]] = true;
				vetkktk[this->edges[i][1]] = true;
				uf.unionClass(this->edges[i][0],this->edges[i][1]);
			} else {
				cout<<"ERROROROROROROROROROROROROROR TREEEE"<<endl;
				return false;
			}
		}
		for(int i=0; i<NUMEROVERTICES; i++){
			if (vetkktk[i]==false) {
				cout<<"ERROROROROROROROROROROROROROR TREEEE"<<endl;
				
				return false;
			}
		}
		if (confereObjetivos()==false) cout<<"OBJECTIVES PROBLEMS"<<endl;
		cout<<"It's a tree! õ/"<<endl;;
		return true;
	}

};

#endif
