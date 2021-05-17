#ifndef PARETOSET_CPP
#define PARETOSET_CPP

/*This code file was kindly provided by Monteiro */


#include <list>
#include <map>
#include "SolucaoEdgeSet.cpp"
#include "param_NSGAIII.h"
#include "Grid.cpp"

using namespace std;

typedef struct {
	double min, max;
} range;

class ParetoSet {
	public: // TODO: deveria ser protected e nao public
	double hyperolume;
	list<SolucaoEdgeSet *> sol;
	range rangeNovo[NUMOBJETIVOS], rangeAtual[NUMOBJETIVOS];
	int frequencia[NUMEROVERTICES][NUMEROVERTICES];
	Grid g;

	int calcularGridPos(SolucaoEdgeSet &s) {
		int bit = 0;
		int gridPos = 0;
		for (int obj=0;obj<NUMOBJETIVOS;obj++) {
			double inicio = rangeAtual[obj].min, fim = rangeAtual[obj].max, meio = (inicio+fim)/2.0;
			for (int k=0;k<PROFUNDIDADEGRID;k++) {
				if (s.getObj(obj) >= meio) {
					gridPos |= (1 << bit);
					inicio = meio;
				} else {
					fim = meio;
				}
				meio = (inicio+fim)/2.0;
				bit++;
			}
		}

		return gridPos;
	}
	void updateGrid() {
		g.clearGrid();

		list<SolucaoEdgeSet *>::iterator it = sol.begin();
		reiniciarRanges();
		while (it != sol.end()) {
			for (int k=0;k<NUMOBJETIVOS;k++) {
				rangeAtual[k].min = rangeNovo[k].min = min(rangeAtual[k].min,(*it)->getObj(k));
				rangeAtual[k].max = rangeNovo[k].max = max(rangeAtual[k].max,(*it)->getObj(k));
			}
			it++;
		}

		it = sol.begin();
		while (it != sol.end()) {
			g.addGrid( calcularGridPos(**it) );
			it++;
		}
	}

	void reiniciarRanges() {
		#define INF 1e9
		for (int k=0;k<NUMOBJETIVOS;k++) {
			rangeAtual[k].min = rangeNovo[k].min = INF;
			rangeAtual[k].max = rangeNovo[k].max = -INF;
		}
		#undef INF
	}

    /* FUNCOES DE FREQUENCIA DE ARESTAS NO PARETO SET */
    void decrementarFrequencia( int a, int b ) {
        if (a > b) std::swap(a,b);
      //  ASS( assert (frequencia[a][b] != 0); )
        frequencia[a][b]--;
    }

    void incrementarFrequencia( int a, int b ) {
        if (a > b) std::swap(a,b);
        frequencia[a][b]++;
    }

	void removerFrequencia( SolucaoEdgeSet *s ) {
	    for (int i=0;i < NUMEROARESTAS;i++)
            decrementarFrequencia( s->edges[i][0], s->edges[i][1] );
	}

	void adicionarFrequencia( SolucaoEdgeSet *s ) {
	    for (int i=0;i < NUMEROARESTAS;i++)
	        incrementarFrequencia ( s->edges[i][0] , s->edges[i][1] );
	}
	/* FIM DAS FUNCOES DE FREQUENCIA NO PARETO SET */

	public:

	double getHypervolume(){
		return hyperolume;
	}

	void setHypervolume(double hyp){
		hyperolume = hyp;
	}
	ParetoSet() {
		hyperolume = 1e9;
		reiniciarRanges();
		memset(frequencia,0,sizeof(frequencia));
	}
	virtual ~ParetoSet() {
		clear();
	}

	int getPositionCount(SolucaoEdgeSet &s) {
		return g.getPositionCount( calcularGridPos(s) );
	}

	int getPositionCount(int p) {
		return g.getPositionCount( p );
	}

	list<SolucaoEdgeSet *> getElementos() {
		return sol;
	}

	/* Complexidade: O(n) */
	virtual bool adicionarSol(SolucaoEdgeSet *s) {
		/* percorre o vetor de solucoes e de valores e, caso exista solucao dominada, retira e retorna true. caso contrario, retorna false */
		list<SolucaoEdgeSet *>::iterator i = sol.begin();
		list< list<SolucaoEdgeSet *>::iterator > remover;
		while (i != sol.end()) {
			if (*s >> **i) { // se *s domina *i
				remover.push_back(i);
				//printf("Dominada -> (%.3lf,%.3lf,%.3lf) por (%.3lf,%.3lf,%.3lf)!\n",(**i).getObj(0),(**i).getObj(1),(**i).getObj(2),s->getObj(0),s->getObj(1),(**i).getObj(2));
			}
			if (**i >> *s || **i == *s)
				return false;
			i++;
		}

		list< list<SolucaoEdgeSet *>::iterator >::iterator j = remover.begin();
		while (j != remover.end()) {
		    // remove do grid
			g.removeGrid( calcularGridPos(***j) );
			// remove a frequencia das arestas
			removerFrequencia( **j );
			delete( **j );
			// remove do conjunto pareto
			sol.erase( *j );
			j++;
		}

		SolucaoEdgeSet *t = new SolucaoEdgeSet(s->nEdges);
		*t = *s;
		// adiciona ao conjunto pareto
		sol.push_front( t );
		// adiciona a frequencia das arestas da solucao
		adicionarFrequencia( t );
		// adiciona ao grid
		g.addGrid( calcularGridPos(*t) );

		for(int k=0;k<NUMOBJETIVOS;k++) {
			rangeNovo[k].min = min(rangeNovo[k].min,t->getObj(k));
			rangeNovo[k].max = max(rangeNovo[k].max,t->getObj(k));
		}

		// se houve uma mudanca grande nos ranges (maior que 10% valor), atualizar o grid
		for (int k=0;k<NUMOBJETIVOS;k++) {
			if (fabs(rangeNovo[k].min-rangeAtual[k].min) > 0.1*rangeAtual[k].min || fabs(rangeNovo[k].max-rangeAtual[k].max) > 0.1*rangeAtual[k].max) {
				//fprintf(stderr,"Atualizando grid!\n");
				updateGrid();
				break;
			}
		}

       // ASS ( assert( confereGrid() ); )
		return true;
	}

	void printSet(FILE *f) {
		list<SolucaoEdgeSet *>::iterator i = sol.begin();
		SolucaoEdgeSet *s;

		fprintf(f,"Ranges = [%.2lf,%.2lf], [%.2lf,%.2lf], [%.2lf,%.2lf]\n",rangeAtual[0].min,rangeAtual[0].max,rangeAtual[1].min,rangeAtual[1].max, rangeAtual[2].min,rangeAtual[2].max);
		while (i != sol.end()) {
			s = *i;
			fprintf(f,"(%.10lf,%.10lf,%.10lf) -> %d\n",s->getObj(0),s->getObj(1),s->getObj(2),calcularGridPos(*s));
			i++;
		}
	}
	void printSetPoints(FILE *f) {
		list<SolucaoEdgeSet *>::iterator i = sol.begin();
		SolucaoEdgeSet *s;
		while (i != sol.end()) {
			s = *i;
			fprintf(f,"%.10lf %.10lf %.10lf\n",s->getObj(0),s->getObj(1), s->getObj(2));
			i++;
		}
	}

	void printAllSolutions(FILE *f) {
		list<SolucaoEdgeSet *>::iterator i = sol.begin();
		SolucaoEdgeSet *s;
		while (i != sol.end()) {
			s = *i;
			fprintf (f, "\n\n\n");
			s->printSolucao (f);
			i++;
		}
	}

	void printOnePoint(list<SolucaoEdgeSet *>::iterator i, FILE *f) {
		SolucaoEdgeSet *s;
		s = *i;
		fprintf(f,"%.10lf %.10lf %.10lf\n",s->getObj(0),s->getObj(1), s->getObj(2));
	}

	int getSize() {
		return sol.size();
	}

	SolucaoEdgeSet *getSolucao(int p) {
		int c = 0;
		list<SolucaoEdgeSet *>::iterator i = sol.begin();
		while (i != sol.end()) {
			if (c == p) return *i;
			i++;
			c++;
		}
		return NULL;
	}

	int getFrequencia(int a, int b) {
	    if (a > b) std::swap(a,b);
	    return frequencia[a][b];
	}


	void clear() {
		list<SolucaoEdgeSet *>::iterator i = sol.begin(), j;
		while (i != sol.end()) {
			j = i;
			i++;
			delete (*j);
		}
		sol.clear();
		g.clearGrid();
	}

	void copy(list<SolucaoEdgeSet *> solucoes){
		clear();
		list<SolucaoEdgeSet *>::iterator i = solucoes.begin();
		while (i!=solucoes.end()){
			SolucaoEdgeSet *t = new SolucaoEdgeSet(NUMEROVERTICES-1);
			*t = **i;
			sol.push_back( t );

			// adiciona a frequencia das arestas da solucao TODO: deletar estas linhas?
			adicionarFrequencia( t );
			// adiciona ao grid
			g.addGrid( calcularGridPos(*t) );

			for(int k=0;k<NUMOBJETIVOS;k++) {
				rangeNovo[k].min = min(rangeNovo[k].min,t->getObj(k));
				rangeNovo[k].max = max(rangeNovo[k].max,t->getObj(k));
			}

			// se houve uma mudanca grande nos ranges (maior que 10% valor), atualizar o grid
			for (int k=0;k<NUMOBJETIVOS;k++) {
				if (fabs(rangeNovo[k].min-rangeAtual[k].min) > 0.1*rangeAtual[k].min || fabs(rangeNovo[k].max-rangeAtual[k].max) > 0.1*rangeAtual[k].max) {
					//fprintf(stderr,"Atualizando grid!\n");
					updateGrid();
					break;
				}
			}
			i++;
		}
	} 
	bool confereGrid() {
	    unsigned s = 0;
	    for (int i=0;i<g.getSize();i++) s += g.getPositionCount(i);
	    return s == sol.size();
	}


	void confere(){
		list<SolucaoEdgeSet *>::iterator i = sol.begin();
		while (i!=sol.end()){
			if ((*i)->isTree()==false) {
				cout<<"ALGO ERRRADOOO"<<endl;
				exit(-1);
			}

			i++;
		}
	}
};

#endif
