#ifndef PARETOSET_CPP
#define PARETOSET_CPP

/*This code file was kindly provided by Monteiro */


#include <list>
#include <map>
#include "Solution.cpp"
#include "Grid.cpp"



using namespace std;

typedef struct {
	double min, max;
} range;

class ParetoSet {
	public:
	double hyperolume;
	list<Solution *> sol; // TODO: deveria ser protected
	range rangeNovo[10], rangeAtual[10]; // 10: quantidade maxima de objetivos. Nossas instâncias vao até 8
	Grid g;

	int calcularGridPos(Solution &s) {
		int bit = 0;
		int gridPos = 0;
		for (int obj=0;obj<n_objetivos;obj++) { // dimension_ = quantiObjetivos
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

		list<Solution *>::iterator it = sol.begin();
		reiniciarRanges();
		while (it != sol.end()) {
			for (int k=0;k<n_objetivos;k++) {
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
		for (int k=0;k<n_objetivos;k++) {
			rangeAtual[k].min = rangeNovo[k].min = INF;
			rangeAtual[k].max = rangeNovo[k].max = -INF;
		}
		#undef INF
	}


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
	}
	virtual ~ParetoSet() {
		clear();
	}

	int getPositionCount(Solution &s) {
		return g.getPositionCount( calcularGridPos(s) );
	}

	int getPositionCount(int p) {
		return g.getPositionCount( p );
	}

	list<Solution *> getElementos() {
		return sol;
	}

	/* Complexidade: O(n) */
	virtual bool adicionarSol(Solution *s) {
		/* percorre o vetor de solucoes e de valores e, caso exista solucao dominada, retira e retorna true. caso contrario, retorna false */
		list<Solution *>::iterator i = sol.begin();
		list< list<Solution *>::iterator > remover;
		while (i != sol.end()) {
			if (*s >> **i) { // se *s domina *i
				remover.push_back(i);
			}
			if (**i >> *s || **i == *s)
				return false;
			i++;
		}

		list< list<Solution *>::iterator >::iterator j = remover.begin();
		while (j != remover.end()) {
		    // remove do grid
			g.removeGrid( calcularGridPos(***j) );
			// remove a frequencia das arestas
			delete( **j );
			// remove do conjunto pareto
			sol.erase( *j );
			j++;
		}

		Solution *t = new Solution(n_locations, n_objetivos);
		*t = *s;
		// adiciona ao conjunto pareto
		sol.push_front( t );
		// adiciona ao grid
		g.addGrid( calcularGridPos(*t) );

		for(int k=0;k<n_objetivos;k++) {
			rangeNovo[k].min = min(rangeNovo[k].min,t->getObj(k));
			rangeNovo[k].max = max(rangeNovo[k].max,t->getObj(k));
		}

		// se houve uma mudanca grande nos ranges (maior que 10% valor), atualizar o grid
		for (int k=0;k<n_objetivos;k++) {
			if (fabs(rangeNovo[k].min-rangeAtual[k].min) > 0.1*rangeAtual[k].min || fabs(rangeNovo[k].max-rangeAtual[k].max) > 0.1*rangeAtual[k].max) {
				//fprintf(stderr,"Atualizando grid!\n");
				updateGrid();
				break;
			}
		}

       // ASS ( assert( confereGrid() ); )
		return true;
	}


	void printSetPoints(FILE *f) {
		list<Solution *>::iterator i = sol.begin();
		Solution *s;
		while (i != sol.end()) {
			s = *i;
			for (int obj=0; obj<n_objetivos; obj++){
		      fprintf(f,"%f ",s->getObj(obj));
		    }
			fprintf(f,"\n");
			i++;
		}
	}


	int getSize() {
		return sol.size();
	}

	Solution *getSolucao(int p) {
		int c = 0;
		list<Solution *>::iterator i = sol.begin();
		while (i != sol.end()) {
			if (c == p) return *i;
			i++;
			c++;
		}
		return NULL;
	}


	void clear() {
		list<Solution *>::iterator i = sol.begin(), j;
		while (i != sol.end()) {
			j = i;
			i++;
			delete (*j);
		}
		sol.clear();
		g.clearGrid();
	}

	// void copy(list<Solution *> solucoes){
	// 	clear();
	// 	list<Solution *>::iterator i = solucoes.begin();
	// 	while (i!=solucoes.end()){ 
	// 		Solution *t = new Solution(my_momkp->size_, n_objetivos);
	// 		*t = **i;
	// 		sol.push_back( t );

		
	// 		// adiciona ao grid
	// 		g.addGrid( calcularGridPos(*t) );

	// 		for(int k=0;k<n_objetivos;k++) {
	// 			rangeNovo[k].min = min(rangeNovo[k].min,t->getObj(k));
	// 			rangeNovo[k].max = max(rangeNovo[k].max,t->getObj(k));
	// 		}

	// 		// se houve uma mudanca grande nos ranges (maior que 10% valor), atualizar o grid
	// 		for (int k=0;k<n_objetivos;k++) {
	// 			if (fabs(rangeNovo[k].min-rangeAtual[k].min) > 0.1*rangeAtual[k].min || fabs(rangeNovo[k].max-rangeAtual[k].max) > 0.1*rangeAtual[k].max) {
	// 				//fprintf(stderr,"Atualizando grid!\n");
	// 				updateGrid();
	// 				break;
	// 			}
	// 		}
	// 		i++;
	// 	}
	// } 
	bool confereGrid() {
	    unsigned s = 0;
	    for (int i=0;i<g.getSize();i++) s += g.getPositionCount(i);
	    return s == sol.size();
	}
};

#endif
