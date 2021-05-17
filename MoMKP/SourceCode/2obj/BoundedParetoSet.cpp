#ifndef BOUNDEDPARETOSET_CPP
#define BOUNDEDPARETOSET_CPP

#include "ParetoSet.cpp"
#include <cstdio>
#include <string>
#include <list>

#include "param_geral.h"

/*This code file was kindly provided by Monteiro */


using namespace std;

class BoundedParetoSet : public ParetoSet {

	public:

	BoundedParetoSet () { 

	}
	~BoundedParetoSet () { }

  
	bool adicionarSol(Knapsack *s) {
		//ASS ( assert( confereGrid() ); )
		list<Knapsack *>::iterator maisPopuloso = sol.begin();
		int maiorPositionCount = -1;
		if (sol.size() > 0) maiorPositionCount = getPositionCount( **sol.begin() );
		
		/* percorre o vetor de solucoes e de valores e, caso exista solucao dominada, retira e retorna true. caso contrario, retorna false */
		list<Knapsack *>::iterator i = sol.begin();
		list< list<Knapsack *>::iterator > remover;
		while (i != sol.end()) {
			// se a solucao que vai entrar domina a solucao i 
			if (*s >> **i) {
				remover.push_back(i); 
			}
			// se a solucao que vai entrar nao domina a solucao i, procuramos a solucao que está no
			// local mais populoso, caso o pareto esteja no tamanho maximo
			if (remover.size() == 0 && getSize()+1 > MAXARCSIZE) {
				int positionCountAtual = getPositionCount(**i); // posiçao de i na grid
				if (maiorPositionCount < positionCountAtual) { // pega o maximo entre os dois. Desta posiçao será removida a soluçao para "sol" entrar
					maiorPositionCount = positionCountAtual;
					maisPopuloso = i;
				}
			}
			
			if (**i >> *s || **i == *s){ // caso algum "i" domine "s", para e retorna false
				return false;
			}
			i++;
		}

		// se a solucao que vai entrar nao domina nenhuma e o tamanho do conjunto pareto ja esta no maximo
		// (se nenhuma solucao vai sair do conjunto), remove a mais populosa
		if (remover.size() == 0 && getSize()+1 > MAXARCSIZE) {
		
			remover.push_back(maisPopuloso);
		}

		//fprintf(stderr,"getSize = %d %d\n",getSize(),sol.size());
			
		list< list<Knapsack *>::iterator >::iterator j = remover.begin();
		while (j != remover.end()) {
			// remove do grid
			g.removeGrid( calcularGridPos(***j) );

			delete( **j ); 
			// remove do conjunto pareto
			sol.erase( *j );
			// insere na lixeira // by felipe
			j++;
		}

		Knapsack *t = new Knapsack(my_momkp->size_, my_momkp->dimension_);
		*t = *s;
		// adiciona ao conjunto pareto
		sol.push_front( t );
		if (sol.size() > MAXARCSIZE) fprintf(stderr,"ERRO!\n");
		
		// adiciona ao grid
		g.addGrid( calcularGridPos(*t) );

		for(int k=0;k<my_momkp->dimension_;k++) {
			rangeNovo[k].min = min(rangeNovo[k].min,t->profit(k));
			rangeNovo[k].max = max(rangeNovo[k].max,t->profit(k));
		}

		// se houve uma mudanca grande nos ranges (maior que 10% valor), atualizar o grid
		for (int k=0;k<my_momkp->dimension_;k++) {
			if (fabs(rangeNovo[k].min-rangeAtual[k].min) > 0.1*rangeAtual[k].min || fabs(rangeNovo[k].max-rangeAtual[k].max) > 0.1*rangeAtual[k].max) {
				//fprintf(stderr,"Atualizando grid!\n");
				updateGrid();
				break;
			}
		}
		return true;
	}

	
	
};

#endif
