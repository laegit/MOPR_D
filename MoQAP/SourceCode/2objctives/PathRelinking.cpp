
#ifndef PATHRELINKING_CPP
#define PATHRELINKING_CPP

#include "Solution.cpp"
#include "BoundedParetoSet.cpp"
#include "param_geral.h"


class PathRelinking{
	public:
		virtual void PR(BoundedParetoSet *retorno, Solution *s1, Solution *s2, int &contAval) = 0; 
		// PathRelinking(){}
		


		// baseado no trabalho de H. Li and D. Landa-Silva (2009) An Elitist GRASP Metaheuristic for the Multi-objective Quadratic Assignment Problem
		//C^k ( π' ) = C^k ( π ) + Δ ( π , k , i , j ) , k = 1 , . . . , m .
		// posiçoes i e j sao trocadas
		// Δ é computado pela equacao 4 do paper supracitao
		// VALE TANTO PARA INSTANCIAS ASSIMETRICAS QUANTO SIMETRICAS  
		// Calcula o valor da funcao para o k-ésimo objetivo. Retorna o valor computado.
		// sol é a SOLUCAO **ANTES** da troca
		double delta(std::vector<int> assignment, double obj_k, int i, int j){ // *******  O(n) ******* //
			double reto = (distances[assignment[j]][assignment[j]] - distances[assignment[i]][assignment[i]])*(flows[obj_k][i][i] - flows[obj_k][j][j])
					    + (distances[assignment[j]][assignment[i]] - distances[assignment[i]][assignment[j]])*(flows[obj_k][i][j] - flows[obj_k][j][i]);
			for (int s = 0; s < n_locations; ++s){
				if (s!=i && s!=j){
					reto += (distances[assignment[s]][assignment[j]] - distances[assignment[s]][assignment[i]])*(flows[obj_k][s][i] - flows[obj_k][s][j])
						 +	(distances[assignment[j]][assignment[s]] - distances[assignment[i]][assignment[s]])*(flows[obj_k][i][s] - flows[obj_k][j][s]);
				}
			}
			return reto;
		}
		
};

#endif
