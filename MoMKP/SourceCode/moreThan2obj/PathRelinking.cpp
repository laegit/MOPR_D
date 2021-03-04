
#ifndef PATHRELINKING_CPP
#define PATHRELINKING_CPP

#include "Knapsack.cpp"

class PathRelinking{
	public:
		virtual void PR(ParetoSet *retorno, Knapsack *s1, Knapsack *s2, int &contAval) = 0; 
			
		/*
			Se for passado um EpislonApproximatePaetoSet, temos a abordagem e-dominância
			Se for passado um BoundedParetoSet, temos a abordagem de dominância classica de pareto
			Se for passado uma ListaSimples, temos a abordagem Intermediate Pool
		*/

		// virtual void PR_intermediate_pool(list<Knapsack *> &retorno, Knapsack *s1, Knapsack *s2, int &contAval) = 0; 
		// virtual void PR_epsilon_dominancia(EpislonApproximatePaetoSet *retorno, Knapsack *s1, Knapsack *s2, int &contAval) = 0; 

		// // guarda itens de s1 que nao estao em s2. 
		// std::vector<size_t> inS1_NotInS2(Knapsack *s1, Knapsack *s2){ // O(n)
		// 	std::vector<size_t> amostral(s1->getSize(), false);
		// 	for (int i=0; i<s1->getSize(); i++){ // varre os itens de s1 que nao estao em s2
		// 		if (s1->at[i]==true && s2->at[i]==false) {
		// 			amostral[i] = true; 
		// 		}
		// 	}
		// 	return amostral;
		// }
};

#endif
