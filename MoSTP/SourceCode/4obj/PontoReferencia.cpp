#ifndef PONTOREFERENCIA_H_
#define PONTOREFERENCIA_H_

#include "param_NSGAIII.h"
#include "SolucaoEdgeSet.cpp"

class PontoReferencia {
	public:
		size_t quantMembros;			// numero de membros associados ao ponto
		vector<pair<SolucaoEdgeSet*, double> > membros; // membros (SolucaoEdgeSet e distancia)
	public:
		vector<double> ponto_Referencia; 
	public:
		PontoReferencia() {
			quantMembros = 0;
		}
		PontoReferencia(size_t n_obj) {
			quantMembros = 0;
			for (size_t i = 0; i < n_obj; ++i)
				ponto_Referencia.push_back(0.0);
		}
		virtual ~PontoReferencia() {}


		size_t getQuantMembros() {
			return quantMembros;
		}
		bool tem_membros() {
			return !membros.empty();
		}

		SolucaoEdgeSet * get_membro_mais_proximo() {

			int id = -1;
			double min_dist = INT_MAX;
			for (size_t i = 0; i < membros.size(); ++i) {
				if (membros[i].second < min_dist) {
					id = i;
					min_dist = membros[i].second;
				}
			}
			if (id!=-1) return membros[id].first;
			else return NULL;
		}

		SolucaoEdgeSet * get_membro_random() {
			return membros[IRandom(0, membros.size()-1)].first;
		}

		void sum_member() {
			quantMembros += 1;
		}

		void add_member(SolucaoEdgeSet * sol, double distance) {
			membros.push_back(make_pair(sol, distance));
		}
		void remove_member(SolucaoEdgeSet * sol) {
			for (size_t i = 0; i < membros.size(); ++i) {
				if (*membros[i].first == *sol) {
					membros.erase(membros.begin() + i);
				//	cout<<"REMOVED"<<endl;
					return;
				}
			}
		}

		void clear_ponto_Referencia() {
			quantMembros = 0;
			membros.clear();
		}

		//void show_rpoint() {
		//	for (size_t i = 0; i < numObjectives; i++)
		//		printf("%f ", ponto_Referencia[i]);
		//	printf("\n");
		//}

};








#endif 