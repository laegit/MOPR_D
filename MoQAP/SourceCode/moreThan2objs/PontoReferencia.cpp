#ifndef PONTOREFERENCIA_H_
#define PONTOREFERENCIA_H_

#include "random.cpp"
#include "Solution.cpp"

class PontoReferencia {
	public:
		size_t quantMembros;			// numero de membros associados ao ponto
		vector<pair<Solution*, double> > membros; // membros (Solution e distancia)
	public:
		vector<double> ponto_Referencia; 
	public:
		PontoReferencia() {
			quantMembros = 0;
		}
		PontoReferencia(size_t num_objetivos) {
			quantMembros = 0;
			for (size_t i = 0; i < num_objetivos; ++i)
				ponto_Referencia.push_back(0.0);
		}
		virtual ~PontoReferencia() {}


		size_t getQuantMembros() {
			return quantMembros;
		}
		bool tem_membros() {
			return !membros.empty();
		}

		Solution * get_membro_mais_proximo() {

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

		Solution * get_membro_random() {
			return membros[pick_a_number(0, membros.size()-1)].first;
		}

		void sum_member() {
			quantMembros += 1;
		}

		void add_member(Solution * sol, double distance) {
			membros.push_back(make_pair(sol, distance));
		}
		void remove_member(Solution * sol) {
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