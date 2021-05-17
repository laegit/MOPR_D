#ifndef SOLUCAO_CPP
#define SOLUCAO_CPP

/*This code file was kindly provided by Monteiro */


// extern int contAvaliacoesObjetivoo;

#include "param_NSGAII.h"

class Solucao {
	protected:
	double f[2];

	public:
	void setObj(int i, double v) {
		f[i] = v;
	}

	public:
	Solucao() {

	}
	bool operator>> (Solucao &d) { // 2 obj 
		double ele[2] = {d.getObj(0),d.getObj(1)}, diff[2] = { ele[0] - f[0] , ele[1] - f[1] };
		if (diff[0] > EPS && (diff[1] > EPS || fabs(diff[1]) < EPS))
			return true;
		if ((diff[0] > EPS || fabs(diff[0]) < EPS) && diff[1] > EPS)
			return true;
		return false;
	}

	void crossover(Solucao *pai, Solucao *mae) {
	}

	const double getObj(int o) {
		return f[o];
	}

	void operator=(Solucao &s) {
		f[0] = s.f[0];
		f[1] = s.f[1];
	}
	bool operator==(Solucao &s) {
		if (fabs(s.getObj(0)-getObj(0)) < EPS && fabs(s.getObj(1)-getObj(1)) < EPS)
			return true;
		return false;
	}
};

#endif
