#ifndef SOLUCAO_CPP
#define SOLUCAO_CPP

/*This code file was kindly provided by Monteiro */


#include "param_NSGAIII.h"

class Solucao {
	protected:
	double f[4];

	public:
	void setObj(int i, double v) {
		f[i] = v;
	}

	public:
	Solucao() {

	}

	// esta rotina funciona melhor quando o f é inteiro (o que é nosso caso)
	bool operator>> (Solucao &d) { // retorna true se f domina d.f
		return (f[0]<=d.f[0] && f[1]<=d.f[1] && f[2]<=d.f[2] && f[3]<=d.f[3] && (f[0]<d.f[0] || f[1]<d.f[1] || f[2]<d.f[2] || f[3]<d.f[3]));
	}

	void crossover(Solucao *pai, Solucao *mae) {
	}

	const double getObj(int o) {
		return f[o];
	}

	void operator=(Solucao &s) {
		f[0] = s.f[0];
		f[1] = s.f[1];
		f[2] = s.f[2];
		f[3] = s.f[3];
	}
	bool operator==(Solucao &s) {
		if (fabs(s.getObj(0)-getObj(0)) < EPS && fabs(s.getObj(1)-getObj(1)) < EPS && fabs(s.getObj(2)-getObj(2)) < EPS    && fabs(s.getObj(3)-getObj(3)) < EPS)
			return true;
		return false;
	}
};

#endif
