
#ifndef UNIONFIND_CPP
#define UNIONFIND_CPP
/*This code file was kindly provided by Monteiro */


#include "param_NSGAIII.h"
#include <cstring>



class UnionFind {
	private:
	short pai[NUMEROVERTICES], altura[NUMEROVERTICES];
	short getClass(short a) {
		if (a != pai[a])
			pai[a] = getClass(pai[a]);
		return pai[a];
	}

	public:
	UnionFind() {
		clear();
	}
	~UnionFind() {
	}
	bool sameClass(int a, int b) {
		return getClass(a) == getClass(b);
	}
	void unionClass(int a, int b) {
		if (a != b) {
			int cl[2] = {getClass(a),getClass(b)};
			if (cl[0] != cl[1]) {
				if (altura[cl[0]] > altura[cl[1]])
					pai[cl[1]] = cl[0];
				else {
					pai[cl[0]] = cl[1];
					if (altura[cl[0]] == altura[cl[1]]) altura[cl[1]]++;
				}
			}
		}
	}
	void clear() {
		memset(altura,0,NUMEROVERTICES*2);
		for (short i=0;i<NUMEROVERTICES;i++) {
			pai[i] = i;
		}
	}
};


#endif
