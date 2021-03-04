#ifndef GRID_CPP
#define GRID_CPP

#include <cstring>

/*This code was provided by Monteiro*/

#include "param_geral.h"

class Grid {
	private:
	int tam; //n_objetivos = quantidade de objetivos
	int *grid;

	public:
		Grid(){
			tam = 1<<(n_objetivos*PROFUNDIDADEGRID);
			grid = new int[1<<(n_objetivos*PROFUNDIDADEGRID)];
		}
	int getPositionCount(int p) {
		if (p < 0 || p >= tam) return -1;
		return grid[p];
	}
	void addGrid(int p) {
		grid[p]++;
	}
	void removeGrid(int p) {
		grid[p]--;
	}
	void clearGrid() {
		memset(grid,0,sizeof(grid[0])*tam);
	}
	int getSize() {
	    return tam;
	}
};

#endif
