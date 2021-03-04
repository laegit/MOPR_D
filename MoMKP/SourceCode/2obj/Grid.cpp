#ifndef GRID_CPP
#define GRID_CPP

#include <cstring>
#include "momkp.cpp"

/*This code was provided by Monteiro*/

#include "param_geral.h"
extern MOMKP *my_momkp;

class Grid {
	private:
	int tam; // my_momkp->dimension_ = quantidade de objetivos
	int *grid;

	public:
		Grid(){
			tam = 1<<(my_momkp->dimension_*PROFUNDIDADEGRID);
			grid = new int[1<<(my_momkp->dimension_*PROFUNDIDADEGRID)];
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
