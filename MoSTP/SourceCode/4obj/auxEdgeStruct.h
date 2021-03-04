
#ifndef AUXEDGESTRUCT_H
#define AUXEDGESTRUCT_H

/*This code file was kindly provided by Monteiro */

struct auxEdgeSet {

	int a; // origem
	int b; // destino
	double fit; // fitness de acordo com um vetor de escalarizaçao

};

inline bool auxEdgeSetComp (auxEdgeSet e1, auxEdgeSet e2) {

	return (e1.fit < e2.fit);
}


struct auxEdgeStruct {

	int a;
	int b;
	double c1;
	double c2;
	double c3;
	double c4;
	double fit;

};


// ordena lexicografiacamente pelo primeiro objetivo
inline bool lexicografica_obj1 (auxEdgeStruct e1, auxEdgeStruct e2) {

	return (e1.c1 < e2.c1); 
}


// ordena lexicografiacamente pelo segundo objetivo
inline bool lexicografica_obj2 (auxEdgeStruct e1, auxEdgeStruct e2) {

	return (e1.c2 < e2.c2); 
}

// ordena lexicografiacamente pelo segundo objetivo
inline bool lexicografica_obj3 (auxEdgeStruct e1, auxEdgeStruct e2) {

	return (e1.c3 < e2.c3); 
}


// ordena lexicografiacamente pelo quarto objetivo
inline bool lexicografica_obj4 (auxEdgeStruct e1, auxEdgeStruct e2) {

	return (e1.c4 < e2.c4); 
}




// ordena lexicografiacamente pelo primeiro objetivo
inline bool lexicografica_max_obj1 (auxEdgeStruct e1, auxEdgeStruct e2) {

	return (e1.c1 > e2.c1); 
}


// ordena lexicografiacamente pelo segundo objetivo
inline bool lexicografica_max_obj2 (auxEdgeStruct e1, auxEdgeStruct e2) {

	return (e1.c2 > e2.c2); 
}


// ordena lexicografiacamente pelo segundo objetivo
inline bool lexicografica_max_obj3 (auxEdgeStruct e1, auxEdgeStruct e2) {

	return (e1.c3 > e2.c3); 
}

// ordena lexicografiacamente pelo segundo objetivo
inline bool lexicografica_max_obj4 (auxEdgeStruct e1, auxEdgeStruct e2) {

	return (e1.c4 > e2.c4); 
}

struct fitVecNode {
	int a, b;
	double c;
};


struct ordenacao {
	int id; //identificador
	double v_esc [NUMOBJETIVOS]; //vetor de escalarização associado
	auxEdgeSet arestas [(NUMEROVERTICES*(NUMEROVERTICES-1))/2]; //conjunto ordenado de arestas de acordo com o v_esc
};

#endif
