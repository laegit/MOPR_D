#ifndef ZG_H
#define ZG_H

//// PARÂMETROS FORNECIDOS PELO IRACE

#include "rand64/mersenne64.c"


// FIXO
#define PI 3.14159265
#define f(k,i,j) custos[k][i][j] // objetivo k, vertice i j 
#define EPS 1e-9 // qualquer coisa menor que esse valor, é considerado 0
#define INF 1e9
#define PROFUNDIDADEGRID 5 //5 ///AQUI///
#define NUMOBJETIVOS 2 // FIXO
#define NUMEROVERTICES 1000 // Varia
#define NUMEROARESTAS (NUMEROVERTICES-1)
#define MAXARCSIZE 300
#define NUM_AVALIACOES 1000000


#define TAMANHOPOPULACAO_NSGAII 300  // tamanho da populaçao
// #define QUANTGERACOES 30  //quantidade de iteracoes do nsgaii (IRACE)
#define TAXADECRUZAMENTO 0.97 // (IRACE para o NSGA-II)
#define TAXADEMUTACAO 0.04 // (IRACE para o NSGA-II)



#define NUMSUMPROBLEMAS  127 
#define NUM_VIZINHOS 20 
#define MAX_SUB_POPULATION 42 

#endif




