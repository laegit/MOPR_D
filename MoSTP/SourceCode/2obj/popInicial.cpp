#ifndef POP_INICIAL
#define POP_INICIAL

#include "param_NSGAII.h"
#include "SolucaoEdgeSet.cpp"

#include <iostream>
using namespace std;



/*
	População inicial é gerada conforme sugerido por Monteiro (2010)
	cria no maximo 90% da populaçao com o rmcprim (NAO CONSIDERA nem individuos repetidos nem dominados)
	cria no minimo 10% da populaçao com o randomWolk (CONSIDERA individuos repetidos e dominados)
*/

// precondiçao : tamanhoPop-init > 0
void gerarPopulacao1(SolucaoEdgeSet **populacao, int tamanhoPop, int init, int &contAval){
	double lambda[NUMOBJETIVOS] = {0.0, 1.0};
	int quant_prim = (int)((tamanhoPop-init)*0.9); // 90% com rmc_prim
	double factor = (quant_prim==0 ? 0 : 1.0/quant_prim); // fator de incremetar o lambda
	int contPop = init;//0; // contador para indexar o vetor populacao
	int quant = 0;
	while (quant<quant_prim && lambda[0]<1.0){
		rmcPrim( *populacao[contPop], lambda, contAval);
		lambda[0]+=factor;
		lambda[1]-=factor;
		bool ha = false;
		for (int j=0;j<contPop && ha==false;j++) {
			if (*populacao[contPop] == *populacao[j]) { // note que aqui vericamos o conteudo, e nao o ponteiro
				ha = true;
			}
		}

		if (ha==false){ 
			bool resp = true; 
			for (int pppp = 0; pppp<contPop && resp==true; pppp++){ // [0 ... contPop-1]. A posicao contPop contem a nova soluçao gerada com o rmcprim
				if (*populacao[pppp] >> *populacao[contPop]) resp=false;
			}
			if (resp==true){ // aceita se pop[contPop] nao for dominada por ninguém no arquivo
				contPop++;
				quant++;
			} else if (genrand64_real3() < 0.6) { // aceita com 60% de chance
				contPop++;
				quant++;
			} 
		}
	}
	
	while (contPop<tamanhoPop){
		populacao[contPop]->doRandomWalk(contAval);
		contPop++;
	}
}

#endif