/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2019)

Multi-objective Quadratic Assignment Problem

=====================================================================================
*/

#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <sys/times.h>
#include <sys/types.h>
#include <iostream>
#include <math.h>       /* acos, cos */
#include <climits>
#include <algorithm>    // std::sort
#include <pthread.h>

#include "Solution.cpp"
#include "random.cpp"
#include "lerInstance.cpp"
#include "nsga3.cpp"
#include "BoundedParetoSet.cpp"
#include "PathRelinking.cpp"
#include "Version1.cpp"
#include "Version2.cpp"
#include "Version3.cpp"
#include "Version4.cpp"
#include "Version6.cpp"
#include "Version7.cpp"
#include "Version5_D.cpp"
#include "hyp_ind.cpp"


using namespace std;


std::vector< std::vector<double> > distances;
// matriz de distâncias. 
// distances[i][j] == distânicia entre a localidade i e j

std::vector< vector< std::vector<double> > > flows;
//flows[obj][i][j] = o obj-ésimo fluxo entre i e j


long init_seed; // semente da instância
int n_locations;// = 30;  // number of facilities/locations in the QAP
int n_objetivos; //= 2;    // number of objectives
double corr; // = 0;   // correlation between 1st and all other objectives
double amp; // relates to the correlation
double offset; // relates to the correlation
int max_flow; //= 100;
int max_dist; // = 100;
// int i,j,k;
// double r1,r2;
// int f1,fk;
// int d1;


void print_output()
{
  int i,j,k;
  printf("facilities: %d objectives: %d max_distances: %d max flows: %d seed %ld\n", n_locations, n_objetivos, max_dist, max_flow, init_seed);
  for(i=0;i<n_locations;i++)
    {
      for(j=0;j<n_locations;j++)
	{
	  printf("%2f ",distances[i][j]);
	}
      printf("\n");
    }
 

  for(k=0;k<n_objetivos;k++)
    {
      printf("\n");
      for(i=0;i<n_locations;i++)
	{
	  for(j=0;j<n_locations;j++)
	    {
	      printf("%2f ",flows[k][i][j]);
	    }
	  printf("\n");
	}
    }
  
}

void getQAP(char *infile)
{
  int i,j,k;
  double dummy,dval;

  FILE *fp;
  if ((fp = fopen(infile, "rb")))
    {
      getint(fp, &n_locations, GDFALSE);
      getint(fp, &n_objetivos, GDFALSE);
      getint(fp, &max_dist, GDFALSE);
      getint(fp, &max_flow, GDFALSE);
      getdouble(fp, &corr, GDFALSE);  // this may be the correlation value (uniform instance) or overlap parameter (real-like instance)
      getlong(fp, &init_seed, GDTRUE);
      do
	{
	  dval=getdouble(fp, &dummy, GDTRUE);
	}while(dval>=1);

 
 
      for(i=0;i<n_locations;i++)
	{
		std::vector<double> linha;
		
	  for(j=0;j<n_locations;j++)
	    {
			int aaux;
	      	getint(fp, &aaux, GDFALSE);
	      	linha.push_back(aaux);
	    }
	    distances.push_back(linha);
	}

    for (k=0;k<n_objetivos;k++){
		std::vector<std::vector<double> > obj_k;
	  	for(i=0;i<n_locations;i++){
	    	std::vector<double> linha;
	      	for(j=0;j<n_locations;j++){
		 		int auxx;
		  		getint(fp, &auxx, GDFALSE);
		  		linha.push_back(auxx);
			}
			obj_k.push_back(linha);
	    }
	    flows.push_back(obj_k);
	}
	      
      fclose(fp);
    }
  else
    {
      fprintf(stderr, "Could not open the fille.\n");
    }
}
  

void getLambdaRandom(std::vector<double> &lambdaAux){
    double sumL = 0.0;
    for (int i = 0; i < n_objetivos-1; ++i){
      double l = pick_a_number(0.0, 1.0-sumL); //Random(0.0,1.0-l1);
      sumL+=l;
      lambdaAux[i] = l;
    }
    lambdaAux[n_objetivos-1] = 1.0 - sumL;
  }






// Ver Arroyo Vieira e Vianna (2008)
void calculaVetoresLambda_SISTEMATICO(vector<double> lambdas[NUMSUMPROBLEMAS], int numObjts, int s){
	int cont =0;
	if (numObjts==2){ // s-1 
		// para atingir 10^6 vetores, necessario s=10^6-1
		double l1, l2;
		for (int i=0; i<=s && cont<NUMSUMPROBLEMAS; i++){
			l1 = i;
			l2 = s - l1;
			std::vector<double> lam = {static_cast<double>(l1/s), static_cast<double>(l2/s)};
			lambdas[cont++] = lam;
		}

	} else if (numObjts==3){ // (s+1)(s+2)/2
		// s = 1417 e s=1418
		double l1, l2, l3;
		for (int i=0; i<=s && cont<NUMSUMPROBLEMAS; i++){
			l1 = i;
			for (int j=0; j<=s-i && cont<NUMSUMPROBLEMAS; j++){
				l2 = j;
				l3 = s-i - l2;
				std::vector<double> lam = {static_cast<double>(l1/s), static_cast<double>(l2/s),  static_cast<double>(l3/s)};
				lambdas[cont++] = lam;
			}
		}

	} else if (numObjts==4){ // quantidade vetores: 1/6 (s^3 + 6 s^2 + 11s + 6)
		// s = 179 ou s=180
		double l1, l2, l3, l4;
		for (int i=0; i<=s  && cont<NUMSUMPROBLEMAS; i++){
			l1 = i;
			for (int j=0; j<=s-i && cont<NUMSUMPROBLEMAS; j++){
				l2 = j;
				for (int z=0; z<=s-i-j && cont<NUMSUMPROBLEMAS; z++){
					l3 = z;
					l4 = s-i-j - l3;
					std::vector<double> lam = {static_cast<double>(l1/s), static_cast<double>(l2/s),  static_cast<double>(l3/s), static_cast<double>(l4/s)};
					lambdas[cont++] = lam;
				}
			}
		}
	} else if (numObjts==5){ // quantidade de vetores: 1/24 (s^4 + 10 s^3 + 35 s^2 + 50 s + 24)
		//para atingir 10^6 de vetores, coloca s=67 ou s=68
		double l1, l2, l3, l4, l5;
		for (int i=0; i<=s  && cont<NUMSUMPROBLEMAS; i++){
			l1 = i;
			for (int j=0; j<=s-i && cont<NUMSUMPROBLEMAS; j++){
				l2 = j;
				for (int z=0; z<=s-i-j && cont<NUMSUMPROBLEMAS; z++){
					l3 = z;
					for (int w=0; w<=s-i-j-z && cont<NUMSUMPROBLEMAS; w++){
						l4 = w;
						l5 = s-i-j-z-l4;
						std::vector<double> lam = {static_cast<double>(l1/s), static_cast<double>(l2/s),  static_cast<double>(l3/s), static_cast<double>(l4/s), static_cast<double>(l5/s)};
						lambdas[cont++] = lam;
					}
				}
			}
		}
	} else if (numObjts==8){ //quantidade de vetores: 1/5040 (5040 + 13068 s + 13132 s^2 + 6769 s^3 + 1960 s^4 + 322 s^5 + 28 s^6 + s^7)
		//para atingir 10^6 de vetores, coloca s=20 ou s=21
		double l1, l2, l3, l4, l5, l6, l7, l8;
		for (int i=0; i<=s && cont<NUMSUMPROBLEMAS; i++){
			l1 = i;
			for (int j=0; j<=s-i && cont<NUMSUMPROBLEMAS; j++){
				l2 = j;
				for (int z=0; z<=s-i-j && cont<NUMSUMPROBLEMAS; z++){
					l3 = z;
					for (int w=0; w<=s-i-j-z && cont<NUMSUMPROBLEMAS; w++){
						l4 = w;
						for (int k=0; k<=s-i-j-z-w && cont<NUMSUMPROBLEMAS; k++){
							l5 = k;
							for (int o=0; o<=s-i-j-z-w-k && cont<NUMSUMPROBLEMAS; o++){
								l6 = o;
								for (int p=0; p<=s-i-j-z-w-k-o && cont<NUMSUMPROBLEMAS; p++){
									l7 = p;
									l8 = s-i-j-z-w-k-o - l7;;
									std::vector<double> lam = {static_cast<double>(l1/s), static_cast<double>(l2/s),  static_cast<double>(l3/s), static_cast<double>(l4/s), static_cast<double>(l5/s), static_cast<double>(l6/s), static_cast<double>(l7/s), static_cast<double>(l8/s)};
									lambdas[cont++] = lam;
								}
							}
						}					
					}
				}
			}
		}
	}
}










/* --------------------------------------------- */
/* -------- FUNCAO *NOVA* PARA O PR ---------- */
/* --------------------------------------------- */

BoundedParetoSet * rodaPR(BoundedParetoSet * conjAprox, PathRelinking *pr, int &contPares){
     int aval = 0;
     contPares = 0;
     BoundedParetoSet *retorno = new BoundedParetoSet();

     // copia previamente
     for (list<Solution *>::iterator it_1=conjAprox->sol.begin(); it_1!=conjAprox->sol.end(); it_1++)
        retorno->adicionarSol(*it_1); // copia, previamente, o conjunto inicial.
            
     // O PR(o,d) adiciona 'o' no conjunto. 
     // Como nós chamamos PR(o,d) e PR(d,o), entao 
     // todos os itens de conjAprox sao adicionados no retorno.

      while(aval<NUM_AVALIACOES){
        list<Solution *>::iterator it_1=conjAprox->sol.begin();
         while (it_1!=conjAprox->sol.end() && aval<NUM_AVALIACOES) {
            list<Solution *>::iterator it_2=conjAprox->sol.begin();
            while (it_2!=conjAprox->sol.end() && aval<NUM_AVALIACOES){
                
                if (it_1!=it_2){
                    pr->PR(retorno, *it_1, *it_2, aval);
                    contPares++;
                }

                it_2++;
            }
            it_1++;
         }
      }
     return retorno;
} 


// 1: nome do arquivo, 2: semente
int main(int argc, char *argv[])
{
	struct tms tempoAntes, tempoDepois;
	int seed = atoi(argv[2]); // seed
	randomize(seed);
	std::srand(seed);

	cout<<"Semente = "<<seed<<endl;
	getQAP(argv[1]); // LER INSTANCIA
	
	
	Hypervolume *hyp = new Hypervolume();




	BoundedParetoSet *resultado = new BoundedParetoSet();
	 
	 times(&tempoAntes);
	
	nsga3 *nii = new nsga3();
	Solution **population = nii->executar();
	//cout<<endl;
	// RESULTADO: apenas as nao-dominadas (primeiro front) da populacao
	for (int iii=0; iii<TAMANHOPOPULACAO_NSGAII; iii++){
		//population[iii]->print_solution(stdout);
		resultado->adicionarSol(population[iii]);
		
	}
	times(&tempoDepois);
//	

	cout << "Tempo_NSGAIII = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
    //fprintf(stdout, "HV_NSGAII = %.9e\n", hyp->getHypervolume(resultado, lbound, ubound));
    cout<<"Num_solucoes_NSGAIII = "<<resultado->getSize()<<endl;
	
	if (resultado->getSize()>1){
			PathRelinking *pr1 = new Version1();
			PathRelinking *pr2 = new Version2(pick_a_number(0, n_objetivos-1));
			PathRelinking *pr3 = new Version3();
			std::vector<double> lambdaAux(n_objetivos, 0.0);
			getLambdaRandom(lambdaAux);
			PathRelinking *pr4 = new Version4(lambdaAux); 
			PathRelinking *pr6 = new Version6();
			PathRelinking *pr7 = new Version7();
			//PathRelinking *pr8 = new Version8();


			std::vector<double> subproblemas[NUMSUMPROBLEMAS];
			if (n_objetivos==2){
		    	calculaVetoresLambda_SISTEMATICO(subproblemas, n_objetivos,NUMSUMPROBLEMAS-1); 
		    } else if (n_objetivos==3){
		    	calculaVetoresLambda_SISTEMATICO(subproblemas, n_objetivos,12); 
		    }else if (n_objetivos==4){
		    	calculaVetoresLambda_SISTEMATICO(subproblemas, n_objetivos,7); 
		    }

			int contPares = 0;

				times(&tempoAntes);
			BoundedParetoSet *newParetoSetAprox_version1 = rodaPR(resultado, pr1, contPares);

			 times(&tempoDepois);
			cout << "Tempo_PR1_NOVO = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
			//fprintf(stdout, "HV_PR1 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version1, lbound, ubound));


			times(&tempoAntes);
			BoundedParetoSet *newParetoSetAprox_version2 = rodaPR(resultado, pr2, contPares);

			 times(&tempoDepois);
			cout << "Tempo_PR2 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
			//fprintf(stdout, "HV_PR2 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version2, lbound, ubound));


			times(&tempoAntes);
			BoundedParetoSet *newParetoSetAprox_version3 = rodaPR(resultado, pr3, contPares);

			 times(&tempoDepois);
			cout << "Tempo_PR3 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
			//fprintf(stdout, "HV_PR3 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version3, lbound, ubound));


			times(&tempoAntes);
			BoundedParetoSet *newParetoSetAprox_version4 = rodaPR(resultado, pr4, contPares);

			 times(&tempoDepois);
			cout << "Tempo_PR4 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
			//fprintf(stdout, "HV_PR4 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version4, lbound, ubound));


			times(&tempoAntes);
			BoundedParetoSet *newParetoSetAprox_version6 = rodaPR(resultado, pr6, contPares);

			 times(&tempoDepois);
			cout << "Tempo_PR6 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
			//fprintf(stdout, "HV_PR6 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version6, lbound, ubound));



			 times(&tempoAntes);
			 BoundedParetoSet *newParetoSetAprox_version7 = rodaPR(resultado, pr7, contPares);

			times(&tempoDepois);
			cout << "Tempo_PR7 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
			//fprintf(stdout, "HV_PR7 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version7, lbound, ubound));

			// times(&tempoAntes);
			//  BoundedParetoSet *newParetoSetAprox_version8 = rodaPR(resultado, pr8, contPares);

			// times(&tempoDepois);
			// cout << "Tempo_PR8 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
			// //fprintf(stdout, "HV_PR8 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version8, lbound, ubound));


			times(&tempoAntes);
				Version5_D *v5_d = new Version5_D(resultado, subproblemas);
			    BoundedParetoSet *newParetoSetAprox_version5_msc = v5_d->MOPR_D(contPares);

			times(&tempoDepois);
			cout << "Tempo_msc = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
			cout <<"Pares_msc = "<<contPares<<endl;
			//fprintf(stdout, "HV_msc = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version5_msc, lbound, ubound));

			std::vector<double> lbound(n_objetivos), ubound(n_objetivos);
			for (int i = 0; i < n_objetivos; ++i){
				lbound[i] = 1e9;
				lbound[i] = 1e9;
				ubound[i] = -1e9;;
				ubound[i] = -1e9;;
			}
			
			for (list<Solution *>::iterator it= resultado->sol.begin() ; it != resultado->sol.end(); ++it) {
				for (int i = 0; i < n_objetivos; ++i){
					if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
					if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
				}
			}

			for (list<Solution *>::iterator it=  newParetoSetAprox_version1->sol.begin() ; it != newParetoSetAprox_version1->sol.end(); ++it) {
				for (int i = 0; i < n_objetivos; ++i){
					if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
					if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
				}
			}

			for (list<Solution *>::iterator it=  newParetoSetAprox_version2->sol.begin() ; it != newParetoSetAprox_version2->sol.end(); ++it) {
				for (int i = 0; i < n_objetivos; ++i){
					if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
					if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
				}
			}

			for (list<Solution *>::iterator it= newParetoSetAprox_version3->sol.begin() ; it != newParetoSetAprox_version3->sol.end(); ++it) {
				for (int i = 0; i < n_objetivos; ++i){
					if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
					if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
				}
			}
			for (list<Solution *>::iterator it= newParetoSetAprox_version4->sol.begin() ; it != newParetoSetAprox_version4->sol.end(); ++it) {
				for (int i = 0; i < n_objetivos; ++i){
					if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
					if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
				}
			}
			for (list<Solution *>::iterator it=  newParetoSetAprox_version6->sol.begin() ; it != newParetoSetAprox_version6->sol.end(); ++it) {
				for (int i = 0; i < n_objetivos; ++i){
					if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
					if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
				}
			}
			for (list<Solution *>::iterator it= newParetoSetAprox_version7->sol.begin() ; it != newParetoSetAprox_version7->sol.end(); ++it) {
				for (int i = 0; i < n_objetivos; ++i){
					if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
					if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
				}
			}
			for (list<Solution *>::iterator it= newParetoSetAprox_version5_msc->sol.begin() ; it != newParetoSetAprox_version5_msc->sol.end(); ++it) {
				for (int i = 0; i < n_objetivos; ++i){
					if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
					if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
				}
			}
			// for (list<Solution *>::iterator it= newParetoSetAprox_version8->sol.begin() ; it != newParetoSetAprox_version8->sol.end(); ++it) {
			// 	for (int i = 0; i < n_objetivos; ++i){
			// 		if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
			// 		if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
			// 	}
			// }

			fprintf(stdout, "HV_NSGAIII = %.9e\n", hyp->getHypervolume(resultado, lbound, ubound));
			fprintf(stdout, "HV_PR1 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version1, lbound, ubound));
			fprintf(stdout, "HV_PR2 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version2, lbound, ubound));
			fprintf(stdout, "HV_PR3 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version3, lbound, ubound));
			fprintf(stdout, "HV_PR4 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version4, lbound, ubound));
			fprintf(stdout, "HV_PR6 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version6, lbound, ubound));
			fprintf(stdout, "HV_PR7 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version7, lbound, ubound));
			// fprintf(stdout, "HV_PR8 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version8, lbound, ubound));
			fprintf(stdout, "HV_msc = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version5_msc, lbound, ubound));

			cout<<"\nResultado_NSGAIII:"<<endl;
			resultado->printSetPoints(stdout);


			cout<<"\nResultado_PR1:"<<endl;
			newParetoSetAprox_version1->printSetPoints(stdout);

			cout<<"\nResultado_PR2:"<<endl;
			newParetoSetAprox_version2->printSetPoints(stdout);

			cout<<"\nResultado_PR3:"<<endl;
			newParetoSetAprox_version3->printSetPoints(stdout);

			cout<<"\nResultado_PR4:"<<endl;
			newParetoSetAprox_version4->printSetPoints(stdout);

			cout<<"\nResultado_PR6:"<<endl;
			newParetoSetAprox_version6->printSetPoints(stdout);


				cout<<"\nResultado_PR7:"<<endl;
			newParetoSetAprox_version7->printSetPoints(stdout);

			// cout<<"\nResultado_PR8:"<<endl;
			// newParetoSetAprox_version8->printSetPoints(stdout);

			cout<<"\nResultado_PR_msc:"<<endl;
			newParetoSetAprox_version5_msc->printSetPoints(stdout);

	} else { ///////// ISSO foi colocado apenas para conter  o erro EXTRAORDINARIO que acontecia quando o NSGA-II retornava apenas uma solucaos
		cout << "Tempo_PR1_NOVO = 0.0 " <<endl;
			
		cout << "Tempo_PR2 = 0.0 " << endl;
					
		cout << "Tempo_PR3 = 0.0 " << endl;
					
		cout << "Tempo_PR4 = 0.0 "  << endl;
					
		cout << "Tempo_PR6 = 0.0 " << endl;
					
		cout << "Tempo_PR7 = 0.0 " << endl;
					
		// cout << "Tempo_PR8 = 0.0 " << endl;
					
		cout << "Tempo_msc = 0.0 " << endl;
			cout <<"Pares_msc = 0 "<<endl;

			std::vector<double> lbound(n_objetivos), ubound(n_objetivos);
			lbound[0] = 1e9;
			lbound[1] = 1e9;
			ubound[0] = -1e9;;
			ubound[1] = -1e9;;
			for (list<Solution *>::iterator i= resultado->sol.begin() ; i != resultado->sol.end(); ++i) {
				if ((*i)->getObj(0)<lbound[0]) lbound[0] = (*i)->getObj(0);
				if ((*i)->getObj(1)<lbound[1]) lbound[1] = (*i)->getObj(1);
				if ((*i)->getObj(0)>ubound[0]) ubound[0] = (*i)->getObj(0);
				if ((*i)->getObj(1)>ubound[1]) ubound[1] = (*i)->getObj(1);
			}
			
		fprintf(stdout, "HV_NSGAIII = %.9e\n", hyp->getHypervolume(resultado, lbound, ubound));
		fprintf(stdout, "HV_PR1 = %.9e\n", hyp->getHypervolume(resultado, lbound, ubound));
		fprintf(stdout, "HV_PR2 = %.9e\n", hyp->getHypervolume(resultado, lbound, ubound));
		fprintf(stdout, "HV_PR3 = %.9e\n", hyp->getHypervolume(resultado, lbound, ubound));
		fprintf(stdout, "HV_PR4 = %.9e\n", hyp->getHypervolume(resultado, lbound, ubound));
		fprintf(stdout, "HV_PR6 = %.9e\n", hyp->getHypervolume(resultado, lbound, ubound));
		fprintf(stdout, "HV_PR7 = %.9e\n", hyp->getHypervolume(resultado, lbound, ubound));
		// fprintf(stdout, "HV_PR8 = %.9e\n", hyp->getHypervolume(resultado, lbound, ubound));
		fprintf(stdout, "HV_msc = %.9e\n", hyp->getHypervolume(resultado, lbound, ubound));

			cout<<"\nResultado_NSGAIII:"<<endl;
			resultado->printSetPoints(stdout);


			cout<<"\nResultado_PR1:"<<endl;
			resultado->printSetPoints(stdout);

			cout<<"\nResultado_PR2:"<<endl;
			resultado->printSetPoints(stdout);

			cout<<"\nResultado_PR3:"<<endl;
			resultado->printSetPoints(stdout);

			cout<<"\nResultado_PR4:"<<endl;
			resultado->printSetPoints(stdout);

			cout<<"\nResultado_PR6:"<<endl;
			resultado->printSetPoints(stdout);


			cout<<"\nResultado_PR7:"<<endl;
			resultado->printSetPoints(stdout);

			// cout<<"\nResultado_PR8:"<<endl;
			// resultado->printSetPoints(stdout);

			cout<<"\nResultado_PR_msc:"<<endl;
			resultado->printSetPoints(stdout);
	}

	// OK OK OK OK
	// int contAval = 0;
	// Solution *s_1 = new Solution(n_locations,n_objetivos);
	// s_1->setAssignment_(0,1);
	// s_1->setAssignment_(1,2);
	// s_1->setAssignment_(2,7);
	// s_1->setAssignment_(3,9);
	// s_1->setAssignment_(4,6);
	// s_1->setAssignment_(5,5);
	// s_1->setAssignment_(6,0);
	// s_1->setAssignment_(7,4);
	// s_1->setAssignment_(8,3);
	// s_1->setAssignment_(9,8);
	// s_1->calculaObjetivos(contAval);
	// s_1->print_solution(stdout);
	
	return 0;
}
