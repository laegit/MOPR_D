#include "binary_heap.cpp"
#include "momkp.cpp"
#include "random.cpp"
#include "BoundedParetoSet.cpp"
#include "Knapsack.cpp"
#include "nsga3.cpp"
// #include "EpislonApproximatePaetoSet.cpp"
// #include "ListaSimples.cpp"
#include "Version1.cpp"
#include "Version2.cpp"
#include "Version3.cpp"
#include "Version4.cpp"
#include "Version6.cpp"
#include "Version7.cpp"
//#include "Version8.cpp"
#include "Version5_D.cpp"
#include "hyp_ind.cpp"


/*


 Toda vez que uma soluçao é criada ou alterada, precisa-se contabilizar a quantidade de avaliaçoes da FO
 Uma nova soluçao é alterada quando os itens da mochila sao modificados


*/


MOMKP *my_momkp;

// para a mochila, o ponto utopico consiste numa
// soluçao que possui todos os itens disponiveis, 
// sem considerar o limite da mochila
void getPontoUtopico(std::vector<double> &ubound, std::vector<double> &lbound){
	for (int j=0; j<my_momkp->dimension_; j++){
		int sum_j = 0;
		for (int i=0; i<my_momkp->size_; i++){
			sum_j+=my_momkp->weight_.at(i).at(j);
		}
		ubound.at(j) = sum_j;
		lbound.at(j) = 0.0;
	}
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
/* ---------- FUNCAO *NOVA* PARA O PR ---------- */
/* --------------------------------------------- */

BoundedParetoSet * rodaPR(BoundedParetoSet * conjAprox, PathRelinking *pr, int &contPares){
    
     int aval = 0;
     BoundedParetoSet *retorno = new BoundedParetoSet();
     // copia previamente
     for (list<Knapsack *>::iterator it_1=conjAprox->sol.begin(); it_1!=conjAprox->sol.end(); it_1++)
        retorno->adicionarSol(*it_1); // copia, previamente, o conjunto inicial.
       
     // O PR(o,d) adiciona 'o' no conjunto. 
     // Como nós chamamos PR(o,d) e PR(d,o), entao 
     // todos os itens de conjAprox sao adicionados no retorno.
     contPares =0;
     // fprintf(stdout, "\t HV = %.9e\n", hyp->getHypervolume(retorno, lbound, ubound));

	while(aval<NUM_AVALIACOES){
		list<Knapsack *>::iterator it_1=conjAprox->sol.begin();
	     while (it_1!=conjAprox->sol.end()  && aval<NUM_AVALIACOES) {
	        list<Knapsack *>::iterator it_2=conjAprox->sol.begin();
	        while (it_2!=conjAprox->sol.end()  && aval<NUM_AVALIACOES){
	            
	            if (it_1!=it_2){
	            	contPares++;
	                pr->PR(retorno, *it_1, *it_2, aval);
	            }

	            it_2++;
	        }
	        it_1++;
	     }
	     //cout<<"aval = "<<aval<<" ---  contPares = "<<contPares<<endl;
	 }

     return retorno;
} 




// 1 seed 2 classe 3 NumItens 4 NumObjetivos
int main(int argc, char const *argv[])
{
	struct tms tempoAntes, tempoDepois;


	int seed = atoi(argv[1]); // seed
	randomize(seed);
	std::srand(seed);
	string classe = argv[2]; // classe
	int NumItens = atoi(argv[3]); // NumItens
	int NumObjetivos = atoi(argv[4]); // NumObjetivos
	// int tam = (int)(1<<(8*3));
	// cout<<tam<<endl;
	my_momkp = new MOMKP(classe, NumItens, NumObjetivos);
	cout<<"Semente = "<<seed<<endl;
	cout<<"Class = "<<classe<<endl;
	cout<<"Size = "<<my_momkp->size_<<endl;
	cout<<"Dimension = "<<my_momkp->dimension_<<endl;



	PathRelinking **pr_versions = new PathRelinking*[8]; 
	int contPares = 0;
    pr_versions[0] = new Version1();
    pr_versions[1] = new Version2(pick_a_number(0,my_momkp->dimension_-1)); 
    pr_versions[2] = new Version3(1); // metodo = 1 ==> otimizar os objetivos 1,2,3... na sequência. Por exemplo, primeiro elemento otimiza objetivo 1, segundo elemento, objetivo 2 e assim por diante...
	std::vector<double> lambdaAux(my_momkp->dimension_, 0.0);
	my_momkp->getLambdaRandom(lambdaAux);
    pr_versions[3] = new Version4(lambdaAux); 

		pr_versions[5] = new Version6();
		pr_versions[6] = new Version7();
		// pr_versions[7] = new Version8();


    std::vector<double> subproblemas[NUMSUMPROBLEMAS];
    if (my_momkp->dimension_==2){
    	calculaVetoresLambda_SISTEMATICO(subproblemas, my_momkp->dimension_,NUMSUMPROBLEMAS-1); 
    } else if (my_momkp->dimension_==3){
    	calculaVetoresLambda_SISTEMATICO(subproblemas, my_momkp->dimension_,11); 
    } else if (my_momkp->dimension_==4){
    	calculaVetoresLambda_SISTEMATICO(subproblemas, my_momkp->dimension_,6); 
    } else if (my_momkp->dimension_==5){
    	calculaVetoresLambda_SISTEMATICO(subproblemas, my_momkp->dimension_,5); 
    }
    

    Hypervolume *hyp = new Hypervolume();
    std::vector<double> ubound(my_momkp->dimension_);
	std::vector<double> lbound(my_momkp->dimension_);
	getPontoUtopico(ubound, lbound);

	times(&tempoAntes);

	BoundedParetoSet *resultado = new BoundedParetoSet();
	nsga3 *nii = new nsga3();
	Knapsack **population = nii->executar();
	
	// RESULTADO: apenas as nao-dominadas (primeiro front) da populacao
	for (int iii=0; iii<TAMANHOPOPULACAO_NSGAII; iii++){
		resultado->adicionarSol(population[iii]);
		
	}

	times(&tempoDepois);
	fprintf(stdout,"Tempo_NSGA3 = %.2lf\n", (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.0 );
	double hyper=hyp->getHypervolume(resultado, lbound, ubound);
	fprintf(stdout, "Hypervolume_NSGA3 = %.9e\n", hyper);

   times(&tempoAntes);
        BoundedParetoSet *newParetoSetAprox_version1 = rodaPR(resultado, pr_versions[0], contPares);
    times(&tempoDepois);
    cout << "Tempo_PR1 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
    cout<<"Pares_PR1 = "<<contPares<<endl;
    fprintf(stdout, "HV_PR1 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version1, lbound, ubound));


    times(&tempoAntes);

        BoundedParetoSet *newParetoSetAprox_version2 = rodaPR(resultado, pr_versions[1], contPares);

    times(&tempoDepois);
    cout << "Tempo_PR2 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
    cout<<"Pares_PR2 = "<<contPares<<endl;
    fprintf(stdout, "HV_PR2 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version2, lbound, ubound));


    times(&tempoAntes);

        BoundedParetoSet *newParetoSetAprox_version3 = rodaPR(resultado, pr_versions[2], contPares);

    times(&tempoDepois);
    cout << "Tempo_PR3 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
    cout<<"Pares_PR3 = "<<contPares<<endl;
    fprintf(stdout, "HV_PR3 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version3, lbound, ubound));


    times(&tempoAntes);

        BoundedParetoSet *newParetoSetAprox_version4 = rodaPR(resultado, pr_versions[3], contPares);

    times(&tempoDepois);
    cout << "Tempo_PR4 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
    cout <<"Pares_PR4 = "<<contPares<<endl;
    fprintf(stdout, "HV_PR4 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version4, lbound, ubound));
    
 
	
	
	times(&tempoAntes);

        BoundedParetoSet *newParetoSetAprox_version6 = rodaPR(resultado, pr_versions[5], contPares);

    times(&tempoDepois);
    cout << "Tempo_PR6 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
    cout <<"Pares_PR6 = "<<contPares<<endl;
    fprintf(stdout, "HV_PR6 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version6, lbound, ubound));

    

	times(&tempoAntes);

        BoundedParetoSet *newParetoSetAprox_version7 = rodaPR(resultado, pr_versions[6], contPares);

    times(&tempoDepois);
    cout << "Tempo_PR7 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
    cout <<"Pares_PR7 = "<<contPares<<endl;
    fprintf(stdout, "HV_PR7 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version7, lbound, ubound));


    // BoundedParetoSet *newParetoSetAprox_version8;
    // if (my_momkp->dimension_==2){
    // 	times(&tempoAntes);

    //     	newParetoSetAprox_version8 = rodaPR(resultado, pr_versions[7], contPares);

	   //  times(&tempoDepois);
	   //  cout << "Tempo_PR8 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
	   //  cout <<"Pares_PR8 = "<<contPares<<endl;
	   //  fprintf(stdout, "HV_PR8 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version8, lbound, ubound));
    // }

    times(&tempoAntes);
		Version5_D *v5_d = new Version5_D(resultado, subproblemas);
	    BoundedParetoSet *newParetoSetAprox_version5_msc = v5_d->MOPR_D(contPares);

	times(&tempoDepois);
	cout << "Tempo_msc = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
	cout <<"Pares_msc = "<<contPares<<endl;
	fprintf(stdout, "HV_msc = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version5_msc, lbound, ubound));
    

	cout<<"\nResultado_NSGA3:"<<endl;
	resultado->printSetPoints(stdout);



    cout<<"\nResultado_PR1:"<<endl;
    newParetoSetAprox_version1->printSetPoints(stdout);

    cout<<"\nResultado_PR2:"<<endl;
    newParetoSetAprox_version2->printSetPoints(stdout);


    cout<<"\nResultado_PR3:"<<endl;
    newParetoSetAprox_version3->printSetPoints(stdout);

    cout<<"\nResultado_PR4:"<<endl;
    newParetoSetAprox_version4->printSetPoints(stdout);

	cout<<"\nResultado_PR_msc:"<<endl;
    newParetoSetAprox_version5_msc->printSetPoints(stdout);

	cout<<"\nResultado_PR6:"<<endl;
 	newParetoSetAprox_version6->printSetPoints(stdout);

 	cout<<"\nResultado_PR7:"<<endl;
 	newParetoSetAprox_version7->printSetPoints(stdout);

 // 	if (my_momkp->dimension_==2){
 // 		cout<<"\nResultado_PR8:"<<endl;
 // 	  	newParetoSetAprox_version8->printSetPoints(stdout);
	// }

	return 0;
}