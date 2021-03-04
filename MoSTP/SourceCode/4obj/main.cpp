/* 
=====================================================================================

Copyright Islame Felipe da COSTA FERNANDES (2018)
	
This code implements a NSGAIII+PLS
The data structure and some functions were kindly provided by Monteiro (2010)


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

#include "BoundedParetoSet.cpp"
#include "SolucaoEdgeSet.cpp"
#include "rmcPrim.cpp"
#include "popInicial.cpp"
#include "nsga3.cpp"
#include "LexKruskal.cpp"
#include "hyp_ind.cpp"
#include "Version1_NOVO.cpp"
#include "Version2.cpp"
#include "Version3.cpp"
#include "Version4.cpp"
#include "Version6.cpp"
#include "Version7.cpp"
#include "Version5_D.cpp"

#include <pthread.h>


using namespace std;


//OK: todos as particulas e seus membros acessam esse array. SOMENTE PARA ACESSAR O VALOR, nao há necessidade de lock()
double custos[NUMOBJETIVOS][NUMEROVERTICES][NUMEROVERTICES]; // global (extern) para todos os arquivos



// global extern  (nao precisa de mutex)
double lbound[NUMOBJETIVOS]; 
double ubound[NUMOBJETIVOS];


void input(){
  int n; // esta leitura de "n" é somente para cumprir o formato da instância. 
  cin>>n; // quantidade de vertices

  int org, dest;
  for (int i=0;i<NUMEROVERTICES-1;i++) {
    for (int j=i+1;j<NUMEROVERTICES;j++) {
      cin>>org;
      cin>>dest;
      if (org!=i) cout<<"ERRO Leitura 1"<<endl;
      if (dest!=j) cout<<"ERRO Leitura 2"<<endl;
      for (int ob = 0; ob<NUMOBJETIVOS; ob++){
        cin>>custos[ob][i][j];
        custos[ob][j][i] = custos[ob][i][j]; // a parte superior da matriz é refletida na inferior
      }
    }
  }
}


/*observacao: como estamos trabalhando com algoritmos heuristicos, e nao exatos, 
entao nao basta ordenar lexicograficamente pelo primeiro objetivo, minizar o primeiro objetivo e obter o segundo objetivo (como sendo o maximo).
Isso serviria para algoritmos exatos. Por exemplo:

(5, 100) --> miniza o primeiro (o exato obteria esse ponto. Aqui, 100 é o "upper bound" para o segundo objetivo)
(5, 105) --> miniza o primeiro (mas o heuristico pode obter esse ponto. Logo, 105 deve ser "upper bound" do segundo objetivo para o conjunto aproximativo em questao )*/
void getLowerAndUpperBound(){
	LexKruskal lex;
	SolucaoEdgeSet s1 (NUMEROVERTICES-1);
	SolucaoEdgeSet s2 (NUMEROVERTICES-1);
  SolucaoEdgeSet s3 (NUMEROVERTICES-1);
  SolucaoEdgeSet s4 (NUMEROVERTICES-1);
	std::vector<auxEdgeStruct> arestas;
	for (int i=0;i<NUMEROVERTICES-1;i++) {
		for (int j=i+1;j<NUMEROVERTICES;j++) {
			arestas.push_back((auxEdgeStruct){i,j,custos[0][i][j],custos[1][i][j],custos[2][i][j], custos[3][i][j],0});
		}
	} 
	int contlxo;
	std::sort(arestas.begin(), arestas.begin()+arestas.size(), lexicografica_obj1);
	lex.executar(s1,arestas, contlxo);; // minimiza o primeiro objetivo 
	std::sort(arestas.begin(), arestas.begin()+arestas.size(), lexicografica_obj2);
	lex.executar(s2,arestas, contlxo);; // minimiza o segundo objetivo 
  std::sort(arestas.begin(), arestas.begin()+arestas.size(), lexicografica_obj3);
  lex.executar(s3,arestas, contlxo);; // minimiza o terceiro objetivo 
  std::sort(arestas.begin(), arestas.begin()+arestas.size(), lexicografica_obj4);
  lex.executar(s4,arestas, contlxo);; // minimiza o terceiro objetivo 

  // s1.isTree();
  //  s2.isTree();
  //   s3.isTree();
  //    s4.isTree();

	cout<<"lower_bound = "<<s1.getObj(0)<<" "<<s2.getObj(1)<<" "<<s3.getObj(2)<<" "<<s4.getObj(3)<<endl;
	lbound[0] = s1.getObj(0);
	lbound[1] = s2.getObj(1);
  lbound[2] = s3.getObj(2);
  lbound[3] = s4.getObj(3);

	std::sort(arestas.begin(), arestas.begin()+arestas.size(), lexicografica_max_obj1);
	lex.executar(s1,arestas, contlxo);; // maximiza o primeiro obejtivo
	std::sort(arestas.begin(), arestas.begin()+arestas.size(), lexicografica_max_obj2);
	lex.executar(s2,arestas, contlxo);; // maximiza o segundo obejtivo
  std::sort(arestas.begin(), arestas.begin()+arestas.size(), lexicografica_max_obj3);
  lex.executar(s3,arestas, contlxo);; // maximiza o terceiro obejtivo
  std::sort(arestas.begin(), arestas.begin()+arestas.size(), lexicografica_max_obj4);
  lex.executar(s4,arestas, contlxo);; // maximiza o terceiro obejtivo

  // s1.isTree();
  //  s2.isTree();
  //   s3.isTree();
	 //   s4.isTree();
	cout<<"upper_bound = "<<s1.getObj(0)<<" "<<s2.getObj(1)<<" "<<s3.getObj(2)<<" "<<s4.getObj(3)<<endl;
	ubound[0] = s1.getObj(0);
	ubound[1] = s2.getObj(1);
  ubound[2] = s3.getObj(2);
  ubound[3] = s4.getObj(3);
}



// Ver Arroyo Vieira e Vianna (2008)
void calculaVetoresLambda_SISTEMATICO(double lambda[][NUMOBJETIVOS], int s){
    int cont =0;
    if (NUMOBJETIVOS==2){ // s-1 
        // para atingir 10^6 vetores, necessario s=10^6-1
        double l1, l2;
        for (int i=0; i<=s; i++){
            l1 = i;
            l2 = s - l1;
            lambda[cont][0] = static_cast<double>(l1/s);
            lambda[cont][1] = static_cast<double>(l2/s);
            cont++;
            // std::vector<double> lam = {static_cast<double>(l1/s), static_cast<double>(l2/s)};
            // lambda.push_back(lam);
        }

    } else if (NUMOBJETIVOS==3){ // quantidade vetores:  (s+1)(s+2)/2
        // s = 1417 e s=1418 para obter 10^6 vetores
        double l1, l2, l3;
        for (int i=0; i<=s && cont<NUMSUMPROBLEMAS; i++){
            l1 = i;
            for (int j=0; j<=s-i && cont<NUMSUMPROBLEMAS; j++){
                l2 = j;
                l3 = s-i - l2;
                lambda[cont][0] = static_cast<double>(l1/s);
                lambda[cont][1] = static_cast<double>(l2/s);
                lambda[cont][2] = static_cast<double>(l3/s);
                cont++;
            }
        }

    } else if (NUMOBJETIVOS==4){ // quantidade vetores: 1/6 (s^3 + 6 s^2 + 11s + 6)
        // s = 179 ou s=180 para obter 10^6 vetores
        double l1, l2, l3, l4;
        for (int i=0; i<=s && cont<NUMSUMPROBLEMAS; i++){
            l1 = i;
            for (int j=0; j<=s-i && cont<NUMSUMPROBLEMAS; j++){
                l2 = j;
                for (int z=0; z<=s-i-j && cont<NUMSUMPROBLEMAS; z++){
                    l3 = z;
                    l4 = s-i-j - l3;
                    lambda[cont][0] = static_cast<double>(l1/s);
                    lambda[cont][1] = static_cast<double>(l2/s);
                    lambda[cont][2] = static_cast<double>(l3/s);
                    lambda[cont][3] = static_cast<double>(l4/s);
                    cont++;
                }
            }
        }
    } 
}

// pega o melhor valor possivel de cada objetivo
void getPontoIdeal(double ideal1[NUMOBJETIVOS]){
    LexKruskal lex;
    SolucaoEdgeSet s1 (NUMEROVERTICES-1);
    std::vector<auxEdgeStruct> arestas;
    for (int i=0;i<NUMEROVERTICES-1;i++) {
        for (int j=i+1;j<NUMEROVERTICES;j++) {
            arestas.push_back((auxEdgeStruct){i,j,custos[0][i][j],custos[1][i][j],custos[2][i][j],custos[3][i][j],0});
        }
    } 
    int contlxo;
    std::sort(arestas.begin(), arestas.begin()+arestas.size(), lexicografica_obj1);
    lex.executar(s1,arestas, contlxo);; // minimiza o primeiro objetivo 
    ideal1[0] = s1.getObj(0);
    
    std::sort(arestas.begin(), arestas.begin()+arestas.size(), lexicografica_obj2);
    lex.executar(s1,arestas, contlxo);; // minimiza o segundo objetivo 
    ideal1[1] = s1.getObj(1);
    
    std::sort(arestas.begin(), arestas.begin()+arestas.size(), lexicografica_obj3);
    lex.executar(s1,arestas, contlxo);; // minimiza o segundo objetivo 
    ideal1[2] = s1.getObj(2);
    
    std::sort(arestas.begin(), arestas.begin()+arestas.size(), lexicografica_obj4);
    lex.executar(s1,arestas, contlxo);; // minimiza o segundo objetivo 
    ideal1[3] = s1.getObj(3);
   
}








/* --------------------------------------------- */
/* -------- FUNCAO *NOVA* PARA O PR ---------- */
/* --------------------------------------------- */

BoundedParetoSet * rodaPR(BoundedParetoSet * conjAprox, PathRelinking *pr, int &contPares){
     int aval = 0;
     contPares = 0;
     BoundedParetoSet *retorno = new BoundedParetoSet();

     // copia previamente
     for (list<SolucaoEdgeSet *>::iterator it_1=conjAprox->sol.begin(); it_1!=conjAprox->sol.end(); it_1++)
        retorno->adicionarSol(*it_1); // copia, previamente, o conjunto inicial.
            
     // O PR(o,d) adiciona 'o' no conjunto. 
     // Como nós chamamos PR(o,d) e PR(d,o), entao 
     // todos os itens de conjAprox sao adicionados no retorno.

      while(aval<NUM_AVALIACOES){
        list<SolucaoEdgeSet *>::iterator it_1=conjAprox->sol.begin();
         while (it_1!=conjAprox->sol.end() && aval<NUM_AVALIACOES) {
            list<SolucaoEdgeSet *>::iterator it_2=conjAprox->sol.begin();
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
      // cout<<"aval = "<<aval<<endl;
     return retorno;
} 


// semente paretofile tempoCPU_USER tempoREAL hypervolume < instância
int main(int argc, char *argv[]){


	struct tms tempoAntes, tempoDepois;

	int seemente = std::atoi(argv[1]);
	init_genrand64(seemente);
	cout<<"========= Estatisticas ========= "<<endl;
	cout<<"Semente: "<<seemente<<endl;
	input(); // ler instância
	cout<<"Instância lida..."<<endl;
	getLowerAndUpperBound();


	times(&tempoAntes);

    

  BoundedParetoSet *resultado = new BoundedParetoSet();
  nsga3 *n3 = new nsga3();
  SolucaoEdgeSet **population = n3->executar();


  times(&tempoDepois);

  // RESULTADO: apenas as nao-dominadas (primeiro front) da populacao
  for (int iii=0; iii<TAMANHOPOPULACAO_NSGAII; iii++){
      resultado->adicionarSol(population[iii]);
  }
  fprintf(stdout,"Tempo_NSGA3 = %.2lf\n", (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.0 );
  fprintf(stdout,"Quantidade_solucoes_NSGA3 = %d\n", resultado->getSize());

  Hypervolume *hyp = new Hypervolume();
  // double hyper=hyp->getHypervolume(resultado, lbound, ubound);


  // fprintf(stdout, "Hypervolume_NSGA3 = %.9e\n", hyper);
  // double ideal[NUMOBJETIVOS];
 //    // getPontoIdeal(conjAprox, ideal);
 //    getPontoIdeal(ideal);

    double lambda_pr[NUMOBJETIVOS];
    getLambdaRandom(lambda_pr);
    double subproblemas[NUMSUMPROBLEMAS][NUMOBJETIVOS];
    calculaVetoresLambda_SISTEMATICO(subproblemas, 8); // NUMSUMPROBLEMAS-1);
     // begin coments
  PathRelinking **pr_versions = new PathRelinking*[7]; 

    pr_versions[0] = new Version1_NOVO();
    pr_versions[1] = new Version2(IRandom(0,NUMOBJETIVOS-1)); 
    pr_versions[2] = new Version3(); 
    pr_versions[3] = new Version4(lambda_pr); 
    pr_versions[5] = new Version6(); 
    pr_versions[6] = new Version7(); 


   int contPares = 0;
   times(&tempoAntes);
        BoundedParetoSet *newParetoSetAprox_version1 = rodaPR(resultado, pr_versions[0], contPares);
    times(&tempoDepois);
    cout << "Tempo_PR1_NOVO = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
    // fprintf(tempoPR1, "%f\n",  (double) (tempoDepois-tempoAntes) / 100.00);
    cout<<"Pares_PR1_NOVO = "<<contPares<<endl;
    //fprintf(stdout, "HV_PR1_NOVO = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version1, lbound, ubound));


    times(&tempoAntes);

        BoundedParetoSet *newParetoSetAprox_version2 = rodaPR(resultado, pr_versions[1], contPares);

    times(&tempoDepois);
    cout << "Tempo_PR2 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
    // fprintf(tempoPR2, "%f\n",  (double) (tempoDepois-tempoAntes) / 100.00);
    cout<<"Pares_PR2 = "<<contPares<<endl;
    //fprintf(stdout, "HV_PR2 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version2, lbound, ubound));


    times(&tempoAntes);

        BoundedParetoSet *newParetoSetAprox_version3 = rodaPR(resultado, pr_versions[2], contPares);

    times(&tempoDepois);
    cout << "Tempo_PR3 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
    // fprintf(tempoPR3, "%f\n",  (double) (tempoDepois-tempoAntes) / 100.00);
    cout<<"Pares_PR3 = "<<contPares<<endl;
    //fprintf(stdout, "HV_PR3 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version3, lbound, ubound));


    times(&tempoAntes);

        BoundedParetoSet *newParetoSetAprox_version4 = rodaPR(resultado, pr_versions[3], contPares);

    times(&tempoDepois);
    cout << "Tempo_PR4 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
    // fprintf(tempoPR4, "%f\n",  (double) (tempoDepois-tempoAntes) / 100.00);
    cout <<"Pares_PR4 = "<<contPares<<endl;
    //fprintf(stdout, "HV_PR4 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version4, lbound, ubound));
    // aqui fim
 

   times(&tempoAntes);
       Version5_D *v5_d = new Version5_D(resultado, subproblemas);
    // times(&tempoDepois);
    // cout << "Tempo_msc_init = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
    // times(&tempoAntes);
       BoundedParetoSet *newParetoSetAprox_version5_msc = v5_d->MOPR_D(contPares);
    times(&tempoDepois);
    cout << "Tempo_msc = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
    cout <<"Pares_msc = "<<contPares<<endl;
   // fprintf(stdout, "HV_msc = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version5_msc, lbound, ubound));
    


    times(&tempoAntes);

        BoundedParetoSet *newParetoSetAprox_version6 = rodaPR(resultado, pr_versions[5], contPares);

    times(&tempoDepois);
    cout << "Tempo_PR6 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
    cout<<"Pares_PR6 = "<<contPares<<endl;
   // fprintf(stdout, "HV_PR6 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version6, lbound, ubound));


    times(&tempoAntes);

        BoundedParetoSet *newParetoSetAprox_version7 = rodaPR(resultado, pr_versions[6], contPares);

    times(&tempoDepois);
    cout << "Tempo_PR7 = " << (double) (tempoDepois.tms_utime - tempoAntes.tms_utime) / 100.00 << endl;
    cout<<"Pares_PR7 = "<<contPares<<endl;
 //   fprintf(stdout, "HV_PR7 = %.9e\n", hyp->getHypervolume(newParetoSetAprox_version7, lbound, ubound));


    std::vector<double> lbound(NUMOBJETIVOS), ubound(NUMOBJETIVOS);
      for (int i = 0; i < NUMOBJETIVOS; ++i){
        lbound[i] = 1e9;
        lbound[i] = 1e9;
        ubound[i] = -1e9;;
        ubound[i] = -1e9;;
      }
      
      for (list<SolucaoEdgeSet *>::iterator it= resultado->sol.begin() ; it != resultado->sol.end(); ++it) {
        for (int i = 0; i < NUMOBJETIVOS; ++i){
          if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
          if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
        }
      }

      for (list<SolucaoEdgeSet *>::iterator it=  newParetoSetAprox_version1->sol.begin() ; it != newParetoSetAprox_version1->sol.end(); ++it) {
        for (int i = 0; i < NUMOBJETIVOS; ++i){
          if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
          if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
        }
      }

      for (list<SolucaoEdgeSet *>::iterator it=  newParetoSetAprox_version2->sol.begin() ; it != newParetoSetAprox_version2->sol.end(); ++it) {
        for (int i = 0; i < NUMOBJETIVOS; ++i){
          if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
          if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
        }
      }

      for (list<SolucaoEdgeSet *>::iterator it= newParetoSetAprox_version3->sol.begin() ; it != newParetoSetAprox_version3->sol.end(); ++it) {
        for (int i = 0; i < NUMOBJETIVOS; ++i){
          if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
          if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
        }
      }
      for (list<SolucaoEdgeSet *>::iterator it= newParetoSetAprox_version4->sol.begin() ; it != newParetoSetAprox_version4->sol.end(); ++it) {
        for (int i = 0; i < NUMOBJETIVOS; ++i){
          if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
          if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
        }
      }
      for (list<SolucaoEdgeSet *>::iterator it=  newParetoSetAprox_version6->sol.begin() ; it != newParetoSetAprox_version6->sol.end(); ++it) {
        for (int i = 0; i < NUMOBJETIVOS; ++i){
          if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
          if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
        }
      }
      for (list<SolucaoEdgeSet *>::iterator it= newParetoSetAprox_version7->sol.begin() ; it != newParetoSetAprox_version7->sol.end(); ++it) {
        for (int i = 0; i < NUMOBJETIVOS; ++i){
          if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
          if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
        }
      }
      for (list<SolucaoEdgeSet *>::iterator it= newParetoSetAprox_version5_msc->sol.begin() ; it != newParetoSetAprox_version5_msc->sol.end(); ++it) {
        for (int i = 0; i < NUMOBJETIVOS; ++i){
          if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
          if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
        }
      }
      // for (list<Solution *>::iterator it= newParetoSetAprox_version8->sol.begin() ; it != newParetoSetAprox_version8->sol.end(); ++it) {
      //  for (int i = 0; i < NUMOBJETIVOS; ++i){
      //    if ((*it)->getObj(i)<lbound[i]) lbound[i] = (*it)->getObj(i);
      //    if ((*it)->getObj(i)>ubound[i]) ubound[i] = (*it)->getObj(i);
      //  }
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

    cout<<"\nResultado_NSGA3:"<<endl;
    resultado->printSetPoints(stdout);


    cout<<"\nResultado_PR1_NOVO:"<<endl;
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

     cout<<"\nResultado_PR_msc:"<<endl;
    newParetoSetAprox_version5_msc->printSetPoints(stdout);


      

    return 0;

}