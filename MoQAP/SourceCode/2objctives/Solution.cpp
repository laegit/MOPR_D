
#ifndef SOLUTION_CPP
#define SOLUTION_CPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <list>
#include <set>
#include <stack>
#include <string>
#include <utility>
#include <vector>
#include <iostream>

#include "random.cpp"
#include "param_geral.h"


using namespace std;


extern std::vector< std::vector<double> > distances;
// matriz de distâncias. 
// distances[i][j] == distânicia entre a localidade i e j

extern std::vector< vector< std::vector<double> > > flows;
//flows[obj][i][j] = o obj-ésimo fluxo entre i e j

extern int n_locations;
extern int n_objetivos;

class Solution
{
 public:
  // Default constructor
  Solution() {;} // ctor

  Solution(const size_t& size, const size_t &numObjetivos)
  : assignment_(size), objetivos_(numObjetivos), distance(0.0)
  {;} // ctor
  // Copy constructor
  explicit Solution(const Solution& rhs)
  : assignment_(rhs.assignment_), objetivos_(rhs.objetivos_), distance(rhs.distance)
  {;} // ctor
  
  // Destructor
  ~Solution()
  {;} // dtor

  // Assignment operator
  const void operator=(const Solution & rhs)
  {
    
    assignment_ = rhs.assignment_; // TODO: testar se, dado s1<-s2, editando s2 nao edita s1
    objetivos_ = rhs.objetivos_;
    rank = rhs.rank;
    distance = rhs.distance;
    posicaoListaNSGAII = rhs.posicaoListaNSGAII;
    
    // return *this; // Allow cascaded assignment.
  }

  const bool operator==(Solution &s) {

    for (size_t i=0; i<objetivos_.size(); i++){
      if (objetivos_[i]!=s.getObj(i)) return false;
    }
    return true;
  }

  // O PQA é de MINIMIZACAO
  bool operator>> (Solution &d) { // retorna true se this domina d
    bool greater = false;
    bool lesser = false;
    // Check if objetivos_ is lesser or greater than d.objetivos_ in any dimension.
    for (size_t i = 0; i < objetivos_.size(); ++i)
    {
      if (objetivos_.at(i) > d.getObj(i))
      {greater = true;}  // "this" perde em pelo menos um objetivo
      else if (objetivos_.at(i)  < d.getObj(i))
      {lesser = true;} // "this" ganha em pelo menos um objetivo
    }

    if (lesser) /// "this" ganha em pelo menos um objetivo
    {
      if (greater==false) // "this" nao perde em nenhum objetivo
        return true; // this domina d
      else return false;
    } else return false; // "this" nao ganha em nenhum objetivo
  }



  void print_solution(FILE *f){

    // fprintf(f,"Ponto objetivo = ");
    for (int i=0; i<objetivos_.size(); i++){
      fprintf(f,"%f ",objetivos_[i]);
    }

    fprintf(f, "\n");
  }

  double getObj(int i){
    return objetivos_[i];
  }

  void setObj(int index, double obj){
    objetivos_[index] = obj;
  }

  int getAssignment(int index){
    return assignment_[index];
  }

  std::vector<int> getAssignment_(){
    return assignment_;
  }

  void setAssignment_(int index, int facility){
    assignment_[index] = facility;
  }

  void calculaObjetivos(int &contAval){
    contAval++;
    for (int obj=0; obj<objetivos_.size(); obj++){ 
      objetivos_[obj] = 0;
      for (int i=0; i<assignment_.size(); i++){ // para cada localidade
        for (int j = 0; j < assignment_.size(); j++){ // para cada localidade
          objetivos_[obj]+=distances[assignment_[i]][assignment_[j]]*flows[obj][i][j]; // OK. Correto!
          //assignment_[i] = facilidade na localidade i
          //flows[obj][i][j] = o k-ésimo fluxo da facilidade i para a facilidade j 
        }
      }
    }
  }



  // One point crossover. Referência MURATA et al.  GENETIC ALGORITHMS FOR FLOWSHOP SCHEDULING PROBLEMS 
  void crossover (Solution *pai, Solution *mae, int &contAval){
    int left = pick_a_number(0,n_locations-1); // sorteia-se um ponto aleatorio
   
    std::vector<bool> auxx(assignment_.size()); // se o job i ja foi assigned para uma posicao, entao auxx[i]=true
    for (int i = 0; i < n_locations; ++i){ auxx[i] = false;}
    int jj = left-1;
    for (int i = 0; i < left; ++i){
      assignment_[i] = pai->getAssignment(jj--);
      auxx[assignment_[i]] = true;
    }
    // sobra exatamente right-left+1 jobs em mae que  nao estao em assignment_
    int cont = left; // de left ... right, inclusive
    for (int i = 0; i < n_locations; ++i){
      if (auxx[mae->getAssignment(i)]==false){
        assignment_[cont++] = mae->getAssignment(i);
        auxx[mae->getAssignment(i)] = true;
      }
    }
    calculaObjetivos(contAval);
  }


  // CX recombination de "Fitness Landscape Analysis and Memetic Algorithms for the Quadratic Assignment Problem"
  // troca dois alelos randômicos
  void crossoverCX(Solution *pai, Solution *mae, int &contAval){
    std::vector<int> auxx(assignment_.size()); // guarda localidades
    int limite = assignment_.size();
    for (int i=0; i<limite; i++) auxx[i] = i; // para cada localidade...
    int i=0;
    while(i<limite){
      if (pai->getAssignment(auxx[i])==mae->getAssignment(auxx[i])){ // se as localides aux[i] do pai e da mae forem iguais
        assignment_[auxx[i]] = pai->getAssignment(auxx[i]);
        for (int j=i+1; j<limite; j++) auxx[j-1] = auxx[j];
        limite--;
      } else i++;
    }
    while(limite>0){
      int pai_ou_mae =  pick_a_number(0,1); // 0 pai; 1 mae
      Solution *escolhida = pai;
      Solution *outra = mae; // apenas aponta para o mae
      if (pai_ou_mae == 1){
        escolhida = mae; // apenas aponta para o mae
        outra = pai; // apenas aponta para o pai
      }
      int indexLocation = pick_a_number(0,limite-1);
      int valorEscolhido = escolhida->getAssignment(auxx[indexLocation]);
      int valorEscolhidoInicial = valorEscolhido;
      assignment_[auxx[indexLocation]] = valorEscolhido;
      int buscarV = outra->getAssignment(auxx[indexLocation]);
      for (int j=indexLocation+1; j<limite; j++) auxx[j-1] = auxx[j];
      limite--;
      do{
        for (int j=0; j<limite; j++){ // for (I)
          if (escolhida->getAssignment(auxx[j])==buscarV){
            indexLocation = j;
            valorEscolhido = buscarV;
            assignment_[auxx[indexLocation]] = valorEscolhido;
            buscarV = outra->getAssignment(auxx[indexLocation]);
            break; // sai do for (I)
          }
        }
        for (int j=indexLocation+1; j<limite; j++) auxx[j-1] = auxx[j];
        limite--;
      } while(valorEscolhidoInicial!=buscarV);
    } 
    calculaObjetivos(contAval);
  }

  // mutaçao adatada de "Fitness Landscape Analysis and Memetic Algorithms for the Quadratic Assignment Problem"
  // troca dois alelos randômicos
  void mutation__a(Solution *sol, int &contAval){
    *this = *sol;
    int f_1 = pick_a_number(0,assignment_.size()-1);
    int f_2 = pick_a_number(0,assignment_.size()-1);
    int aux = assignment_[f_1];
    assignment_[f_1] = assignment_[f_2];
    assignment_[f_2] = aux;
    calculaObjetivos(contAval);
  }


  // insertion operation (baseado em "design-of-multiobjective-evolutionary-algorithms-application-to-...")
  void mutation(Solution *sol, int &contAval){
    *this = *sol;
    int f_1 = pick_a_number(0,assignment_.size()-1);
    int f_2 = pick_a_number(0,assignment_.size()-1);
    int init = f_1, fim=f_2;
    if (f_2<f_1){
       init = f_2; 
       fim = f_1;
    } 
    int aux = assignment_[fim];
    for (int i = fim; i > init; --i){
      assignment_[i] = assignment_[i-1];
    }
    assignment_[init] = aux;
    calculaObjetivos(contAval);
  }


  void  generateRandomSolution(int &contAval){
    std::vector<int> aux(assignment_.size());
    int limite = assignment_.size();
    for (int i=0; i<assignment_.size(); i++) aux[i] = i;

    for (int i=0; i<assignment_.size(); i++){
      int random = pick_a_number(0, limite-1);
      assignment_[i] = aux[random];
      for (int j=random+1; j<assignment_.size(); j++) aux[j-1] = aux[j];
      limite--;
    }
    calculaObjetivos(contAval);
  }

 private:
  
  // solution representation. 
  // CUIDADO: os índices sao as localidades. Os valores sao as facilities
  //assignment_[i] = j. Na localidade i temos a facilidade j
  // size = number of locations/facilities

  std::vector<int> assignment_;

  // objective vector
  // size = number of objectives
  std::vector<double> objetivos_;

  

  public:
   // @brief
  // crownd distance used by NSGAA-II 
  double distance; 
  int rank;
  int posicaoListaNSGAII; // guarda o index onde a soluçao é guardada na popupacao NUMPOPULACAO*2 do NSGA-II

}; // class Knapsack

#endif 
