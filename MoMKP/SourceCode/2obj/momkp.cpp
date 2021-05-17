// @file momkp.cpp
// @author Denis Felipe
// @version 0.1
//
// @brief
//   Multiobjective Multidimensional Knapsack problem solver.
//
// @example
//   #include <iostream>
//   #include "momst.h"
//
//   int main(int argc, char** argv)
//   {
//     algorithm::MOMKP momkp("ZMKP", 250, 2);
//     std::cout << "Runtime: " << momkp.mosca_d("pareto_approximation_file.sol") << "s" << std::endl;
//
//     return 0;
//   } // function main
//
// @section LICENSE
// Copyright (c) 2017, Denis Felipe
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met: 
// 
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution. 
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// The views and conclusions contained in the software and documentation are those
// of the authors and should not be interpreted as representing official policies, 
// either expressed or implied, of the FreeBSD Project.

// #define guard to prevent multiple inclusion.
#ifndef MOMKP_CPP
#define MOMKP_CPP

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

#include "binary_heap.cpp"
#include "random.cpp"
#include "Knapsack.cpp"

using namespace std;


 inline bool decrescente(pair<double, int> a, pair<double, int> b){
    return a.first>b.first; // ordem decrescente
  }

class MOMKP
{
 public:

  // Constructor
  // @param string $type
  //   Source of instances.
  // @param size_t $size
  //   Number of vertices.
  // @param size_t $dimension
  //   Number of objectives.
  MOMKP(const std::string& type, const size_t& size, const size_t& dimension)
  : profit_(size), weight_(size), capacity_(dimension), type_(type), size_(size), dimension_(dimension)
  {
    for (size_t i = 0; i < size_; ++i)
    {
      profit_.at(i).resize(dimension_);
      weight_.at(i).resize(dimension_);
    }

    const std::string instance_file = "Casos de Teste/" + type_ + "/" + std::to_string(size_) + "_" + std::to_string(dimension_) + ".txt";
    std::ifstream in_file;
    in_file.exceptions(std::ifstream::badbit);
    try
    {
      in_file.open(instance_file);

      for (size_t j = 0; j < dimension_; ++j)
      {
        in_file.ignore(std::numeric_limits<std::streamsize>::max(), ':');
        in_file.ignore(std::numeric_limits<std::streamsize>::max(), ':');
        in_file >> capacity_.at(j);

        for (size_t i = 0; i < size_; ++i)
        {
          in_file.ignore(std::numeric_limits<std::streamsize>::max(), ':');
          in_file.ignore(std::numeric_limits<std::streamsize>::max(), ':');
          in_file >> weight_.at(i).at(j);
          in_file.ignore(std::numeric_limits<std::streamsize>::max(), ':');
          in_file >> profit_.at(i).at(j);
        }
      }
      in_file.close();
    }
    catch (const std::ifstream::failure& e)
    {
      std::cerr << "Error reading file at location [" << instance_file << "]." << std::endl;
      std::exit(1);
    }

    // calcula upper_bound
    for (int j=0; j<dimension_; j++){
      int sum_j = 0;
      for (int i=0; i<size_; i++){
        sum_j+=weight_.at(i).at(j);
      }
      upper_bound.push_back(sum_j);
      lower_bound.push_back(0);
    }


  } // ctor
  // Destructor
  ~MOMKP()
  {;} // dtor

  
 public:
  
  // @brief
  //   Ensure a solution feasability.
  // @param Knapsack $knapsack
  //   Knapsack to be repaired.
  // @param vector<bool> $theme
  //   Theme of research.
  // @param vector<pair<double, size_t>> $sorted_item
  //   Preference of items.
  void repair(Knapsack *knapsack, std::vector<std::pair<double, size_t>>& sorted_item, int &intContAval)
  {
    size_t tmp = 0;
    ++intContAval;  // UMA AVALICAO DA FO

    // Try to remove the worst items from the knapsack.
    for (size_t i = 0; knapsack->overload(capacity_) && i < size_; ++i)
    {
      tmp = sorted_item.at(i).second;
      if (knapsack->at(tmp))
      {knapsack->remove_item(tmp, profit_.at(tmp), weight_.at(tmp));}
    }
   
  } // function repair
  


  void getLambdaRandom(std::vector<double> &lambdaAux){
    double sumL = 0.0;
    for (int i = 0; i < dimension_-1; ++i){
      double l = pick_a_number(0.0, 1.0-sumL); //Random(0.0,1.0-l1);
      sumL+=l;
      lambdaAux[i] = l;
    }
    lambdaAux[dimension_-1] = 1.0 - sumL;
  }


  void preferenceItens(std::vector<std::pair<double, size_t>>& sorted_item, Knapsack *solution){
      // Sort items by prefference.
        double ratio_a;
        double ratio_b;
        std::vector<double> lambdaAux(dimension_,0.0);
        getLambdaRandom(lambdaAux);
        for (size_t i = 0; i < size_; ++i)
        {
          ratio_a = 0.0;
          ratio_b = 0.0;
          for (size_t j = 0; j < dimension_; ++j)
          {
            ratio_a += lambdaAux.at(j) * profit_.at(i).at(j);
            ratio_b += static_cast<double>(weight_.at(i).at(j)) / (1.01 * capacity_.at(j) - solution->weight(j));
          }
          sorted_item.at(i).first = ratio_a / ratio_b;
          sorted_item.at(i).second = i;
        }
        sort(sorted_item.begin(), sorted_item.end());

    }

  // crossover utilizado por: Lust et Teghem (2008) "MEMOTS: A MEMETIC ALGORITHM INTEGRATING TABU SEARCH FOR COMBINATORIAL MULTIOBJECTIVE OPTIMIZATION"
   Knapsack * crossover(Knapsack *pai, Knapsack *mae, int &intContAval){

      Knapsack *filho = new Knapsack(size_, dimension_);
      for (int i=0; i<size_/2; i++){
        if (pai->at(i)==true){
            filho->insert_item(i,  profit_.at(i), weight_.at(i));
        }
      }
      for (int i=size_/2; i<size_; ++i)
      {
        if (mae->at(i)==true){
          filho->insert_item(i,  profit_.at(i), weight_.at(i));
        }
      }

      if (filho->overload(capacity_)){
        std::vector<std::pair<double, size_t>> sorted_item(size_);
        preferenceItens(sorted_item, filho);
        repair(filho, sorted_item, intContAval);
      } else  intContAval++;
      return filho;
   }
  
 
   void mutation(Knapsack *solution, int &intContAval){
      int it = pick_a_number(0,size_-1);
      if (solution->at(it)==true){
        solution->remove_item(it, profit_.at(it), weight_.at(it)); // nao precisa de reparo
        intContAval++;
      }else{
        solution->insert_item(it, profit_.at(it), weight_.at(it));
        
        if (solution->overload(capacity_)){
          std::vector<std::pair<double, size_t>> sorted_item(size_);
          preferenceItens(sorted_item, solution);
          repair(solution, sorted_item, intContAval);
        }
      }
   }


 

  // baseado no algoritmo de Vianna e Arroyo (2004) 
  //alfa deve ser entre 0 e 1
  // Quanto mais alfa for proximo de 1, mais randomizado é a fase construtiva
  // quanto mais proximo de 0, mais guloso 
  void BuildSolution(Knapsack *solution, double alfa, int &contAval){
    contAval++;
    std::vector< pair<double, int> > CL;
    for (size_t i = 0; i < size_; ++i){// para cada item
      if (solution->at(i)==false){//se esse item ainda nao está contido na soluçao
        double ratio_a = 0.0;
        double ratio_b = 0.0;
        for (size_t j = 0; j < dimension_; ++j){
          
            ratio_a += solution->lambda(j) * static_cast<double>(profit_.at(i).at(j));
            ratio_b += static_cast<double>(weight_.at(i).at(j));
          
        }
        CL.push_back(make_pair(ratio_a/ratio_b, i));
      }
    }
    
    // ORDENAR
    sort(CL.begin(), CL.end(), decrescente); // ordem decrescente

    int tamRCL=max((int)(alfa*(double)CL.size()), 1);
    //RCL consiste nos alfa% primeiros elementos de CL
    // cout<<"tamRCL = "<<tamRCL<<endl;
    do{
      int index1 = pick_a_number(0, tamRCL-1);
      solution->insert_item(CL[index1].second, profit_.at(CL[index1].second), weight_.at(CL[index1].second));
      if (solution->overload(capacity_)) {
        solution->remove_item(CL[index1].second, profit_.at(CL[index1].second), weight_.at(CL[index1].second));
        break; // sai do laço
      } else { // se vier pro else, entao o item foi inserido na mochila
        CL.erase(CL.begin()+index1);
        tamRCL--;
      }
    } while(tamRCL>0); // enquanto a inserçao do elemento "e" nao viola "x" e enquanto RCL é diferente de vazia

    for (int i=0; i<CL.size(); i++){
      int item = CL[i].second;
      solution->insert_item(item, profit_.at(item), weight_.at(item));
      if (solution->overload(capacity_)) {
        solution->remove_item(item, profit_.at(item), weight_.at(item));
      }
    }

  }












  // @brief
  //   Generate initial population of solutions.
  // @param vector<Knapsack> $population
  //   Population.
  // @param int $population_size
  //   Size of the population.
  void generate_population(std::vector<Knapsack* > &population, const int population_size, const int method, double alfa, int &contAval)
  {

    if (dimension_ == 2) // Two dimensions
    {
      for (int i = 0; i < population_size; ++i)
      {
        
        Knapsack *solution = new Knapsack(size_, dimension_);
        
        // Generate scalarization vector.
        solution->set_lambda(0, static_cast<double>(i) / (population_size - 1));
        solution->set_lambda(1, 1.0 - solution->lambda(0));
        // Generate solution.
        generate_solution(solution,method, contAval);
        // BuildSolution(solution,alfa,contAval);
        population.push_back(solution);
      }
    }
    else if (dimension_ == 3) // Three dimensions
    {
      int range = sqrt(2 * population_size);
      for (int i = 0; i < range; ++i)
      {
        for (int j = 0; j < range - i; ++j)
        {
           Knapsack *solution = new Knapsack(size_, dimension_);
          // Generate scalarization vectors.
          solution->set_lambda(0, static_cast<double>(i) / (range - 1));
          solution->set_lambda(1, static_cast<double>(j) / (range - 1));
          solution->set_lambda(2, 1.0 - solution->lambda(0) - solution->lambda(1));
          // Generate the solution.
          generate_solution(solution,method, contAval);
           // BuildSolution(solution,alfa,contAval);
          population.push_back(solution);
        }
      }
    }
    else if (dimension_ == 4) // Four dimensions
    {
      int range = cbrt(6 * population_size);
      for (int i = 0; i < range; ++i)
      {
        for (int j = 0; j < range - i; ++j)
        {
          for (int k = 0; k < range - j - i; ++k)
          {
             Knapsack *solution = new Knapsack(size_, dimension_);
            // Generate scalarization vectors.
            solution->set_lambda(0, static_cast<double>(i) / (range - 1));
            solution->set_lambda(1, static_cast<double>(j) / (range - 1));
            solution->set_lambda(2, static_cast<double>(k) / (range - 1));
            solution->set_lambda(3, 1.0 - solution->lambda(0) - solution->lambda(1) - solution->lambda(2));
            // Generate the solution.
            generate_solution(solution,method, contAval);
             // BuildSolution(solution,alfa,contAval);
            population.push_back(solution);
          }
        }
      }
    }
    else // More than four dimensions
    {
      for (int i = 0; i < population_size; ++i)
      {
         Knapsack *solution = new Knapsack(size_, dimension_);
        // Generate scalarization vector.
        for (size_t j = 0; j < dimension_; ++j)
        {solution->set_lambda(j, pick_a_number(0.0, 1.0));}
        // Generate solution.
        generate_solution(solution,method, contAval);
       // BuildSolution(solution,alfa,contAval);
        population.push_back(solution);
      }
    }
  } // function generate_population
  // @brief
  //   Generate a new solution using Greed Algorithm.
  // @param Knapsack $solution
  //   Reference to the clean solution.
  // @param int &method
  //   1 = Random, 0 = Greedy
  void generate_solution(Knapsack *solution, const int& method, int &contAval)
  {

    contAval++;
    if (method == 0)
    {
      BinaryHeap<std::pair<double, size_t>> binary_heap;
      double ratio_a;
      double ratio_b;

      for (size_t i = 0; i < size_; ++i)
      {
        ratio_a = 0.0;
        ratio_b = 0.0;
        for (size_t j = 0; j < dimension_; ++j)
        {
          ratio_a -= solution->lambda(j) * static_cast<double>(profit_.at(i).at(j));
          ratio_b += static_cast<double>(weight_.at(i).at(j)) / capacity_.at(j);
        }
        binary_heap.insert(std::make_pair(ratio_a / ratio_b, i));
      }

      while (!binary_heap.empty())
      {
        solution->insert_item(binary_heap.find_min().second, profit_.at(binary_heap.find_min().second), weight_.at(binary_heap.find_min().second));
        if (solution->overload(capacity_))
        {solution->remove_item(binary_heap.find_min().second, profit_.at(binary_heap.find_min().second), weight_.at(binary_heap.find_min().second));}
        binary_heap.delete_min();
      }
    }
    if (method == 1)
    {
      std::vector<size_t> item;
      for (size_t i = 0; i < size_; ++i) {item.push_back(i);}

      do{

        if (item.size()>0){
          int i = pick_a_number(0, item.size()-1); // 0 ... item.size()-1
          solution->insert_item(item[i], profit_.at(item[i]), weight_.at(item[i]));
          if (solution->overload(capacity_)) {solution->remove_item(item[i], profit_.at(item[i]), weight_.at(item[i]));}
        
          item[i] = item[item.size()-1];
          item.pop_back();
        }
      }while(item.size()>0);
    }
  } // function generate_solution

  // @brief
  //   Item profits.
  std::vector<std::vector<int>> profit_;
  // @brief
  //   Item weights.
  std::vector<std::vector<int>> weight_;
  // @brief
  //   Knapsack capacities.
  std::vector<int> capacity_;
  // @brief
  //   Source of instances.
  const std::string type_;
  // @brief
  //   Number of items.
  const size_t size_;
  // @brief
  //   Number of objectives.
  const size_t dimension_;

  // @brief
  // lower_bound[i] is less or equal than the i-th dimenssion of any objective point
  // where i = 0...dimension_-1
  std::vector<int> lower_bound;
  // no caso da mochila, o lower_bound pode ser 0 para todas as dimenssoes, ou seja, a mochila pode estar vazia.
  // mas isso nao vale para todos os problemas (por exemplo, pra AGMO e PQA é diferente)

  // @brief
  // lower_bound[i] is greater or equal than the i-th dimenssion of any objective point
  // where i = 0...dimension_-1
  std::vector<int> upper_bound; 
  // para a mochila, o upper_bound (ponto utopico) consiste numa
  // soluçao que possui todos os itens disponiveis, 
  // sem considerar o limite da mochila


}; // class MOMKP


#endif 

