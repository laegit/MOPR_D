// @file Knapsack.cpp
// @author Denis Felipe
// @version 0.1
//
// @brief
//   Multiobjective Multidimensional Knapsack problem solver.
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
#ifndef KNAPSACK_CPP
#define KNAPSACK_CPP

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
#include "param_geral.h"
#include "momkp.cpp"


using namespace std;




// @struct Knapsack
// @brief
//   Solution for the Multiobjective Multidimensional Knapsack problem.
class Knapsack
{
 public:
  // Default constructor
  Knapsack() {;} // ctor
  // Constructor
  // @param size_t $size
  //   Number of items.
  // @param size_t $dimension
  //   Number of objectives.
  Knapsack(const size_t& size, const size_t& dimension)
  : knapsack_(size, false), profit_(dimension), weight_(dimension), lambda_(dimension), distance(0.0),profit_NORMALIZADO(dimension)
  {;} // ctor
  // Copy constructor
  explicit Knapsack(const Knapsack& rhs)
  : knapsack_(rhs.knapsack_), profit_(rhs.profit_), weight_(rhs.weight_), lambda_(rhs.lambda_), distance(rhs.distance), profit_NORMALIZADO(rhs.profit_NORMALIZADO)
  {;} // ctor
  // Move constructor
  explicit Knapsack(const Knapsack&& rhs)
    : knapsack_(std::move(rhs.knapsack_)), profit_(std::move(rhs.profit_)), weight_(std::move(rhs.weight_)), lambda_(std::move(rhs.lambda_)), distance(0.0), profit_NORMALIZADO(std::move(rhs.profit_NORMALIZADO))
  {;} // ctor
  // Destructor
  ~Knapsack()
  {;} // dtor
  // Assignment operator
  const void operator=(const Knapsack& rhs)
  {
    
    knapsack_ = rhs.knapsack_; // TODO: testar se, dado s1<-s2, editando s2 nao edita s1
    profit_ = rhs.profit_;
    weight_ = rhs.weight_;
    lambda_ = rhs.lambda_;
    rank = rhs.rank;
    distance = rhs.distance;
    posicaoListaNSGAII = rhs.posicaoListaNSGAII;
    profit_NORMALIZADO = rhs.profit_NORMALIZADO;
    
    // return *this; // Allow cascaded assignment.
  }

  const bool operator==(Knapsack &s) {

    for (size_t i=0; i<profit_.size(); i++){
      if (profit_[i]!=s.profit(i)) return false;
    }
    return true;
  }

  // esta rotina funciona melhor quando o profit_ é inteiro (o que é nosso caso)
  bool operator>> (Knapsack &d) { // retorna true se profit_ domina d.profit()
    bool greater = false;
    bool lesser = false;
    // Check if profit_ is lesser or greater than d.profit in any dimension.
    for (size_t i = 0; i < profit_.size(); ++i)
    {
      if (profit_.at(i) > d.profit(i))
      {greater = true;}
      else if (profit_.at(i)  < d.profit(i))
      {lesser = true;}
    }

    if (greater)
    {
      // if (lesser) // Non-dominated points
      // {return 0;}
      if (lesser==false) // $point_a dominates $point_b
        return true;
      else return false;
    } else return false;
  }



  // // a epsilin_dominance considera que os objetivos (profits) estao normalizados entre 1 e 2
  // // this operator implements the (aditive) épislon-dominance as defined as in expression (7) in Laumanns et al (2002) 
  // // Then profit_ is said to epislon-dominate d.profit() for some  epislon > 0, iff
  // // for all i ∈ {1,...,dimenssion} epsilon + profit_.at(i) >= d.profit(i) 
  // bool epsilon_dominance (Knapsack &d, double epsilon, std::vector<int> lb, std::vector<int> ub) {
  //   bool greater = false;
  //   bool lesser = false;
  //   // Check if profit_ is lesser or greater than d.profit in any dimension.
  //   for (size_t i = 0; i < profit_.size(); ++i)
  //   {
  //     double profit_at_i = 1.0+((double)profit_.at(i)-lb[i])/(ub[i]-lb[i]);
  //     double d_profit_i =  1.0+((double)  d.profit(i)-lb[i])/(ub[i]-lb[i]);
  //     // se o problema fosse de minimizaçao, esse epsilon+profit_.at(i) seria epsilon-profit_.at(i)
  //     if (epsilon+profit_at_i >= d_profit_i) //TODO: colocar tolerância de 1e-9 
  //     // if (epsilon+profit_at_i - d_profit_i >= 1e-9)
  //     {greater = true;}
  //     else 
  //     {lesser = true;}
  //   }

  //   if (greater)
  //   {
  //     // if (lesser) // Non-dominated points
  //     // {return 0;}
  //     if (lesser==false) // $point_a epislon_dominates $point_b
  //       return true;
  //     else return false;
  //   } else return false;
  // }

  // @brief
  //   Add item to knapsack.
  // @param size_t $a
  //   knapsack.
  // @param vector<int>
  //   Profit.
  // @param vector<int>
  //   Weight.
  void insert_item(const size_t& a, const std::vector<int>& profit, const std::vector<int>& weight)
  { 
    knapsack_.at(a) = true;
    for (size_t i = 0; i < weight_.size(); ++i)
    {
      profit_.at(i) += profit.at(i);
      weight_.at(i) += weight.at(i);
    }
  } // function insert_item
  // @brief
  //   Remove item from knapsack.
  // @param size_t $a
  //   knapsack.
  // @param vector<int>
  //   Profit.
  // @param vector<int>
  //   Weight.
  void remove_item(const size_t& a, const std::vector<int>& profit, const std::vector<int>& weight)
  {
    knapsack_.at(a) = false;
    for (size_t i = 0; i < weight_.size(); ++i)
    {
      profit_.at(i) -= profit.at(i);
      weight_.at(i) -= weight.at(i);
    }
  } // function remove_item


  // Get item
  // @param size_t $a
  //   Item.
  // @return bool
  //   True if knapsack contains item.
  const bool at(const size_t& a) const
  {return knapsack_.at(a);} // function at


  // Get profit
  // @param size_t $objective
  //   Objective index.
  // @return int
  //   Profit of objective.
  const int profit(const size_t& objective) const
  {return profit_.at(objective);} // function profit
  // Get profit
  // @return vector<int>
  //   Profit vector.
  const std::vector<int>& profit() const
  {return profit_;} // function profit
  // Get weight
  // @param size_t $objective
  //   Objective index.
  // @return int
  //   Weight of objective.
  const int weight(const size_t& objective) const
  {return weight_.at(objective);} // function weight
  // Get weight
  // @return vector<int>
  //   Weight vector.
  const std::vector<int>& weight() const
  {return weight_;} // function weight
  // Get lambda
  // @param size_t $objective
  //   Objective index.
  // @return double
  //   Factor of objective.
  const double lambda(const size_t& objective) const
  {return lambda_.at(objective);} // function lambda
  // Set lambda
  // @param size_t $objective
  //   Objective index.
  // @param double $value
  //   Factor of objective.
  void set_lambda(const size_t& objective, const double& value)
  {lambda_.at(objective) = value;} // function set_lambda
  // @brief
  //   Use weighted sum method to decompose objectives into single profit.
  // @return double
  //   Decomposed weight.
  const double decompose() const
  {
    double total_profit = 0.0;
    for (size_t i = 0; i < profit_.size(); ++i)
    {total_profit += static_cast<double>(profit(i)) * lambda(i);}
    return total_profit;
  } // function decompose
  const bool overload(const std::vector<int>& capacity)
  {
    for (size_t i = 0; i < weight_.size(); ++i)
    {
      if (weight(i) > capacity.at(i))
      {return true;}
    }
    return false;
  }

  size_t getSize(){
    return knapsack_.size();
  }

  void print_solution(FILE *f){

    // fprintf(f,"Ponto objetivo = ");
    for (int i=0; i<profit_.size(); i++){
      fprintf(f,"%i ",profit_[i]);
    }
    // fprintf(f,"\nItens na mochila = ");
    // for (int i=0; i<knapsack_.size(); i++){
    //   if (knapsack_[i]==true){
    //     fprintf(f,"%i ",i);
    //   }
    // }
    // fprintf(f,"\nPesos dimensionais da mochila = ");
    // for (int i=0; i<weight_.size(); i++){
    //   fprintf(f,"%i ",weight_[i]);
    // }
    fprintf(f, "\n");
  }



 private:
  
  // @brief
  //   Item list.
  std::vector<bool> knapsack_;
  // @brief
  //   Profit for each objective.
  std::vector<int> profit_;
  // @brief
  //   Weight for each objective.
  std::vector<int> weight_;
  // @brief
  //   Scalarization vector.
  std::vector<double> lambda_;

  public:
   // @brief
  // crownd distance used by NSGAA-II 
  double distance; 
  int rank;
  int posicaoListaNSGAII; // guarda o index onde a soluçao é guardada na popupacao NUMPOPULACAO*2 do NSGA-II

  std::vector<double> profit_NORMALIZADO;


}; // class Knapsack

#endif 
