// @file random.h
// @author Denis Felipe
// @version 0.1
// 
// @brief
//   Pseudo-random number generator.
// @warning
//   Contains C++11 features.
// 
// @example
//   #include <iostream>
//   #include "random.h"
//   int main()
//   {
//     rng::randomize(); // Select a random seed for the pseudo-random number generator.
//     std::cout << rng::pick_a_number(1, 6) << std::endl; // Roll a die.
//
//     return 0;
//   } // Funcion main
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
#ifndef RANDOM_CPP
#define RANDOM_CPP

// Random devices and distributions
#include <random>


// @brief
//   Access the random number engine used to generate pseudo-random numbers.
//
// @return default_random_engine
//   Random number engine.
std::mt19937& global_urng()
{
  static std::mt19937 u{};
  return u;
} // function global_urng

// @brief
//   Use a true random number generator to re-initialize the internal state value of the random number engine.
void randomize()
{
  static std::random_device rd{};
  global_urng().seed(rd());
} // function randomize

// @brief
//   Set a seed to re-initialize the internal state value of the random number engine.
// 
// @param unsigned int $seed
//   Seeding value.
void randomize(const unsigned int& seed)
{global_urng().seed(seed);} // function randomize

// @brief
//   Generate a pseudo-random integer number in the range [$from, $thru] according to a uniform discrete distribution.
//
// @param int $from
//   Lower bound of range.
// @param int $thru
//   Upper bound of range.
const int pick_a_number(const int& from, const int& thru)
{
  static std::uniform_int_distribution<> d{};
  using parm_t = decltype(d)::param_type;
  return d(global_urng(), parm_t{from, thru});
} // function pick_a_number

// @brief
//   Generate a pseudo-random real number in the range [$from, $upto) according to a uniform distribution.
// 
// @param int $from
//   Lower bound of range.
// @param int $upto
//   Upper bound of range.
const double pick_a_number(const double& from, const double& upto)
{
  static std::uniform_real_distribution<> d{};
  using parm_t = decltype(d)::param_type;
  return d(global_urng(), parm_t{from, upto});
} // function pick_a_number


#endif 

