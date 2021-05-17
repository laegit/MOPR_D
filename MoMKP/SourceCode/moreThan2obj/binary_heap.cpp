// @file binary_heap.cpp
// @author Denis Felipe
// @version 0.1
//
// @brief
//   Binary heap data structure.
//
// @example
//   #include <iostream>
//   #include "binary_heap.cpp"
//
//   int main(int argc, char** argv)
//   {
//    data_structure::BinaryHeap<int> binary_heap;
//
//    binary_heap.insert(2);
//    binary_heap.insert(0);
//    binary_heap.insert(4);
//    binary_heap.insert(1);
//    binary_heap.decrease_key(0, 3);
//
//    while (!binary_heap.empty())
//    {
//      std::cout << binary_heap.find_min() << std::endl;
//      binary_heap.delete_min();
//    }
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
#ifndef BINARY_HEAP_CPP
#define BINARY_HEAP_CPP

#include <algorithm>
#include <vector>


// @class BinaryHeap
// @brief
//   Binary Heap data structure.
template <class T>
class BinaryHeap
{
 public:
  // Constructor
  BinaryHeap()
  {;} // ctor
  // Constructor
  // @param size_t $capacity
  //   Initial capacity of the heap.
  BinaryHeap(const size_t& capacity)
  {binary_heap_.reserve(capacity);} // ctor
  // Destructor
  ~BinaryHeap()
  {;} // dtor

  // @brief
  //   Check if heap is empty.
  // @return bool
  //   True if heap has zero elements, false otherwise.
  const bool empty() const
  {return binary_heap_.empty();} // function empty
  // @brief
  //   Access number of elements in heap.
  // @return size_t
  //   Number of elements.
  const size_t size() const
  {return binary_heap_.size();} // function size
  // @brief
  //   Access first element in heap.
  // @return T
  //   First element.
  const T& find_min() const
  {return binary_heap_.front();} // function find_min
  // @brief
  //   Add a new element in heap, effectively increasing the container size by one.
  // @param T $element
  //   New element.
  void insert(const T& element)
  {
    binary_heap_.push_back(element);
    std::push_heap(binary_heap_.begin(), binary_heap_.end(), [](const T& a, const T& b){return a > b;});
  } // function insert
  // @brief
  //   Remove the minimum element of heap, effectively decreasing the container size by one.
  void delete_min()
  {
    std::pop_heap(binary_heap_.begin(), binary_heap_.end(), [](const T& a, const T& b){return a > b;});
    binary_heap_.pop_back();
  } // function delete_min
  // @brief
  //   Random access to element in heap.
  // @param size_t $index
  //   Index of element.
  // @return T
  //   Element at index $index.
  const T& at(const size_t& index)
  {return binary_heap_.at(index);} // function at
  // @brief
  //   Update an element in heap to a lower key.
  // @param size_t $index
  //   Index of element.
  // @param T $element
  //   New value with a lower key.
  void increase_key(const size_t& index, const T& element)
  {
    binary_heap_.at(index) = element;

    size_t current = index;
    size_t parent = current / 2;
    while (current > 0 && binary_heap_.at(current) < binary_heap_.at(parent))
    {
      std::swap(binary_heap_.at(current), binary_heap_.at(parent));
      current = parent;
      parent = current / 2;
    }
  } // function increase_key
  // @brief
  //   Update an element in heap to a lower key.
  // @param size_t $index
  //   Index of element.
  // @param T $element
  //   New value with a lower key.
  void decrease_key(const size_t& index, const T& element)
  {
    binary_heap_.at(index) = element;

    size_t current;
    size_t left;
    size_t right;
    size_t lesser = index;

    do
    {
      current = lesser;
      left = current * 2;
      right = (current * 2) + 1;

      if (left < binary_heap_.size() && binary_heap_.at(left) < binary_heap_.at(lesser))
      {lesser = left;}
      if (right < binary_heap_.size() && binary_heap_.at(right) < binary_heap_.at(lesser))
      {lesser = right;}
      if (lesser != current)
      {std::swap(binary_heap_.at(current), binary_heap_.at(lesser));}
    }
    while (lesser != current);
  } // function decrease_key
 private:
  // @brief
  //   Heap container.
  std::vector<T> binary_heap_;
}; // class BinaryHeap


#endif

