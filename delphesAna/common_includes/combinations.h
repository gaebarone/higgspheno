#ifndef COMBINATIONS_H
#define COMBINATIONS_H

//  (C) Copyright Howard Hinnant 2005-2011.
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).
//
//  See http://www.boost.org/libs/type_traits for most recent version including documentation.

//  Details are in namespace detail.  Every effort has been made to make
//  combine_discontinuous and permute as fast as possible.  They minimize the number
//  of swaps that are performed. Everything else builds on these two primitives.
//  The most complicated algorithm is for_each_reversible_permutation.  But it
//  builds on combine_discontinuous and permute and I believe represents a minimum
//  number of swaps.  Without care, algorithms such as for_each_reversible_permutation
//  will take longer than for_each_permutation instead of the intended half the time.

//  Speed is everything.  Lest you could just use std::next_permutation and manually
//  eliminate duplicate permutations.  If the implementation fails in being orders
//  of magnitude faster than that, then it has failed miserably.

#include <iterator>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include "other_combinations.h"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <iomanip>



using namespace std;

std::vector<std::vector<int>> combinationsNoRepetitionAndOrderDoesNotMatter (int subsetSize, std::vector<int> setOfNumbers){
  std::vector<std::vector<int> > subsets{};
  subsets.reserve (count_each_combination (setOfNumbers.begin (), setOfNumbers.begin () + subsetSize, setOfNumbers.end ()));
  for_each_combination (setOfNumbers.begin (), setOfNumbers.begin () + subsetSize, setOfNumbers.end (), [&subsets] (auto first, auto last) {
      subsets.push_back (std::vector<int>{ first, last });
      return false;
    });
  return subsets;
}

std::vector <std::vector <int>> getAllCombinations(vector<int> inputVector, int k){
  std::vector<vector<int>> combinations; 
  std::vector<int> selector(inputVector.size());
  std::fill(selector.begin(), selector.begin() + k, 1);

  do {
    std::vector<int> selectedIds;
    std::vector<int> selectedVectorElements;
    for (int i = 0; i < inputVector.size(); i++) {
      if (selector[i]) {
	selectedIds.push_back(i);
      }
    }
    for (auto& id : selectedIds) {
      selectedVectorElements.push_back(inputVector[id]);
    }
    combinations.push_back(selectedVectorElements);
  } while (std::prev_permutation(selector.begin(), selector.end()));

  return combinations;
}


#endif
