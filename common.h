#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include <vector>
#include <utility> 
#include <map>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <stack>
#include <bitset>
#include <ctime>
#include <unordered_map>
#include <unordered_set>

using namespace std;

#define SPACESIZE	30

typedef  int* Point; 
typedef  pair<int,int> Order;
typedef int Space;
typedef int id;

void listeAttributsPresents(Space subspace, Space d, vector<Space> &result);

struct pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
  }
};

#endif // DECLARATIONS_H_INCLUDED