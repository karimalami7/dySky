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
#include <list>
#include "config.h"
#include <sys/sysinfo.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

using namespace std;

#define SPACESIZE	30

typedef  int* Point; 
typedef  pair<int,int> Order;
typedef int Space;
typedef int id;
typedef vector<int> chain;

void listeAttributsPresents(Space subspace, Space d, vector<Space> &result);

long fact(const int n);

struct pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
  }
};

struct order_tree{
	unordered_map<Order, order_tree*, pairhash> order_child;
	vector<id> ids;
};

struct chain_tree{
	map<chain, chain_tree*> chain_child;
	vector<id> ids;
};


void generate_all_orders(Config *cfg, vector<vector<Order>> &all_orders);

void printConfig(Config *cfg);

#endif // DECLARATIONS_H_INCLUDED