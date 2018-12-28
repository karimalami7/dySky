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
typedef  std::pair<int,int> Order;
typedef int  Space;

void listeAttributsPresents(Space subspace, Space d, vector<Space> &result);