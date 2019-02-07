/*
 * arg.h
 *
 *  Created on: January ??, 2019
 *      Author: karim
 *
 * Functionalities:  
 * 		+ Decompose a preference into chains
 * 		+ Modify the dataset
 *		+ Compute skyline by BSkyTree
 */

#include "common.h"
#include "config.h"
#include "generator/generateur.h"
#include "BSkyTree/bskytree.h"
using namespace std;

class Arg{

public:

	vector<Point> to_dataset;
	vector<int> po_dataset;
	vector<Graph<int>> cached_preferences;
	vector<vector<id>> cached_skylines;

};