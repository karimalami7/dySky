/*
 * tos.h
 *
 *  Created on: Feb 2, 2019
 *      Author: karim
 *
 * Functionalities:  
 * 		+ 
 */

#include "common.h"
#include "config.h"
#include "generator/generateur.h"
#include "BSkyTree/bskytree.h"
using namespace std;

class Tos{

public:

	vector<Point> to_dataset;
	vector<int> po_dataset;
	vector<Graph<int>> cached_preferences;
	vector<vector<id>> cached_skylines;

};