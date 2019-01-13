/*
 * dySky.h
 *
 *  Created on: January 13, 2019
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

class Cps{
	vector<Graph<int>> chains;
public: 
	void decompose_preference(Preference p);
	Graph<int> compute_transitive_closure(Preference p);
};

void Cps::decompose_preference(Preference p){

Graph<int> transitive_preference=Cps::compute_transitive_closure(p);

}

Graph<int> Cps::compute_transitive_closure(Preference p){

	Graph<int> g;

}