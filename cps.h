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
	void decompose_preference(Graph<int> p);
	
};

void Cps::decompose_preference(Graph<int> p){

	Graph<int> transitive_preference;

	// compute transitive

	transitive_preference.compute_transitive_closure(p);

	// find missing order

	vector<Order> missing_orders;
	cout << "missing orders"<<endl;
	for (auto it = transitive_preference.vertices.begin(); it != transitive_preference.vertices.end(); it++){
		if (transitive_preference.out_edges.find(*it)!=transitive_preference.out_edges.end()){
			for (int i=0;i<transitive_preference.vertices.size();i++){
				if (transitive_preference.out_edges[*it].find(i)==transitive_preference.out_edges[*it].end() && (*it)!=i){
					missing_orders.push_back(pair<id,id>((*it),i));
					cout <<(*it)<<i<<endl;
				}
			}
		}
		else{
			for (int i=0;i<transitive_preference.vertices.size();i++){
				if ((*it)!=i){
					missing_orders.push_back(pair<id,id>((*it),i));
					cout <<(*it)<<i<<endl;
				}
			}
		}
	}

	// remove induced orders

}


