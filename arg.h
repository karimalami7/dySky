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
	vector<pair<Preference, vector<id>>> cached_views;
	vector<id> skyline_result;

	void compute_views(Config *cfg);

	void compute_skyline(Config *cfg, Query q);

};

// compute and store a set of views
void Arg::compute_views(Config *cfg){
	
	cout << "Arg::compute_views"<<endl;

	int number_views_stored=0;
	while (number_views_stored< fact(cfg->dyDim_val)){
		Query q;
		q.generate_preference(cfg);
		Cps cps;
		cps.decompose_preference(q.preference,cfg);
		cps.to_dataset=this->to_dataset;
		cps.po_dataset=this->po_dataset;
		cps.compute_skyline(cfg);
		this->cached_views.push_back(pair<Preference, vector<id>>(q.preference, cps.skyline_result));
		number_views_stored++;
	}
}

void Arg::compute_skyline(Config *cfg, Query q){

	cout << "Arg::compute_skyline"<<endl;
	cout << "number of cached views: "<<cached_views.size()<<endl;
	int i;
	bool refinement_found=false;
	for (i=0; i<cached_views.size(); i++){
		if(cached_views[i].first.is_subgraph(q.preference)){
			refinement_found=true;
			cout <<"refinement found: view"<< i<<endl;
			break;
		}
	}
	Cps cps;
	cps.decompose_preference(q.preference,cfg);
	if(refinement_found){
		for (int j=0; j<cached_views[i].second.size();j++){
			cps.to_dataset.push_back(this->to_dataset[cached_views[i].second[j]]);
			cps.po_dataset.push_back(this->po_dataset[cached_views[i].second[j]]);
		}
		cerr << "Refinement, dataset size: "<<cps.to_dataset.size()<<endl; 
	}
	else{
		cps.to_dataset=this->to_dataset;
		cps.po_dataset=this->po_dataset;
		cerr << "No view found, dataset size: "<<cps.to_dataset.size()<<endl; 
	}

	cps.compute_skyline(cfg);
	this->skyline_result=cps.skyline_result;
}