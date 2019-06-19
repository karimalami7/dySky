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

struct preference_tree{
	Preference p;
	vector <preference_tree*> preference_child;
	vector<id> ids;
};

class Arg{

public:

	vector<Point> to_dataset;
	vector<vector<int>> po_dataset;
	vector<pair<Preference, vector<id>>> cached_views;
	vector<preference_tree*> sky_view;
	vector<id> skyline_result;

	void compute_views(Config *cfg);
	void compute_view_recursively(Config *cfg, int niveau, vector<Preference> preference_stack, vector<preference_tree*> &view_node);
	void compute_view_1d(Config *cfg, vector<preference_tree*> &view_node);
	void compute_skyline(Config *cfg, Query q);

};

// compute and store a set of views
void Arg::compute_views(Config *cfg){
	
	//cout << "Arg::compute_views"<<endl;

	if(cfg->dyDim_size==1){
		compute_view_1d(cfg, this->sky_view);
	}
	else{
		vector<vector<Preference>> preference_stack(fact(cfg->dyDim_val));
		this->sky_view=vector<preference_tree*>(fact(cfg->dyDim_val));
		#pragma omp parallel for schedule(dynamic)
		for (int number_view_stored=0; number_view_stored< fact(cfg->dyDim_val); number_view_stored++){
			Query q;
			q.generate_preference(cfg);		
			preference_stack[number_view_stored].push_back(q.preference[0]);
			preference_tree *pt =new preference_tree;
			pt->p=q.preference[0];
			this->sky_view[number_view_stored]=pt;
	  		this->compute_view_recursively(cfg, 1, preference_stack[number_view_stored], this->sky_view[number_view_stored]->preference_child);		
		}
	}

}

void Arg::compute_view_1d(Config *cfg, vector<preference_tree*> &view_node){

	#pragma omp parallel for schedule(dynamic)
	for (int number_views_stored=0; number_views_stored< fact(cfg->dyDim_val); number_views_stored++){
		Query q;
		q.generate_preference(cfg);
		Cps cps(cfg);
		cps.decompose_preference(q.preference[0],cfg,0);	
		cps.encoding(cfg);
		cps.to_dataset=this->to_dataset;
		cps.po_dataset=this->po_dataset;
		cps.compute_skyline(cfg);
		preference_tree *pt =new preference_tree;
		pt->p=q.preference[0];
		pt->ids=cps.skyline_result;
		view_node.push_back(pt);
	}
}

void Arg::compute_view_recursively(Config *cfg, int niveau, vector<Preference> preference_stack, vector<preference_tree*> &view_node){
	//cout << "Arg::compute_view_recursively niveau: "<< niveau<<endl; 
	
	for (int number_views_stored=0; number_views_stored< fact(cfg->dyDim_val); number_views_stored++){
		Query q;
		q.generate_preference(cfg);
		// q.preference[0].print_edges();
		// cout << endl;

		if (niveau<cfg->dyDim_size -1){
			preference_stack.push_back(q.preference[0]);
			preference_tree *pt =new preference_tree;
			pt->p=q.preference[0];
			view_node.push_back(pt);
			this->compute_view_recursively(cfg,  niveau+1, preference_stack, view_node[number_views_stored]->preference_child);
		}
		else{
			preference_stack.push_back(q.preference[0]);
			Cps cps(cfg);
			for (int i=0;i<cfg->dyDim_size;i++){
				cps.decompose_preference(preference_stack[i],cfg,i);	
			}
			cps.encoding(cfg);
			cps.to_dataset=this->to_dataset;
			cps.po_dataset=this->po_dataset;
			cps.compute_skyline(cfg);
			preference_tree *pt =new preference_tree;
			pt->p=q.preference[0];
			pt->ids=cps.skyline_result;
			view_node.push_back(pt);
		}
		preference_stack.pop_back();
	}

}

void Arg::compute_skyline(Config *cfg, Query q){

	//cout << "Arg::compute_skyline"<<endl;
	//cout << "number of cached views: "<<this->sky_view.size()<<endl;
	int i;
	vector<bool> refinement_found(cfg->dyDim_size,false);
	vector<preference_tree*> pt=this->sky_view;
	for (int d=0;d<cfg->dyDim_size;d++){
		for (i=0; i<pt.size(); i++){
			if(pt[i]->p.is_subgraph(q.preference[d])){
				refinement_found[d]=true;
				cout <<"refinement found: view"<< i<<endl;
				// avancer au deuxieme niveau
				pt=pt[i]->preference_child;
				break;
			}
		}
		if (refinement_found[d]==false){
			break;
		}	
	}

	//for (int d=0;d<cfg->dyDim_size;d++) cout << refinement_found[d]<<endl;
	
	Cps cps(cfg);
	for (int d=0;d<cfg->dyDim_size;d++){
		cps.decompose_preference(q.preference[d],cfg,d);	
	}

	if ( adjacent_find( refinement_found.begin(), refinement_found.end(), not_equal_to<bool>() ) == refinement_found.end() && refinement_found[0]==true){
		cerr << "Refinement found"<<endl; 
		for (int j=0; j<pt[i]->ids.size();j++){
			cps.to_dataset.push_back(this->to_dataset[pt[i]->ids[j]]);
			cps.po_dataset.push_back(this->po_dataset[pt[i]->ids[j]]);
		}
		cerr <<"dataset size: "<<cps.to_dataset.size()<<endl; 
	}
	else{
		cerr << "No view found"<<endl;
		cps.to_dataset=this->to_dataset;
		cps.po_dataset=this->po_dataset; 
	}
	

	cps.encoding(cfg);
	cps.compute_skyline(cfg);
	this->skyline_result=cps.skyline_result;
	
}