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

	void compute_views(Config *cfg, uint64_t *storage);
	void compute_view_recursively(Config *cfg, int niveau, vector<Preference> preference_stack, vector<preference_tree*> &view_node, uint64_t *storage);
	void compute_view_1d(Config *cfg, vector<preference_tree*> &view_node, uint64_t *storage);
	void compute_skyline(Config *cfg, Query q);

};

// compute and store a set of views
void Arg::compute_views(Config *cfg, uint64_t *storage){
	
	cout << "Arg::compute_views"<<endl;

	if(cfg->dyDim_size==1){
		compute_view_1d(cfg, this->sky_view, storage);
	}
	else{
		vector<vector<Preference>> preference_stack(fact(cfg->dyDim_val));
		this->sky_view=vector<preference_tree*>(fact(cfg->dyDim_val));
		#pragma omp parallel for schedule(dynamic)
		for (int id_view=0; id_view< fact(cfg->dyDim_val); id_view++){
			cout << "id_view: "<< id_view<<endl;
			Query q;
			q.generate_preference(cfg);		
			preference_stack[id_view].push_back(q.preference[0]);
			preference_tree *pt =new preference_tree;
			pt->p=q.preference[0];
			this->sky_view[id_view]=pt;
	  		this->compute_view_recursively(cfg, 1, preference_stack[id_view], this->sky_view[id_view]->preference_child, storage);		
		}
	}

}

void Arg::compute_view_1d(Config *cfg, vector<preference_tree*> &view_node, uint64_t *storage){

	#pragma omp parallel for schedule(dynamic)
	for (int number_views_stored=0; number_views_stored< fact(cfg->dyDim_val); number_views_stored++){
		Query q;
		q.generate_preference(cfg);
		Cps *cps_arg=new Cps(cfg);
		cps_arg->decompose_preference(q.preference[0],cfg,0);	
		cps_arg->encoding(cfg);
		cps_arg->compute_skyline(cfg, false);
		preference_tree *pt =new preference_tree;
		pt->p=q.preference[0];
		*storage=*storage+cps_arg->skyline_result.size();
		//pt->ids.swap(cps_arg->skyline_result);
		view_node.push_back(pt);
		delete cps_arg;
	}
}

void Arg::compute_view_recursively(Config *cfg, int niveau, vector<Preference> preference_stack, vector<preference_tree*> &view_node, uint64_t *storage){
	// cout << "Arg::compute_view_recursively niveau: "<< niveau<<endl; 
	
	//************************************
	// memory
	// struct sysinfo sys_info;
	// uint64_t info_ram;
	// if (!(sysinfo(&sys_info) == -1)) {
	// 	info_ram=sys_info.totalram - sys_info.freeram;
	// 	info_ram = (info_ram * sys_info.mem_unit)/1024;
	// 	printf("memory %d\n", info_ram);
	// }
	
	for (int id_view=0; id_view< fact(cfg->dyDim_val); id_view++){
		// cout << "id_view: "<< id_view<<endl;
		Query q;
		q.generate_preference(cfg);
		// q.preference[0].print_edges();
		// cout << endl;

		if (niveau<cfg->dyDim_size -1){
			preference_stack.push_back(q.preference[0]);
			preference_tree *pt =new preference_tree;
			pt->p=q.preference[0];
			view_node.push_back(pt);
			this->compute_view_recursively(cfg,  niveau+1, preference_stack, view_node[id_view]->preference_child, storage);
		}
		else{
			preference_stack.push_back(q.preference[0]);
			Cps *cps_arg=new Cps(cfg);
			for (int i=0;i<cfg->dyDim_size;i++){
				cps_arg->decompose_preference(preference_stack[i],cfg,i);	
			}
			cps_arg->encoding(cfg);
			cps_arg->compute_skyline(cfg, false);
			preference_tree *pt =new preference_tree;
			pt->p=q.preference[0];
			*storage=*storage+cps_arg->skyline_result.size();
			//pt->ids.swap(cps_arg->skyline_result);
			view_node.push_back(pt);
			delete cps_arg;
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
	
	Cps cps_arg(cfg);
	for (int d=0;d<cfg->dyDim_size;d++){
		cps_arg.decompose_preference(q.preference[d],cfg,d);	
	}

	if ( adjacent_find( refinement_found.begin(), refinement_found.end(), not_equal_to<bool>() ) == refinement_found.end() && refinement_found[0]==true){
		cerr << "Refinement found"<<endl; 
		for (int j=0; j<pt[i]->ids.size();j++){
			cps_arg.to_dataset.push_back(cfg->to_dataset[pt[i]->ids[j]]);
			cps_arg.po_dataset.push_back(cfg->po_dataset[pt[i]->ids[j]]);
		}
		cerr <<"dataset size: "<<cps_arg.to_dataset.size()<<endl; 
		cps_arg.encoding(cfg);
		cps_arg.compute_skyline(cfg, true);
		this->skyline_result=cps_arg.skyline_result;
	}
	else{
		cerr << "No view found"<<endl;
		cps_arg.encoding(cfg);
		cps_arg.compute_skyline(cfg, false);
		this->skyline_result=cps_arg.skyline_result;
	}
	
}