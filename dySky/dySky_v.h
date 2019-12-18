/*
 * dySky.h
 *
 *  Created on: December 19, 2018
 *      Author: karim
 *
 * Functionalities:  
 * 		+ Decide and compute relevant views, i.e. skyline wrt a a giver preference R
 * 		+ Given a query, decide the views to compute candidate skyline (minimal views, minimal number of candidate)
 *		+ Compute skyline from candidate skyline if necessary.
 */

#include "../common/common.h"
#include "../common/config.h"
#include "../generator/generateur.h"
#include "../BSkyTree/bskytree.h"
#include <omp.h>
using namespace std;

class dySky_v {
	public:
	vector<Point> to_dataset;
	vector<vector<int>> po_dataset;
	vector<map<int, vector<id>>> partition_ids;
	vector<id> never_sky;
	vector<id> candidates;

	dySky_v(Config *cfg);

	int compute_candidates(Config* cfg);
	void compute_views(Config* cfg, vector<vector<Order>> preference_orders, bool* notSkyline);
	void compute_view_1d(Config* cfg, vector<Point> &dataset, vector<vector<Order>> preference_orders, bool* notSkyline);
	void compute_view_recursively_md(Config* cfg, int niveau, vector<Point> &dataset, vector<vector<Order>> preference_orders, bool* notSkyline);
	vector<id> compute_skyline(Config* cfg, vector<vector<Order>> preference_orders);
	void find_dominated(Config* cfg, vector<Point> &dominating, vector<Point> &dominated, vector<id> &notSky);

};

dySky_v::dySky_v(Config *cfg){
	this->po_dataset=vector<vector<int>>(cfg->dataset_size);
}


int dySky_v::compute_candidates(Config* cfg){
	cout << "dySky::compute_candidates"<<endl;

	// cluster dataset depending on the value in the po dimenion
	map<vector<int>, vector<Point>> partitions;
	for (int i=0; i<cfg->dataset_size; i++){
		partitions[this->po_dataset[i]].push_back(this->to_dataset[i]);
	}
	vector<vector<Point>> partitions_vec(partitions.size());
	int i=0;
	for (auto it_partition=partitions.begin(); it_partition!=partitions.end(); it_partition++){
		partitions_vec[i]=it_partition->second;
		i++;
	}

	// compute skyline for every partition -> always skyline + candidates
	int All = (1<<cfg->statDim_size)-1;
	vector<Space> full_Space;
	listeAttributsPresents(All, cfg->statDim_size, full_Space);

	bool* candidates_bool=new bool[cfg->dataset_size];
	#pragma omp parallel for
    for (int i=0; i<partitions_vec.size(); ++i){
    	vector<id> Skyline;
    	Skyline=subspaceSkylineSize_TREE(full_Space, partitions_vec[i]);
    	for (auto elm: Skyline) candidates_bool[elm]=true;
    }
	
	for (int i=0; i<cfg->dataset_size; ++i){
		if(candidates_bool[i]==true){
			this->candidates.push_back(i);
		}else {
			this->never_sky.push_back(i);
		}
	}
    
    //sort (this->candidates.begin(),this->candidates.end());   
    cout << "Candidate set size: "<<this->candidates.size()<<endl;
  	cout << "never_sky set size: "<<this->never_sky.size()<<endl;
}


void dySky_v::compute_view_1d(Config* cfg, vector<Point> &dataset, vector<vector<Order>> preference_orders, bool* notSkyline){
	
	#pragma omp parallel for schedule(dynamic)
	for(int i=0;i<preference_orders[0].size();i++){
		int best_value=preference_orders[0][i].first;
		int worst_value=preference_orders[0][i].second;
		//cout << best_value << ":" << worst_value <<endl;

		
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// partitionner les donnees par rapport aux valeurs best_value, worst_value
		vector<Point> branch_dataset;
		for (int i=0; i<dataset.size(); i++){
			if (dataset[i][cfg->statDim_size+1]==best_value && !notSkyline[dataset[i][0]]){
				Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
				memcpy(p, dataset[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
				p[cfg->statDim_size+1]=0;
				branch_dataset.push_back(p);
			}else if (dataset[i][cfg->statDim_size+1]==worst_value && !notSkyline[dataset[i][0]]){
				Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
				memcpy(p, dataset[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
				p[cfg->statDim_size+1]=1;
				branch_dataset.push_back(p);
			}
		}

		//*************************************
		// calculer le skyline
		int All = (1<<(cfg->statDim_size+cfg->dyDim_size))-1;
		vector<Space> full_Space;
		listeAttributsPresents(All, cfg->statDim_size+cfg->dyDim_size, full_Space);
		//vector<id> sky=subspaceSkylineSize_TREE(full_Space, branch_dataset);
		ExecuteBSkyTree_bis(full_Space, branch_dataset, notSkyline);

	  	//******************************************
	  	//sort ids
	 //  	std::sort (sky.begin(),sky.end());

		// //******************************************
		// // compute dominated points
		// vector<id> branch_dataset_ids(branch_dataset.size());
		// for (int b=0;b<branch_dataset.size();b++){
		// 	branch_dataset_ids[b]=branch_dataset[b][0];
		// }
		// //sort(branch_dataset_ids.begin(),branch_dataset_ids.end());
		// std::vector<id> notSky(branch_dataset_ids.size());   
	 //  	std::vector<id>::iterator it4;
	 //  	it4=std::set_difference(branch_dataset_ids.begin(), branch_dataset_ids.end(), sky.begin(), sky.end(), notSky.begin());                        
	 //  	notSky.resize(it4-notSky.begin());
	 //  	for (auto id: notSky) {
	 //  		notSkyline[id]=true;
	 //  	}
  		//////
		//destroy Point pointers
		for ( auto p : branch_dataset)
		delete p;
		//	
	}
}



void dySky_v::compute_views(Config* cfg, vector<vector<Order>> preference_orders, bool* notSkyline){
	cout << "dySky::compute_views" <<endl;

	vector<Point> candidates_tuples;
	for (int i=0; i<this->candidates.size(); i++){
		Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
		for (int j=0;j<=cfg->statDim_size;j++){
			p[j]=to_dataset[this->candidates[i]][j];
		}
		for (int j=0;j<cfg->dyDim_size;j++){
			p[j+(cfg->statDim_size+1)]=po_dataset[this->candidates[i]][j];
		}
		candidates_tuples.push_back(p);
	}

	if (cfg->dyDim_size==1){
		compute_view_1d(cfg, candidates_tuples, preference_orders, notSkyline);
	}
	else{ // case of multiple dimensions

		int niveau=0;
		#pragma omp parallel for schedule(dynamic)
		for(int s=0;s<preference_orders[0].size();s++){
			int best_value=preference_orders[0][s].first;
			int worst_value=preference_orders[0][s].second;
			//cout <<"niveau: "<<niveau<< " --> "<< best_value<<" : "<<worst_value<<endl;		
			
			vector<Point> branch_dataset;

			if (best_value!=worst_value){
				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				// partitionner les donnees par rapport à cette dimension
				for (int i=0; i<candidates_tuples.size(); i++){
					if (candidates_tuples[i][cfg->statDim_size+1]==best_value && !notSkyline[candidates_tuples[i][0]]){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, candidates_tuples[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						p[1+cfg->statDim_size+niveau]=0;
						branch_dataset.push_back(p);
					}else if (candidates_tuples[i][cfg->statDim_size+1]==worst_value && !notSkyline[candidates_tuples[i][0]]){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, candidates_tuples[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						p[1+cfg->statDim_size+niveau]=1;
						branch_dataset.push_back(p);
					}
				}
			}	
			else{
				for (auto tuple : candidates_tuples){
					if(tuple[1+cfg->statDim_size+niveau]==best_value && !notSkyline[tuple[0]]){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, tuple, (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						branch_dataset.push_back(p);
					}
				}
			}
			compute_view_recursively_md(cfg, niveau+1, branch_dataset, preference_orders, notSkyline);
			
			/////
			//destroy Point pointers
			for ( auto p : branch_dataset)
			delete p;				
			
		}
	}
}

void dySky_v::compute_view_recursively_md(Config* cfg, int niveau, vector<Point> &dataset, vector<vector<Order>> preference_orders, bool* notSkyline){
	
	//#pragma omp parallel for schedule(dynamic) if (niveau==1) 
	for(int s=0;s<preference_orders[niveau].size();s++){
		int best_value=preference_orders[niveau][s].first;
		int worst_value=preference_orders[niveau][s].second;
		// cout <<"niveau: "<<niveau<< " --> "<< best_value<<" : "<<worst_value<<endl;		


		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//if we are not processing last dimension			
		if (niveau<cfg->dyDim_size -1){
			vector<Point> branch_dataset;
			
			if (best_value!=worst_value){
				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				// partition data wrt worst and best value in dimension niveau
				for (int i=0; i<dataset.size(); i++){
					if (dataset[i][cfg->statDim_size+1+niveau]==best_value && !notSkyline[dataset[i][0]]){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, dataset[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						p[cfg->statDim_size+1+niveau]=0;
						branch_dataset.push_back(p);
					}else if (dataset[i][cfg->statDim_size+1+niveau]==worst_value){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, dataset[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						p[cfg->statDim_size+1+niveau]=1;
						branch_dataset.push_back(p);
					}
				}
			}else{
				for (auto tuple : dataset){
					if(tuple[1+cfg->statDim_size+niveau]==best_value && !notSkyline[tuple[0]]){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, tuple, (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						branch_dataset.push_back(p);
					}
				}
			}

			this->compute_view_recursively_md(cfg,  niveau+1, branch_dataset, preference_orders, notSkyline);
			/////
			//destroy Point pointers
			for ( auto p : branch_dataset)
			delete p;
		}

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//if we are processing last dimension
		if (niveau==cfg->dyDim_size-1){

			vector<Point> branch_dataset;
			if (best_value!=worst_value){
				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				// partitionner les donnees par rapport à cette dimension
				for (int i=0; i<dataset.size(); i++){
					if (dataset[i][cfg->statDim_size+1+niveau]==best_value && !notSkyline[dataset[i][0]]){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, dataset[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						p[cfg->statDim_size+1+niveau]=0;
						branch_dataset.push_back(p);
					}else if (dataset[i][cfg->statDim_size+1+niveau]==worst_value && !notSkyline[dataset[i][0]]){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, dataset[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						p[cfg->statDim_size+1+niveau]=1;
						branch_dataset.push_back(p);
					}
				}
			}else{
				for (auto tuple : dataset){
					if(tuple[1+cfg->statDim_size+niveau]==best_value){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, tuple, (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						branch_dataset.push_back(p);
					}
				}
			}

			//*************************************
			// calculer le skyline
			// vector<id> sky;
			if(branch_dataset.size()>0){
				int All = (1<<(cfg->statDim_size+cfg->dyDim_size))-1;
				vector<Space> full_Space;
				listeAttributsPresents(All, cfg->statDim_size+cfg->dyDim_size, full_Space);
				//sky=subspaceSkylineSize_TREE(full_Space, branch_dataset);
				ExecuteBSkyTree_bis(full_Space, branch_dataset, notSkyline);
			}

			// sort(sky.begin(),sky.end());
			
			// //MINUS : get non skyline points in branch_dataset
			// vector<id> branch_dataset_ids(branch_dataset.size());
			// for (int b=0;b<branch_dataset.size();b++){
			// 	branch_dataset_ids[b]=branch_dataset[b][0];
			// }
			// sort(branch_dataset_ids.begin(),branch_dataset_ids.end());

			// std::vector<id> notSky(branch_dataset_ids.size());   
		 //  	std::vector<id>::iterator it4;
		 //  	it4=std::set_difference(branch_dataset_ids.begin(), branch_dataset_ids.end(), sky.begin(), sky.end(), notSky.begin());                        
		 //  	notSky.resize(it4-notSky.begin());
		 //  	for (auto id: notSky) {
		 //  		notSkyline[id]=true;
		 //  	}
			/////
			//destroy Point pointers
			for ( auto p : branch_dataset)
			delete p;
			//	

		}
		
	}	

}



vector<id> dySky_v::compute_skyline(Config* cfg, vector<vector<Order>> preference_orders){
	cout << "dySky::compute_skyline" <<endl;

	bool* notSkyline=new bool[cfg->dataset_size];
	for(int i=0;i<cfg->dataset_size;++i){
		notSkyline[i]=false;
	}

	this->compute_views(cfg, preference_orders, notSkyline);

	for (int id_tuple : this->never_sky) notSkyline[id_tuple]=true;			

	// compute positive skyline

	vector<id> skyline_result=vector<id>(cfg->dataset_size);	
	int index=0;
	for (int i=0; i<cfg->dataset_size; i++){
		if (notSkyline[i]==false){ // means i is skyline
			skyline_result[index]=i;
			index++;
		}
	}
	delete[] notSkyline;
	skyline_result.resize(index);

	return skyline_result;
}

// void dySky_v::compute_view_1d(Config* cfg, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, uint64_t *storage){
	
// 	#pragma omp parallel for schedule(dynamic)
// 	for(int i=0;i<preference_orders[0].size();i++){
// 		int best_value=preference_orders[0][i].first;
// 		int worst_value=preference_orders[0][i].second;
// 		//cout << best_value << ":" << worst_value <<endl;

		
// 		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 		// partitionner les donnees par rapport aux valeurs best_value, worst_value
// 		vector<Point> dominating;
// 		vector<Point> dominated;
// 		for (int i=0; i<dataset.size(); i++){
// 			if (dataset[i][cfg->statDim_size+1]==best_value){
// 				Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
// 				memcpy(p, dataset[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
// 				p[cfg->statDim_size+1]=0;
// 				dominating.push_back(p);
// 			}else if (dataset[i][cfg->statDim_size+1]==worst_value){
// 				Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
// 				memcpy(p, dataset[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
// 				p[cfg->statDim_size+1]=1;
// 				dominated.push_back(p);
// 			}
// 		}


// 		//*************************************
// 		// calculer le skyline

// 		vector<id> notSky;
// 		find_dominated(cfg, dominating, dominated, notSky);

// 	  	//*********************************
// 	  	// store notSky

// 	  	#pragma omp critical
// 	  	{
// 	  		sky_view[Order(best_value,worst_value)]=new order_tree;
// 	  	}
// 	  	sky_view[Order(best_value,worst_value)]->ids=notSky;
// 	  	*storage=*storage+notSky.size();
// 	}
// }

bool check_dominance(Config* cfg, Point tuple_dominated, Point tuple_dominating, int dim){
	if (tuple_dominating[dim]>tuple_dominated[dim]){
		return false;
	}
	if (dim==cfg->statDim_size+1){
		return true;
	}
	return check_dominance(cfg, tuple_dominated, tuple_dominating, dim+1);
}

void dySky_v::find_dominated(Config* cfg, vector<Point> &dominating, vector<Point> &dominated, vector<id> &notSky){

	for (auto tuple1: dominated){
		for (auto tuple2: dominating){
			if (check_dominance(cfg, tuple1, tuple2, 1)){
				notSky.push_back(tuple1[0]);
				break;
			}
		}
	}

}

