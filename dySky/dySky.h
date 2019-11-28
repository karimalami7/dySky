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

class dySky {
	public:
	vector<Point> to_dataset;
	vector<vector<int>> po_dataset;
	vector<map<int, vector<int>>> partition_ids;
	vector<id> never_sky;
	vector<id> candidates;
	unordered_map<Order, order_tree*, pairhash> sky_view;

	void generate_to_data(Config* cfg); 
	void generate_po_data(Config* cfg);
	int compute_candidates(Config* cfg);
	
	void print_dataset (Config* cfg);
	void print_dataset (vector<Point> data, string name, int d);
	dySky(Config *cfg);
	void compute_views(Config* cfg, vector<vector<Order>> preference_orders, uint64_t *storage);
	void compute_view_1d(Config* cfg, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, uint64_t *storage);
	void compute_view_recursively_md(Config* cfg, int niveau, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, uint64_t *storage);
	vector<id> compute_skyline(Config* cfg, vector<vector<Order>> preference);

};

dySky::dySky(Config *cfg){
	this->po_dataset=vector<vector<int>>(cfg->dataset_size);
}

void dySky::generate_to_data(Config* cfg){
	loadData("ANTI","",cfg->dataset_size, cfg->statDim_size, cfg->statDim_val, this->to_dataset);

	
} 

void dySky::generate_po_data(Config* cfg){
	for (int j=0;j<cfg->dyDim_size;j++) {
		for (int i=0; i<cfg->dataset_size; i++){
			this->po_dataset[i].push_back(rand()%cfg->dyDim_val);
		}
	}
}


void dySky::print_dataset (Config* cfg){

    string const nomFichier1("datasets/INDE-"+to_string(cfg->statDim_size)+"-"+to_string(cfg->dataset_size)+".csv");
    ofstream monFlux1(nomFichier1.c_str());

    for (int i = 0; i < cfg->dataset_size ; i++)
    {
        for (int j = 0; j <= cfg->statDim_size ; j++){
            monFlux1 << this->to_dataset[i][j]<<";";  
        }
        for (int j = 0; j < cfg->dyDim_size ; j++){
            monFlux1 << this->po_dataset[i][j]<< ";";  
        } 
        monFlux1 << endl;
    }
}

void dySky::print_dataset (vector<Point> data, string name, int d){

    //string  nomFichier1(string);
    ofstream monFlux1(name);

    for (int i = 0; i < data.size() ; i++)
    {
        for (int j = 0; j <= d ; j++){
            monFlux1 << data[i][j]<<";";  
        }

        monFlux1 << endl;
    }
}

int dySky::compute_candidates(Config* cfg){
	cout << "dySky::compute_candidates"<<endl;

	// cluster dataset depending on the value in the po dimenion
	map<vector<int>, vector<Point>> partitions;
	for (int i=0; i<cfg->dataset_size; i++){
		partitions[this->po_dataset[i]].push_back(this->to_dataset[i]);
	}

	// compute skyline for every partition -> always skyline + candidates
	int All = (1<<cfg->statDim_size)-1;
	vector<Space> full_Space;
	listeAttributsPresents(All, cfg->statDim_size, full_Space);
    for (auto it_partition=partitions.begin(); it_partition!=partitions.end(); it_partition++){
    	vector<id> Skyline;
    	Skyline=subspaceSkylineSize_TREE(full_Space, it_partition->second);
    	this->candidates.insert(this->candidates.end(),
    		Skyline.begin(), Skyline.end());
    }
    sort (this->candidates.begin(),this->candidates.end());   
    cout << "Candidate set size: "<<this->candidates.size()<<endl;

    vector<id> all_ids=vector<id>(cfg->dataset_size);
    for (int i=0; i<cfg->dataset_size;i++){
    	all_ids[i]=i;
    }
	this->never_sky=vector<id>(cfg->dataset_size);  
  	std::vector<id>::iterator it4;
  	it4=std::set_difference(all_ids.begin(), all_ids.end(), candidates.begin(), candidates.end(), this->never_sky.begin());                        
  	this->never_sky.resize(it4-this->never_sky.begin());
  	cout << "never_sky set size: "<<this->never_sky.size()<<endl;

}


void dySky::compute_view_1d(Config* cfg, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, uint64_t *storage){
	
	#pragma omp parallel for schedule(dynamic)
	for(int i=0;i<preference_orders[0].size();i++){
		int best_value=preference_orders[0][i].first;
		int worst_value=preference_orders[0][i].second;
		//cout << best_value << ":" << worst_value <<endl;

		
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// partitionner les donnees par rapport aux valeurs best_value, worst_value
		vector<Point> branch_dataset;
		for (int i=0; i<dataset.size(); i++){
			if (dataset[i][cfg->statDim_size+1]==best_value){
				Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
				memcpy(p, dataset[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
				p[cfg->statDim_size+1]=0;
				// p[0]=0;
				// memcpy(p[0], dataset[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
				branch_dataset.push_back(p);
			}else if (dataset[i][cfg->statDim_size+1]==worst_value){
				Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
				memcpy(p, dataset[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
				p[cfg->statDim_size+1]=1;
				// p[0]=1;
				// memcpy(p[0], dataset[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
				branch_dataset.push_back(p);
			}
		}

		//*************************************
		// calculer le skyline
		int All = (1<<(cfg->statDim_size+cfg->dyDim_size))-1;
		vector<Space> full_Space;
		listeAttributsPresents(All, cfg->statDim_size+cfg->dyDim_size, full_Space);
		vector<id> sky=subspaceSkylineSize_TREE(full_Space, branch_dataset);

	  	//******************************************
	  	// sort ids
	  	std::sort (sky.begin(),sky.end());

		//******************************************
		// compute dominated points
		vector<id> branch_dataset_ids(branch_dataset.size());
		for (int b=0;b<branch_dataset.size();b++){
			branch_dataset_ids[b]=branch_dataset[b][0];
		}
		//sort(branch_dataset_ids.begin(),branch_dataset_ids.end());
		std::vector<id> notSky(branch_dataset_ids.size());   
	  	std::vector<id>::iterator it4;
	  	it4=std::set_difference(branch_dataset_ids.begin(), branch_dataset_ids.end(), sky.begin(), sky.end(), notSky.begin());                        
	  	notSky.resize(it4-notSky.begin());

	  	//*********************************
	  	// store notSky

	  	#pragma omp critical
	  	{
	  		sky_view[Order(best_value,worst_value)]=new order_tree;
	  	}
	  	sky_view[Order(best_value,worst_value)]->ids=notSky;
	  	*storage=*storage+notSky.size();
	}
}

void dySky::compute_views(Config* cfg, vector<vector<Order>> preference_orders, uint64_t *storage){
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
		compute_view_1d(cfg, candidates_tuples, this->sky_view, preference_orders, storage);
	}
	else{ // case of multiple dimensions

		for(int s=0;s<preference_orders[0].size();s++){
			this->sky_view[Order(preference_orders[0][s].first,preference_orders[0][s].second)]=new order_tree;
		}
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
					if (candidates_tuples[i][cfg->statDim_size+1]==best_value){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, candidates_tuples[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						p[1+cfg->statDim_size+niveau]=0;
						branch_dataset.push_back(p);
					}else if (candidates_tuples[i][cfg->statDim_size+1]==worst_value){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, candidates_tuples[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						p[1+cfg->statDim_size+niveau]=1;
						branch_dataset.push_back(p);
					}
				}
			}	
			else{
				for (auto tuple : candidates_tuples){
					if(tuple[1+cfg->statDim_size+niveau]==best_value){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, tuple, (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						branch_dataset.push_back(p);
					}
				}
			}
			compute_view_recursively_md(cfg, niveau+1, branch_dataset, this->sky_view[Order(best_value,worst_value)]->order_child, preference_orders, storage);
			
			/////
			//destroy Point pointers
			for ( auto p : branch_dataset)
			delete p;				
			
		}
	}
}

void dySky::compute_view_recursively_md(Config* cfg, int niveau, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, uint64_t *storage){
	
	#pragma omp parallel for schedule(dynamic) //if (niveau==1)
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
					if (dataset[i][cfg->statDim_size+1+niveau]==best_value){
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
					if(tuple[1+cfg->statDim_size+niveau]==best_value){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, tuple, (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						branch_dataset.push_back(p);
					}
				}
			}

			if (sky_view.find(Order(best_value,worst_value))==sky_view.end()){
				sky_view[Order(best_value,worst_value)]=new order_tree;;
			}
			this->compute_view_recursively_md(cfg,  niveau+1, branch_dataset, sky_view[Order(best_value,worst_value)]->order_child, preference_orders, storage);
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
					if (dataset[i][cfg->statDim_size+1+niveau]==best_value){
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
					if(tuple[1+cfg->statDim_size+niveau]==best_value){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, tuple, (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						branch_dataset.push_back(p);
					}
				}
			}

			//*************************************
			// calculer le skyline
			vector<id> sky;
			if(branch_dataset.size()>0){
				int All = (1<<(cfg->statDim_size+cfg->dyDim_size))-1;
				vector<Space> full_Space;
				listeAttributsPresents(All, cfg->statDim_size+cfg->dyDim_size, full_Space);
				sky=subspaceSkylineSize_TREE(full_Space, branch_dataset);
			}

			sort(sky.begin(),sky.end());
			
			//MINUS : get non skyline points in branch_dataset
			vector<id> branch_dataset_ids(branch_dataset.size());
			for (int b=0;b<branch_dataset.size();b++){
				branch_dataset_ids[b]=branch_dataset[b][0];
			}
			sort(branch_dataset_ids.begin(),branch_dataset_ids.end());

			std::vector<id> notSky(branch_dataset_ids.size());   
		  	std::vector<id>::iterator it4;
		  	it4=std::set_difference(branch_dataset_ids.begin(), branch_dataset_ids.end(), sky.begin(), sky.end(), notSky.begin());                        
		  	notSky.resize(it4-notSky.begin());

			/////
			//destroy Point pointers
			for ( auto p : branch_dataset)
			delete p;
			//	

	  		sky_view[Order(best_value,worst_value)]=new order_tree;
	  		sky_view[Order(best_value,worst_value)]->ids.swap(notSky);
		  	
		  	*storage=*storage+sky_view[Order(best_value,worst_value)]->ids.size();
		}
		
	}	

}



vector<id> dySky::compute_skyline(Config* cfg, vector<vector<Order>> preference_cross){
	cout << "dySky::compute_skyline" <<endl;
	// cerr << "dySky::compute_skyline" <<endl;
	// for (int i=0;i<preference_cross.size();i++){
	// 	for (int j=0;j<preference_cross[i].size();j++){
	// 		cout << preference_cross[i][j].first<< ":"<< preference_cross[i][j].second <<endl;
	// 	}
	// 	cout<<endl;
	// }


	bool* cSkyline=new bool[cfg->dataset_size];
	for (int i=0;i<cfg->dataset_size;i++) cSkyline[i]=true;

	for (int i=0;i<preference_cross.size();i++){ // loop on orders
		order_tree *ot=sky_view[preference_cross[i][0]];
		//cout << "niv: 0 " <<preference_cross[i][0].first<< ":"<< preference_cross[i][0].second <<endl;
		int j=1;
		while (j<cfg->dyDim_size){ // enter here when multi po dimension
			ot=ot->order_child[preference_cross[i][j]];
			//cout << "niv:"<<j <<" "<< preference_cross[i][j].first<< ":"<< preference_cross[i][j].second <<endl;
			j++;
		}
		
		for ( int id_tuple : ot->ids) cSkyline[id_tuple]=false;		
	}
	// remove never_sky 
	
	for (int id_tuple : this->never_sky) cSkyline[id_tuple]=false;			

	// compute positive skyline

	vector<id> skyline_result=vector<id>(cfg->dataset_size);	
	int index=0;

	for (int i=0; i<cfg->dataset_size; i++){
		if (cSkyline[i]==true){
			skyline_result[index]=i;
			index++;
		}
		
	}
	delete[] cSkyline;
	skyline_result.resize(index);

	return skyline_result;
}



