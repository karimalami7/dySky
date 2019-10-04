/*
 * dySky_h.h
 *
 *  Created on: septembre 20, 2019
 *      Author: karim
 *
 * Functionalities:  hybrid dySky
 */

#include "../common/common.h"
#include "../common/config.h"
#include "../generator/generateur.h"
#include "../BSkyTree/bskytree.h"
#include <omp.h>
using namespace std;

class dySky_h: public dySky {
	public:

	map<vector<Order>, vector<id>> selected_spos;

	dySky_h(Config *cfg);
	
	void views_selection(Config* cfg, uint64_t max_storage, vector<Query> &workload);

	vector<id> hybrid_compute_skyline(Config* cfg, vector<vector<Order>> preference_cross);
	void hybrid_compute_view_1d(Config* cfg, vector<Point> &dataset, Order o, bool* cSkyline);
	void compute_skyline_wrt_missing_spos(Config* cfg, vector<vector<Order>> missing_spos, bool* cSkyline);
	void hybrid_compute_view_recursively_md(Config* cfg, int niveau, vector<Point> &dataset, vector<Order> spo, bool* cSkyline);
	int knapSack(int W, int wt[], int val[], int n, vector <int> &selected_spos); 
	void comb(int N, vector<vector<int>> &all_combination);
};



dySky_h::dySky_h(Config *cfg) : dySky(cfg){}



void dySky_h::hybrid_compute_view_recursively_md(Config* cfg, int niveau, vector<Point> &dataset, vector<Order> spo, bool* cSkyline){

	int best_value=spo[niveau].first;
	int worst_value=spo[niveau].second;
	//cout <<"niveau: "<<niveau<< " --> "<< best_value<<" : "<<worst_value<<endl;		


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//si c'est pas la dernière dimension 				
	if (niveau<cfg->dyDim_size -1){
		vector<Point> branch_dataset;
		if(best_value!=worst_value){
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
		}
		else{
			for (auto tuple : dataset){
				if(tuple[1+cfg->statDim_size+niveau]==best_value){
					Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
					memcpy(p, tuple, (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
					branch_dataset.push_back(p);
				}
			}

		}

		this->hybrid_compute_view_recursively_md(cfg, niveau+1, branch_dataset, spo, cSkyline);
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//si c'est la dernière dimension
	if (niveau==cfg->dyDim_size-1){

		vector<Point> branch_dataset;
		if(best_value!=worst_value){
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
		}
		else{
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

	  	//*********************************
	  	// fill cSkyline

	  	for (int id_tuple : notSky) cSkyline[id_tuple]=false;	

	}	
	
}





void dySky_h::hybrid_compute_view_1d(Config* cfg, vector<Point> &dataset, Order o, bool* cSkyline){
	

	int best_value=o.first;
	int worst_value=o.second;
	//cout << best_value << ":" << worst_value <<endl;

	
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// partitionner les donnees par rapport aux valeurs best_value, worst_value
	vector<Point> branch_dataset;
	for (int i=0; i<dataset.size(); i++){
		if (dataset[i][cfg->statDim_size+1]==best_value){
			Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
			memcpy(p, dataset[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
			p[cfg->statDim_size+1]=0;
			branch_dataset.push_back(p);
		}else if (dataset[i][cfg->statDim_size+1]==worst_value){
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
  	// fill cSkyline

  	for (int id_tuple : notSky) cSkyline[id_tuple]=false;	

}









void dySky_h::compute_skyline_wrt_missing_spos(Config* cfg, vector<vector<Order>> missing_spos, bool* cSkyline){

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

	if(cfg->dyDim_size==1){
		#pragma omp parallel for schedule(dynamic)
		for (int i=0; i<missing_spos.size(); i++){
			hybrid_compute_view_1d(cfg, candidates_tuples, missing_spos[i][0], cSkyline);
		}		
	}
	else{
		#pragma omp parallel for schedule(dynamic)
		for (int i=0; i<missing_spos.size(); i++){
			int best_value=missing_spos[i][0].first;
			int worst_value=missing_spos[i][0].second;
			//cout <<"niveau: "<<niveau<< " --> "<< best_value<<" : "<<worst_value<<endl;		
			vector<Point> branch_dataset;
			if (best_value!=worst_value){
				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				// partitionner les donnees par rapport à cette dimension
				for (int i=0; i<candidates_tuples.size(); i++){
					if (candidates_tuples[i][cfg->statDim_size+1]==best_value){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, candidates_tuples[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						p[1+cfg->statDim_size]=0;
						branch_dataset.push_back(p);
					}else if (candidates_tuples[i][cfg->statDim_size+1]==worst_value){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, candidates_tuples[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						p[1+cfg->statDim_size]=1;
						branch_dataset.push_back(p);
					}
				}	
			}else {
				for (auto tuple : candidates_tuples){
					if(tuple[1+cfg->statDim_size]==best_value){
						Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
						memcpy(p, tuple, (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
						branch_dataset.push_back(p);
					}
				}
			}

			hybrid_compute_view_recursively_md(cfg, 1, branch_dataset, missing_spos[i], cSkyline);
				
			/////
			//destroy Point pointers
			for ( auto p : branch_dataset)
			delete p;				
		}
	}
}

vector<id> dySky_h::hybrid_compute_skyline(Config* cfg, vector<vector<Order>> preference_cross){
	cout << "dySky::compute_skyline" <<endl;

	// for (int i=0;i<preference_cross.size();i++){
	// 	for (int j=0;j<preference_cross[i].size();j++){
	// 		cout << preference_cross[i][j].first<< ":"<< preference_cross[i][j].second <<endl;
	// 	}
	// 	cout<<endl;
	// }
	vector<vector<Order>> missing_spos;

	bool* cSkyline=new bool[cfg->dataset_size];
	for (int i=0;i<cfg->dataset_size;i++) cSkyline[i]=true;

	for (auto query_spo : preference_cross){ // loop on orders

		auto it=this->selected_spos.find(query_spo);
		if(it!=this->selected_spos.end()){
			for (int id_tuple : it->second) cSkyline[id_tuple]=false;
		}else{
			missing_spos.push_back(query_spo);
		}				
	}

	//here compute c-skyline wrt missing spos

	cout << "# missing spos to compute: "<< missing_spos.size()<<endl;

	// remove never_sky 

	compute_skyline_wrt_missing_spos(cfg, missing_spos,cSkyline);
	
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

void dySky_h::views_selection(Config* cfg, uint64_t max_storage, vector<Query> &workload){
	
	
	//generate workload

	// vector<Query> workload(cfg->dyDim_val*3);
	// for (int q=0; q<workload.size(); q++){
	// 	workload[q].generate_preference(cfg);
	// 	workload[q].graph_to_orderPairs(cfg);
	// 	workload[q].cross_orders_over_dimensions(cfg);
	// }

	// compute gain of each spo
	map<vector<Order>, int> gain;
	for (int q=0; q<workload.size(); q++){
		// boucle sur les spos
		for (auto spo: workload[q].preference_orders_cross){
			if (gain.find(spo)!=gain.end()){
				gain[spo]++;
			}
			else{
				gain[spo]=1;
			}
		} 
	}

	cout <<"total number of spos: "<<gain.size()<<endl;

	// fill weight and value arrays

	int *wt = new int[gain.size()];
	int *val = new int[gain.size()];

	int i=0;
	for (auto it=gain.begin(); it!=gain.end(); it++){

		auto ot=sky_view[it->first[0]];
		int j=1;
		while (j<cfg->dyDim_size){
			ot=ot->order_child[it->first[j]];
			j++;
		}
			wt[i]=ot->ids.size();
			val[i]=it->second;
		i++;
	}

	int total_required_memory=0;
	int total_gain=0;

	cout << "spos, their weight and gain"<<endl; 
	for (int i=0; i<gain.size();i++){
		cout << "spo " <<i <<" : "<<wt[i] << "  " << val[i] <<endl;
		total_required_memory+=wt[i];
		total_gain+=val[i];
	}

	cout << "Memory required to store all spos: "<< total_required_memory<<endl;
	cout << "Maximum Gain: "<< total_gain<<endl;

	
	vector<int> selected_spos_id;
	knapSack(max_storage, wt, val, gain.size(), selected_spos_id) ;


	// keep selected spos


	for (auto id_spo : selected_spos_id){
		
		auto it=gain.begin();
		for(int i=0;i<id_spo;i++) it++;

		auto ot=sky_view[it->first[0]];
		int j=1;
		while (j<cfg->dyDim_size){
			ot=ot->order_child[it->first[j]];
			j++;
		}

		this->selected_spos.insert(pair<vector<Order>, vector<id>>(it->first,ot->ids));
		delete ot;
	}
}



// // A utility function that returns maximum of two integers 
// int max(int a, int b) { return (a > b)? a : b; } 

// Returns the maximum value that can be put in a knapsack of capacity W 
int dySky_h::knapSack(int W, int wt[], int val[], int n, vector<int> &spo_ids) 
{ 	
	cout << "****************************"<<endl;
	cout << "knapsack Algo starts "<<endl;
	cout << "Available memory resources: "<< W<<endl;

	int i, w; 
	//int K[n+1][W+1]; 
	//bool keep[n+1][W+1];

	int* K= new int[W+1];
	int* K_temp= new int[W+1];
	for (i=0; i< W+1; i++) K_temp[i]=0;
	//for(int i=0; i<n+1; i++) K[i]=new int[W+1];
	// vector<bool*> keep= vector<bool*>(n+1);
	// for(int i=0; i<n+1; i++) keep[i]=new bool[W+1];

	// for (i = 0; i <= n; i++) for (w = 0; w <= W; w++) keep[i][w]=0;
	vector<set<int>> keep(n+1);
	// Build table K[][] in bottom up manner 
	for (i = 0; i <= n; i++) 
	{ 
		for (w = 0; w <= W; w++) 
		{ 
			if (w==0) 
				K[w] = 0; 
			else if (wt[i-1] <= w && val[i-1] + K_temp[w-wt[i-1]] > K_temp[w]){ 
				K[w] = val[i-1] + K_temp[w-wt[i-1]];
				keep[i].insert(w); 
			}	
			else{
				K[w] = K_temp[w];
				//keep[i][w]=0;
			}
		} 
		for (int j=0; j< W+1; j++) K_temp[j]=K[j];
	} 

	delete K; delete K_temp;

	//for (auto k : K) delete k;

	cout << "Selected spo: "<<endl;

	int total_gain=0;
	int total_poids=0;

	int X=W;
	for (i=n; i>=1; i--){
		if (keep[i].find(X)!=keep[i].end()){
			cout << "spo: " << i-1 << ",gain: "<< val[i-1]<< ",weight: "<< wt[i-1]<<endl;
			spo_ids.push_back(i-1);
			total_poids+=wt[i-1];
			total_gain+=val[i-1];
			X=X-wt[i-1];
		}
	}

	cout <<"Gain total: "<<total_gain<<", Poids total: " <<total_poids<<endl;

	return K[W]; 
} 

void dySky_h::comb(int N, vector<vector<int>> &all_combination)
{	
	for (int K=1; K<N; K++){

	    std::string bitmask(K, 1); // K leading 1's
	    bitmask.resize(N, 0); // N-K trailing 0's

	    // print integers and permute bitmask
	    
	    do {
	    	vector <int> combin;
	        for (int i = 0; i < N; ++i) // [0..N-1] integers
	        {
	            if (bitmask[i]) combin.push_back(i);
	        }
	        all_combination.push_back(combin);
	    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
	}
}