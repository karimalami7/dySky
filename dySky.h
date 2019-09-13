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

#include "common.h"
#include "config.h"
#include "generator/generateur.h"
#include "BSkyTree/bskytree.h"
using namespace std;

class dySky {
	public:
	vector<Point> to_dataset;
	vector<vector<int>> po_dataset;
	vector<map<int, vector<int>>> partition_ids;
	//vector<id> always_sky;
	vector<id> never_sky;
	vector<id> candidates;
	unordered_map<Order, order_tree*, pairhash> sky_view;

	void generate_to_data(Config* cfg); 
	void generate_po_data(Config* cfg);
	int compute_candidates(Config* cfg);
	
	void print_dataset (Config* cfg);
	void print_dataset (vector<Point> data, string name, int d);
	dySky(Config *cfg);
	void recur_intersect(Config *cfg, vector<id> &dataset, int niveauf, vector<int> f, vector<int> l, vector<id> &nsky_for_all_orders_global, vector<Order> orders_stack);
	void recur_preparation(Config *cfg, vector<Point> &to, vector<id> &nsky_for_all_orders_global, int niveaul, vector<int> l, vector<Order> orders_stack);
	void views_selection(Config* cfg);
	void compute_views(Config* cfg, vector<vector<Order>> preference_orders, int *storage);
	void compute_view_1d(Config* cfg, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, int *storage);
	void compute_view_recursively_md(Config* cfg, int niveau, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, int *storage);
	void compute_other_nsky(Config *cfg, vector<id> &nsky_for_all_orders_global, vector<id> &dataset, vector<Order> orders_stack);
	vector<id> compute_skyline(Config* cfg, vector<vector<Order>> preference);
	void compute_to_skyline(Config *cfg, vector<Point> &to, vector<id> &nsky_for_all_orders_global, int niveaul);
};

int knapSack(int W, int wt[], int val[], int n); 
void comb(int N, vector<vector<int>> &all_combination);

dySky::dySky(Config *cfg){
	this->po_dataset=vector<vector<int>>(cfg->dataset_size);
}

void dySky::generate_to_data(Config* cfg){
	loadData("INDE","",cfg->dataset_size, cfg->statDim_size, cfg->statDim_val, this->to_dataset);

	
} 

void dySky::generate_po_data(Config* cfg){
	for (int j=0;j<cfg->dyDim_size;j++) {
		for (int i=0; i<cfg->dataset_size; i++){
			this->po_dataset[i].push_back(rand()%cfg->dyDim_val);
		}
	}
}


void dySky::print_dataset (Config* cfg){

    string const nomFichier1("INDE-"+to_string(cfg->statDim_size)+"-"+to_string(cfg->dataset_size)+".csv");
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
	// cout << "print ids" <<endl; 
	// for (auto tuple : this->candidates){
	// 	cout << tuple << endl;
	// }
    vector<id> all_ids=vector<id>(cfg->dataset_size);
    for (int i=0; i<cfg->dataset_size;i++){
    	all_ids[i]=i;
    }
	this->never_sky=vector<id>(cfg->dataset_size);  
  	std::vector<id>::iterator it4;
  	it4=std::set_difference(all_ids.begin(), all_ids.end(), candidates.begin(), candidates.end(), this->never_sky.begin());                        
  	this->never_sky.resize(it4-this->never_sky.begin());
  	cout << "never_sky set size: "<<this->never_sky.size()<<endl;

	//cout << "print ids of never_sky" <<endl; 
	//for (auto tuple : this->never_sky){
	//	cout << tuple << endl;
	//}

 	// string const nomFichier1("candidates.csv");
  //   ofstream monFlux1(nomFichier1.c_str());

  //   for (int i = 0; i < this->candidates.size() ; i++)
  //   {
  //       for (int j = 0; j <= cfg->statDim_size ; j++){
  //           monFlux1 << this->to_dataset[this->candidates[i]][j]<<";";  
  //       } 
  //       monFlux1 << this->po_dataset[this->candidates[i]][0]<< endl;
  //   }

}

void dySky::views_selection(Config* cfg){
	
	cout <<sky_view.size()<<endl;


	//generate workload

	vector<Query> workload(cfg->dyDim_val);
	for (int q=0; q<workload.size(); q++){
		workload[q].generate_preference(cfg);
		workload[q].graph_to_orderPairs(cfg);	
	}

	// compute gain of each order
	unordered_map<Order, int, pairhash> gain;
	for (int q=0; q<workload.size(); q++){
		for (int i=0;i<workload[q].preference_orders[0].size();i++){
			gain[workload[q].preference_orders[0][i]]++;
		}
	}

	// fill weight and value arrays

	int *wt = new int[gain.size()];
	int *val = new int[gain.size()];

	int i=0;
	for (auto it=gain.begin(); it!=gain.end(); it++){
		auto it2=sky_view.find(it->first);
		if(it2!=sky_view.end()){
			wt[i]=it2->second->ids.size();
			val[i]=it->second;
		}
		else{
			cout<<"error, view not found"<<endl;
		}
		i++;
	}

	for (int i=0; i<gain.size();i++){
		cout << i <<" : "<<wt[i] << "  " << val[i] <<endl;
	}

	cout << knapSack(10000, wt, val, gain.size()) << endl;

}



void dySky::compute_view_1d(Config* cfg, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, int *storage){
	
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
	  	// store notSky

	  	#pragma omp critical
	  	{
	  		sky_view[Order(best_value,worst_value)]=new order_tree;
	  	}
	  	sky_view[Order(best_value,worst_value)]->ids=notSky;
	  	*storage=*storage+notSky.size();
	}
}

void dySky::compute_views(Config* cfg, vector<vector<Order>> preference_orders, int *storage){
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
	else{ // case of multi-dimensions

		//*******************************************
		// compute all partitions for all dimensions, l * m vector

		this->partition_ids=vector<map<int, vector<int>>>(cfg->dyDim_size);
		for (int i=0; i<this->candidates.size(); i++){
			for (int j=0; j<cfg->dyDim_size; j++){
				partition_ids[j][this->po_dataset[this->candidates[i]][j]].push_back(this->candidates[i]);
			}	
		}


		int niveau=0;
		//#pragma omp parallel for schedule(dynamic)
		for(int s=0;s<preference_orders[0].size();s++){
			int best_value=preference_orders[0][s].first;
			int worst_value=preference_orders[0][s].second;
			//cout <<"niveau: "<<niveau<< " --> "<< best_value<<" : "<<worst_value<<endl;		
			
			
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			// partitionner les donnees par rapport à cette dimension
			vector<Point> branch_dataset;
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
			#pragma omp critical
			{
				this->sky_view[Order(best_value,worst_value)]=new order_tree;
			}		
			compute_view_recursively_md(cfg, niveau+1, branch_dataset, this->sky_view[Order(best_value,worst_value)]->order_child, preference_orders, storage);
			
			/////
			//destroy Point pointers
			for ( auto p : branch_dataset)
			delete p;
			//	
		}

		//here, we parition by values in dimension 0 and we call the recursive function Order (-1,-1)
 		
		//cout <<"niveau: "<<niveau<< " --> "<< "-1"<<" : "<<"-1"<<endl;

		map<int,vector<Point>> partitions;
		for (auto tuple : candidates_tuples){
			Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
			memcpy(p, tuple, (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
			partitions[p[1+cfg->statDim_size+niveau]].push_back(p);
		}

		this->sky_view[Order(-1,-1)]=new order_tree;


		for (auto partition : partitions){
			// cout << "partition size: "<< partition.second.size()<<endl;
	  // 	 	cout << "print ids of the partition: " <<endl; 
			// for (auto tuple : partition.second){
			// 	cout << tuple[0] << endl;
			// }
			compute_view_recursively_md(cfg, niveau+1, partition.second, this->sky_view[Order(-1,-1)]->order_child, preference_orders, storage);
		}



		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// add other nsky
		// vector<id> nsky_for_all_orders_global;
		// compute_other_nsky(cfg, nsky_for_all_orders_global, this->candidates, orders_stack[0]);
	}

}

void dySky::compute_view_recursively_md(Config* cfg, int niveau, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, int *storage){
	
	//#pragma omp parallel for schedule(dynamic) if (niveau==0)
	for(int s=0;s<preference_orders[niveau].size();s++){
		int best_value=preference_orders[niveau][s].first;
		int worst_value=preference_orders[niveau][s].second;
		//cout <<"niveau: "<<niveau<< " --> "<< best_value<<" : "<<worst_value<<endl;		

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//si c'est pas la dernière dimension 				
		if (niveau<cfg->dyDim_size -1){
			vector<Point> branch_dataset;
			
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
			if (sky_view.find(Order(best_value,worst_value))==sky_view.end()){
				sky_view[Order(best_value,worst_value)]=new order_tree;;
			}
			this->compute_view_recursively_md(cfg,  niveau+1, branch_dataset, sky_view[Order(best_value,worst_value)]->order_child, preference_orders, storage);
		}

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//si c'est la dernière dimension
		if (niveau==cfg->dyDim_size-1){

			vector<Point> branch_dataset;
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

			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			// add other nsky
			// vector<id> nsky_for_all_orders_global;
			// compute_other_nsky(cfg, nsky_for_all_orders_global, this->candidates, orders_stack);
			
			//UNION 
			// std::vector<id> h(nsky_for_all_orders_global.size()+notSky.size());   
		 //  	std::vector<id>::iterator it3;
		 //  	it3=std::set_union(notSky.begin(), notSky.end(), nsky_for_all_orders_global.begin(), nsky_for_all_orders_global.end(), h.begin());                        
		 //  	h.resize(it3-h.begin());   
		 //  	notSky=h; 

		  	if (sky_view.find(Order(best_value,worst_value))==sky_view.end()){
		  		//cout << "order not found" <<endl;
		  		sky_view[Order(best_value,worst_value)]=new order_tree;
		  		sky_view[Order(best_value,worst_value)]->ids=notSky;
		  	}
		  	else{
		  		// cout << "order found" <<endl;
		  		sky_view[Order(best_value,worst_value)]->ids.insert(sky_view[Order(best_value,worst_value)]->ids.end(),notSky.begin(),notSky.end());
		  		sort(sky_view[Order(best_value,worst_value)]->ids.begin(),sky_view[Order(best_value,worst_value)]->ids.end());
		  	}

	  // 	 	cout << "print ids of notSky: " <<endl; 
			// for (auto tuple : sky_view[Order(best_value,worst_value)]->ids){
			// 	cout << tuple << endl;
			// }
		  	
		  	*storage=*storage+notSky.size();
		}
	}	

	//here, we parition by values in dimension 0 and we call the recursive function Order (-1,-1)

	// cout <<"niveau: "<<niveau<< " --> "<< "-1"<<" : "<<"-1"<<endl;

	if (niveau<cfg->dyDim_size -1){ // if not last dimension
		map<int,vector<Point>> partitions;
		for (auto tuple : dataset){
			Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
			memcpy(p, tuple, (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
			partitions[p[1+cfg->statDim_size+niveau]].push_back(p);
		}


		

		for (auto partition : partitions){

			// cout << "partition size: "<< partition.second.size()<<endl;
	  // 	 	cout << "print ids of the partition: " <<endl; 
			// for (auto tuple : partition.second){
			// 	cout << tuple[0] << endl;
			// }
			if (sky_view.find(Order(-1,-1))==sky_view.end()){
				sky_view[Order(-1,-1)]=new order_tree;
				compute_view_recursively_md(cfg, niveau+1, partition.second, sky_view[Order(-1,-1)]->order_child, preference_orders, storage);
			}
			else
			{
				compute_view_recursively_md(cfg, niveau+1, partition.second, sky_view[Order(-1,-1)]->order_child, preference_orders, storage);
			}
		}
	}

	if (niveau==cfg->dyDim_size-1){ // if last dimension
		map<int,vector<Point>> partitions;
		for (auto tuple : dataset){
			Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
			memcpy(p, tuple, (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
			partitions[p[1+cfg->statDim_size+niveau]].push_back(p);
		}
		vector<id> notSky;
		for (auto partition : partitions){ // compute skyline for each partition

			//*************************************
			// calculer le skyline
			vector<id> sky_partition;
			if(partition.second.size()>0){
				int All = (1<<(cfg->statDim_size+cfg->dyDim_size))-1;
				vector<Space> full_Space;
				listeAttributsPresents(All, cfg->statDim_size+cfg->dyDim_size, full_Space);
				sky_partition=subspaceSkylineSize_TREE(full_Space, partition.second);				
			}

			sort(sky_partition.begin(),sky_partition.end());

			//MINUS : get non skyline points in branch_dataset
			vector<id> partition_ids(partition.second.size());
			for (int b=0;b<partition.second.size();b++){
				partition_ids[b]=partition.second[b][0];
			}
			sort(partition_ids.begin(),partition_ids.end());

			std::vector<id> notSky_partition(partition_ids.size());   
		  	std::vector<id>::iterator it4;
		  	it4=std::set_difference(partition_ids.begin(), partition_ids.end(), sky_partition.begin(), sky_partition.end(), notSky_partition.begin());                        
		  	notSky_partition.resize(it4-notSky_partition.begin());


			notSky.insert(notSky.end(), notSky_partition.begin(), notSky_partition.end());	
		}

		sort(notSky.begin(),notSky.end());
		
		if (sky_view.find(Order(-1,-1))==sky_view.end()){
			sky_view[Order(-1,-1)]=new order_tree;
			sky_view[Order(-1,-1)]->ids=notSky;			
		}
		else{
	  		sky_view[Order(-1,-1)]->ids.insert(sky_view[Order(-1,-1)]->ids.end(),notSky.begin(),notSky.end());
	  		sort(sky_view[Order(-1,-1)]->ids.begin(),sky_view[Order(-1,-1)]->ids.end());				
		}

  	 	
  	 	// cout << "print ids of notSky: " <<endl; 
		// for (auto tuple : notSky){
		// 	cout << tuple << endl;
		// }
		
		*storage=*storage+notSky.size();
	}

}

void dySky::recur_intersect(Config *cfg, vector<id> &dataset, int niveauf, vector<int> f, vector<int> l, vector<id> &nsky_for_all_orders_global, vector<Order> orders_stack){
	
	
	for (int i=0 ; i<this->partition_ids[f[niveauf]].size(); i++){

		// compute intersection between dataset and the current partition
		vector<id> intersection(dataset.size()+partition_ids[f[niveauf]][i].size());
		vector<id>::iterator it;
		it=std::set_intersection(dataset.begin(), dataset.end(), partition_ids[f[niveauf]][i].begin(), partition_ids[f[niveauf]][i].end(), intersection.begin());
		intersection.resize(it-intersection.begin());

		if (niveauf+1 < f.size()){ // if there is again another dimension, call recur_intersect
			recur_intersect(cfg, intersection, niveauf+1, f, l, nsky_for_all_orders_global, orders_stack);	
		} 
		else{
			int niveaul=0;
			vector<Point> to;
			for (int e=0; e<intersection.size(); e++){

				if (this->po_dataset[intersection[e]][l[niveaul]]==orders_stack[l[niveaul]].first){
					Point p=(int*)malloc((cfg->statDim_size+2)*sizeof(int));
					for (int j=0;j<=cfg->statDim_size;j++) p[j]=this->to_dataset[intersection[e]][j];
					p[cfg->statDim_size+1]=0;
					to.push_back(p);

				}else if (this->po_dataset[intersection[e]][l[niveaul]]==orders_stack[l[niveaul]].second){
					Point p=(int*)malloc((cfg->statDim_size+2)*sizeof(int));
					for (int j=0;j<=cfg->statDim_size;j++) p[j]=this->to_dataset[intersection[e]][j];
					p[cfg->statDim_size+1]=1;
					to.push_back(p);
				}
			}
			if (niveaul+1 < l.size()){
				recur_preparation(cfg, to, nsky_for_all_orders_global, niveaul, l, orders_stack);
			}
			else{
				compute_to_skyline(cfg, to, nsky_for_all_orders_global, niveaul);
			}
		}
	}
}

void dySky::recur_preparation(Config *cfg, vector<Point> &to_dataset, vector<id> &nsky_for_all_orders_global, int niveaul, vector<int> l, vector<Order> orders_stack){

	vector<Point> to;
	for (int i=0; i<to_dataset.size(); i++){
		
		id tuple_id= to_dataset[i][0];
		int dim=l[niveaul];

		if (this->po_dataset[tuple_id][dim]==orders_stack[dim].first){
			Point tuple=(int*)malloc((1+cfg->statDim_size+1+niveaul)*sizeof(int));
			for (int j=0;j<cfg->statDim_size+1+niveaul;j++) tuple[j]=to_dataset[i][j];
			tuple[cfg->statDim_size+1+niveaul]=0;
			to.push_back(tuple);

		}
		else if (this->po_dataset[tuple_id][dim]==orders_stack[dim].second){
			Point tuple=(int*)malloc((1+cfg->statDim_size+1+niveaul)*sizeof(int));
			for (int j=0;j<cfg->statDim_size+1+niveaul;j++) tuple[j]=to_dataset[i][j];
			tuple[cfg->statDim_size+1+niveaul]=1;
			to.push_back(tuple);
		}
	}

	/////
	//destroy Point pointers
	for ( auto p : to_dataset)
	delete p;
	//	
	
	if (niveaul+1 < l.size()){
		recur_preparation(cfg, to, nsky_for_all_orders_global, niveaul+1, l, orders_stack);
	}
	else{
		compute_to_skyline(cfg, to, nsky_for_all_orders_global, niveaul);
	}
}

void dySky::compute_other_nsky(Config *cfg, vector<id> &nsky_for_all_orders_global, vector<id> &dataset, vector<Order> orders_stack){	

		vector<int> f,l;

	vector<vector<int>> all_combination;
	comb(cfg->dyDim_size, all_combination);
	vector<int> all_dim(cfg->dyDim_size); for (int i=0; i<cfg->dyDim_size; i++) all_dim[i]=i;

	for (auto combination : all_combination){

		// **********************************************
		// for each combination of the set of dynamic dimensions, compute
		//
		f=combination;
		l=vector<int>(all_dim.size());
		vector<int>::iterator it;
		it=std::set_difference(all_dim.begin(), all_dim.end(), combination.begin(), combination.end(), l.begin());
		l.resize(it-l.begin());

        //***********************************************
		//  f= the set of dimensions to partition by equal values, l= the set of dimensions to partition by the order

		int niveauf=0;// start with first dimension in 
		for (int i=0 ; i<this->partition_ids[f[niveauf]].size(); i++){ // loop on partitions with respect to f[niveauf]
			vector<id> dataset_id= this->partition_ids[f[niveauf]][i];
			if (niveauf+1 < f.size()){ // if there is more than, then compute intersection with the others
				recur_intersect(cfg, dataset_id, niveauf+1, f, l, nsky_for_all_orders_global, orders_stack);	
			}else{ // go to second step

				int niveaul=0;
				vector<Point> to;
				for (int e=0; e<dataset_id.size(); e++){
					if (this->po_dataset[dataset_id[e]][l[niveaul]]==orders_stack[l[niveaul]].first){
						Point p=(int*)malloc((1+cfg->statDim_size+1)*sizeof(int));
						for (int j=0;j<cfg->statDim_size+1;j++) p[j]=this->to_dataset[dataset_id[e]][j];
						p[cfg->statDim_size+1]=0;
						to.push_back(p);
					}else if (this->po_dataset[dataset_id[e]][l[niveaul]]==orders_stack[l[niveaul]].second){
						Point p=(int*)malloc((cfg->statDim_size+2)*sizeof(int));
						for (int j=0;j<cfg->statDim_size+1;j++) p[j]=this->to_dataset[dataset_id[e]][j];
						p[cfg->statDim_size+1]=1;
						to.push_back(p);
					}
				}
				if (niveaul+1 < l.size()){
					recur_preparation(cfg, to, nsky_for_all_orders_global, niveaul+1, l, orders_stack);
				}
				else{
					compute_to_skyline(cfg, to, nsky_for_all_orders_global, niveaul);
				}
			}
		}
	}
}

void dySky::compute_to_skyline(Config *cfg, vector<Point> &to, vector<id> &nsky_for_all_orders_global, int niveaul){
	
	if (to.size()!=0){

		//******************************************
		// 1/ compute skyline for every partition
		int All = (1<<cfg->statDim_size+1+niveaul)-1;
		vector<Space> full_Space;
		listeAttributsPresents(All, cfg->statDim_size+1+niveaul, full_Space);
		vector<id> sky_for_all_orders;
    	vector<id> Skyline;
    	Skyline=subspaceSkylineSize_TREE(full_Space, to);
    	sort(Skyline.begin(),Skyline.end());

    	//**************************************
    	// 2 Skyline minus it_partition->second
    	vector<id> partition_ids(to.size());
    	for (int i=0; i<partition_ids.size();i++ ){
    		partition_ids[i]=to[i][0];
    	}
    	sort(partition_ids.begin(), partition_ids.end());
    	vector<id> notSky(to.size());  
	  	vector<id>::iterator it;
	  	it=std::set_difference(partition_ids.begin(), partition_ids.end(), Skyline.begin(), Skyline.end(), notSky.begin());                        
	  	notSky.resize(it-notSky.begin());  

	  	//**************************************
	  	// 3 notSky union nsky_for_all_orders_global

	  	vector<id> unionNsky(nsky_for_all_orders_global.size()+notSky.size());
	  	it=std::set_union(nsky_for_all_orders_global.begin(), nsky_for_all_orders_global.end(), notSky.begin(), notSky.end(), unionNsky.begin());
	  	unionNsky.resize(it-unionNsky.begin());
	  	nsky_for_all_orders_global=unionNsky;

	    sort(nsky_for_all_orders_global.begin(),nsky_for_all_orders_global.end());

    	/////
		//destroy Point pointers
		for ( auto p : to)
		delete p;
		//	
	}    
}


vector<id> dySky::compute_skyline(Config* cfg, vector<vector<Order>> preference_cross){
	cout << "dySky::compute_skyline" <<endl;

	// for (int i=0;i<preference_cross.size();i++){
	// 	for (int j=0;j<preference_cross[i].size();j++){
	// 		cout << preference_cross[i][j].first<< ":"<< preference_cross[i][j].second <<endl;
	// 	}
	// 	cout<<endl;
	// }

	vector<id> skyline_result=vector<id>(cfg->dataset_size);
    for (int i=0; i<cfg->dataset_size;i++){
    	skyline_result[i]=i;
    }

	for (int i=0;i<preference_cross.size();i++){ // loop on orders
		order_tree *ot=sky_view[preference_cross[i][0]];
		//cout << preference_cross[i][0].first<< ":"<< preference_cross[i][0].second <<endl;
		int j=1;
		while (j<cfg->dyDim_size){ // enter here when multi po dimension
			ot=ot->order_child[preference_cross[i][j]];
			//cout << preference_cross[i][j].first<< ":"<< preference_cross[i][j].second <<endl;
			j++;
		}

	 // 	cout << "print ids of the spo: " <<endl; 
		// for (auto tuple : ot->ids){
		// 	cout << tuple << endl;
		// }

		std::vector<id> v(skyline_result.size());   
	  	std::vector<id>::iterator it;
	  	it=std::set_difference(skyline_result.begin(), skyline_result.end(), ot->ids.begin(), ot->ids.end(), v.begin());                        
	  	v.resize(it-v.begin());   
	  	skyline_result=v; 			
	}
	// remove never_sky 
	std::vector<id> v(skyline_result.size());   
  	std::vector<id>::iterator it;
  	it=std::set_difference(skyline_result.begin(), skyline_result.end(), this->never_sky.begin(), this->never_sky.end(), v.begin());                        
  	v.resize(it-v.begin());   
  	skyline_result=v; 	

	//result.insert(result.end(), this->never_sky.begin(), this->never_sky.end());
 	
 // 	cout << "print ids of skyline result: " <<endl; 
	// for (auto tuple : skyline_result){
	// 	cout << tuple << endl;
	// }
	return skyline_result;
}

// A utility function that returns maximum of two integers 
int max(int a, int b) { return (a > b)? a : b; } 

// Returns the maximum value that can be put in a knapsack of capacity W 
int knapSack(int W, int wt[], int val[], int n) 
{ 
int i, w; 
int K[n+1][W+1]; 
bool keep[n+1][W+1];
for (i = 0; i <= n; i++) for (w = 0; w <= W; w++) keep[i][w]=0;


// Build table K[][] in bottom up manner 
for (i = 0; i <= n; i++) 
{ 
	for (w = 0; w <= W; w++) 
	{ 
		if (i==0 || w==0) 
			K[i][w] = 0; 
		else if (wt[i-1] <= w && val[i-1] + K[i-1][w-wt[i-1]] > K[i-1][w]){ 
			K[i][w] = val[i-1] + K[i-1][w-wt[i-1]];
			keep[i][w]=1; 
		}	
		else{
			K[i][w] = K[i-1][w];
			keep[i][w]=0;
		}
	} 
} 

int X=W;
for (i=n; i>=1; i--){
	if (keep[i][X]==1){
		cout << "view: " << i-1 << ",gain: "<< val[i-1]<< ",weight: "<< wt[i-1]<<endl;
		X=X-wt[i-1];
	}
}

return K[n][W]; 
} 

void comb(int N, vector<vector<int>> &all_combination)
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