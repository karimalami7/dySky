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
	
	void views_selection(Config* cfg);
	void compute_views(Config* cfg, vector<vector<Order>> preference_orders, int *storage);
	void compute_view_1d(Config* cfg, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, int *storage);
	void compute_view_recursively_md(Config* cfg, int niveau, vector<Point> &dataset,vector<Order> orders_stack, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, int *storage);
	void compute_other_sky(Config *cfg, vector<id> &sky_for_all_orders_global, vector<Order> orders_stack);
	vector<id> compute_skyline(Config* cfg, vector<vector<Order>> preference);
	
};

int knapSack(int W, int wt[], int val[], int n); 

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

	  	//***************************************
	  	//add other candidates
	  	for (int i=0; i< this->candidates.size();i++){
	  		if( po_dataset[this->candidates[i]][0]!=best_value && po_dataset[this->candidates[i]][0]!=worst_value){
	  			sky.push_back(this->candidates[i]);
	  		}
	  	}

	  	//******************************************
	  	// sort ids
	  	std::sort (sky.begin(),sky.end());

	  	#pragma omp critical
	  	{
	  		sky_view[Order(best_value,worst_value)]=new order_tree;
	  	}
	  	sky_view[Order(best_value,worst_value)]->ids=sky;
	  	*storage=*storage+sky.size();
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
	else{
		vector<vector<Order>> orders_stack(preference_orders[0].size());
		#pragma omp parallel for schedule(dynamic)
		for(int s=0;s<preference_orders[0].size();s++){
			int best_value=preference_orders[0][s].first;
			int worst_value=preference_orders[0][s].second;
			//cout <<"niveau: "<<0<< " --> "<< best_value<<" : "<<worst_value<<endl;		
			orders_stack[s].push_back(Order(best_value,worst_value));
			
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			// partitionner les donnees par rapport à cette dimension
			vector<Point> branch_dataset;
			for (int i=0; i<candidates_tuples.size(); i++){
				if (candidates_tuples[i][cfg->statDim_size+1]==best_value){
					Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
					memcpy(p, candidates_tuples[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
					p[cfg->statDim_size+1]=0;
					branch_dataset.push_back(p);
				}else if (candidates_tuples[i][cfg->statDim_size+1]==worst_value){
					Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
					memcpy(p, candidates_tuples[i], (cfg->statDim_size+cfg->dyDim_size+1) * sizeof(int));
					p[cfg->statDim_size+1]=1;
					branch_dataset.push_back(p);
				}
			}
			#pragma omp critical
			{
				this->sky_view[Order(best_value,worst_value)]=new order_tree;;
			}		
			compute_view_recursively_md(cfg, 1, branch_dataset, orders_stack[s], this->sky_view[Order(best_value,worst_value)]->order_child, preference_orders, storage);
		}
	}

	// int views_total_storage=0;
	// for (auto it=this->sky_view.begin();it!=this->sky_view.end();it++){
	// 	views_total_storage+=(it->second)->ids.size();
	// }
	// cout << "Views total storage: " << views_total_storage <<endl;

}

void dySky::compute_view_recursively_md(Config* cfg, int niveau, vector<Point> &dataset, vector<Order> orders_stack, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, int *storage){
	
	//#pragma omp parallel for schedule(dynamic) if (niveau==0)
	for(int s=0;s<preference_orders[niveau].size();s++){
		int best_value=preference_orders[niveau][s].first;
		int worst_value=preference_orders[niveau][s].second;
		//cout <<"niveau: "<<niveau<< " --> "<< best_value<<" : "<<worst_value<<endl;		
		orders_stack.push_back(Order(best_value,worst_value));

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

			sky_view[Order(best_value,worst_value)]=new order_tree;;
			this->compute_view_recursively_md(cfg,  niveau+1, branch_dataset, orders_stack, sky_view[Order(best_value,worst_value)]->order_child, preference_orders, storage);
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
			int All = (1<<(cfg->statDim_size+cfg->dyDim_size))-1;
			vector<Space> full_Space;
			listeAttributsPresents(All, cfg->statDim_size+cfg->dyDim_size, full_Space);
			vector<id> view=subspaceSkylineSize_TREE(full_Space, branch_dataset);

			sort(view.begin(),view.end());
			
			//MINUS
			vector<id> branch_dataset_ids(branch_dataset.size());
			for (int b=0;b<branch_dataset.size();b++){
				branch_dataset_ids[b]=branch_dataset[b][0];
			}
			sort(branch_dataset_ids.begin(),branch_dataset_ids.end());
			std::vector<id> notSky(branch_dataset_ids.size());   
		  	std::vector<id>::iterator it4;
		  	it4=std::set_difference(branch_dataset_ids.begin(), branch_dataset_ids.end(), view.begin(), view.end(), notSky.begin());                        
		  	notSky.resize(it4-notSky.begin());

		  	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		  	// add points not included in the partition and remove those dominated
		  	vector<id> sky_for_all_orders_global=this->candidates;
		  	compute_other_sky(cfg, sky_for_all_orders_global, orders_stack);

			//UNION
			std::vector<id> h(this->to_dataset.size());   
		  	std::vector<id>::iterator it3;
		  	it3=std::set_union(view.begin(), view.end(), sky_for_all_orders_global.begin(), sky_for_all_orders_global.end(), h.begin());                        
		  	h.resize(it3-h.begin());   
		  	view=h; 

		  	//MINUS
			std::vector<id> g(view.size());   
		  	std::vector<id>::iterator it5;
		  	it5=std::set_difference(view.begin(), view.end(), notSky.begin(), notSky.end(), g.begin());                        
		  	g.resize(it5-g.begin());   
		  	view=g; 

		  	sky_view[Order(best_value,worst_value)]=new order_tree;
		  	sky_view[Order(best_value,worst_value)]->ids=view;
		  	*storage=*storage+view.size();
		}
		orders_stack.pop_back();
	}	
}

void dySky::compute_other_sky(Config *cfg, vector<id> &sky_for_all_orders_global, vector<Order> orders_stack){

	for(int f=0;f<cfg->dyDim_size;f++)
  	{		  		
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// partitioner par rapport aux valeurs best_value, worst_value dans la dimension "niveau"

		vector<Point> to;
		vector<vector<int>> po;
		vector<id> directly_sky;
		for (int i=0; i<this->candidates.size(); i++){
			//cout <<this->candidates[i]<<" : "<<po_dataset[this->candidates[i]][niveau]<< " : "<<best_value <<" : "<<worst_value<<endl;
			if (this->po_dataset[this->candidates[i]][f]==orders_stack[f].first){
				//cout << this->candidates[i]<<endl;
				Point p=(int*)malloc((cfg->statDim_size+2)*sizeof(int));
				vector<int> t;
				for (int j=0;j<=cfg->statDim_size;j++) p[j]=this->to_dataset[this->candidates[i]][j];
				p[cfg->statDim_size+1]=0;
				for (int j=0;j<cfg->dyDim_size;j++) if (j!=f) t.push_back(this->po_dataset[this->candidates[i]][j]);
				to.push_back(p);
				po.push_back(t);

			}else if (this->po_dataset[this->candidates[i]][f]==orders_stack[f].second){
				//cout << this->candidates[i]<<endl;
				Point p=(int*)malloc((cfg->statDim_size+2)*sizeof(int));
				vector<int> t;
				for (int j=0;j<=cfg->statDim_size;j++) p[j]=this->to_dataset[this->candidates[i]][j];
				p[cfg->statDim_size+1]=1;
				for (int j=0;j<cfg->dyDim_size;j++) if (j!=f) t.push_back(this->po_dataset[this->candidates[i]][j]);
				to.push_back(p);
				po.push_back(t);
			}
			else{
				directly_sky.push_back(this->candidates[i]);
			}
		}

		//**************************************************************************
		// partitionner par rapport aux valeurs égales sur les attributes dynamiques
		map<vector<int>, vector<Point>> partitions;
		auto it_to=to.begin();
		auto it_po=po.begin();
		while(it_to!=to.end()){
			partitions[(*it_po)].push_back(*it_to);
			it_to++;
			it_po++;
		}

		// compute skyline for every partition
		int All = (1<<cfg->statDim_size+1)-1;
		vector<Space> full_Space;
		listeAttributsPresents(All, cfg->statDim_size+1, full_Space);
		vector<id> sky_for_all_orders;
	    for (auto it_partition=partitions.begin(); it_partition!=partitions.end(); it_partition++){
	    	vector<id> Skyline;
	    	Skyline=subspaceSkylineSize_TREE(full_Space, it_partition->second);
	    	sky_for_all_orders.insert(sky_for_all_orders.end(),
	    		Skyline.begin(), Skyline.end());
	    }
	    sky_for_all_orders.insert(sky_for_all_orders.begin(), directly_sky.begin(),directly_sky.end());

	    sort(sky_for_all_orders.begin(),sky_for_all_orders.end());

	    //****************************************************************
	  	//INTERSECTION
	  	std::vector<id> l(this->to_dataset.size());   
	  	std::vector<id>::iterator it2;
	  	it2=std::set_intersection(sky_for_all_orders.begin(), sky_for_all_orders.end(), sky_for_all_orders_global.begin(), sky_for_all_orders_global.end(), l.begin());                        
	  	l.resize(it2-l.begin());   
	  	sky_for_all_orders_global=l; 

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

	vector<id> result;

	for (int i=0;i<preference_cross.size();i++){ // loop on orders
		order_tree *ot=sky_view[preference_cross[i][0]];
		int j=1;
		while (j<cfg->dyDim_size){ // enter here when multi po dimension
			ot=ot->order_child[preference_cross[i][j]];
			j++;
		}
		if(i==0){ // if there is one order
			result=ot->ids;
		}
		else{
			std::vector<id> v(result.size());   
		  	std::vector<id>::iterator it;
		  	it=std::set_intersection(result.begin(), result.end(), ot->ids.begin(), ot->ids.end(), v.begin());                        
		  	v.resize(it-v.begin());   
		  	result=v; 			
		}
  
	}

	return result;
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
