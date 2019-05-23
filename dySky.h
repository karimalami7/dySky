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
	vector<id> always_sky;
	vector<id> never_sky;
	vector<id> candidates;
	unordered_map<Order, order_tree*, pairhash> sky_view;

	void generate_to_data(Config* cfg); 
	void generate_po_data(Config* cfg);
	void compute_always_skyline(Config* cfg);
	int compute_candidates(Config* cfg);
	
	void print_dataset (Config* cfg);
	void print_dataset (vector<Point> data, string name, int d);
	dySky(Config *cfg);
	
	void compute_views(Config* cfg, vector<vector<Order>> preference_orders);
	void compute_view_recursively(Config* cfg, int niveau, vector<Point> &dataset, vector<Order> orders, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders);
	vector<id> compute_skyline(Config* cfg, vector<vector<Order>> preference);

};

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

void dySky::compute_always_skyline(Config* cfg){
	cout << "dySky::compute_always_skyline"<<endl;

	int All = (1<<cfg->statDim_size)-1;
	vector<Space> full_Space;
	listeAttributsPresents(All, cfg->statDim_size, full_Space);
    this->always_sky=subspaceSkylineSize_TREE(full_Space, to_dataset);
    cout << "Always Skyline size: "<<this->always_sky.size()<<endl;
    // cerr << "::" <<endl;
    // for (int k=0; k< this->always_sky.size(); k++) cerr << this->always_sky[k]<<endl; 
	// int All = (1<<cfg->statDim_size+1)-1;
	// vector<Space> full_Space;
	// listeAttributsPresents(All, cfg->statDim_size+1, full_Space);
	// vector<Point> temp_dataset;
	// for (int i=0; i< cfg->dataset_size; i++){
	// 	Point p=(int*)malloc((cfg->statDim_size+2)*sizeof(int));
	// 		for (int j=0;j<=cfg->statDim_size;j++){
	// 			p[j]=to_dataset[i][j];
	// 		}
	// 		p[cfg->statDim_size+1]=po_dataset[i];
	// 		temp_dataset.push_back(p);
	// }
	// // gooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooood      results
 //    this->always_sky=subspaceSkylineSize_TREE(full_Space, temp_dataset);
 //    cout << "Always Skyline size: "<<this->always_sky.size()<<endl;
 //    cerr << "::" <<endl;
 //    for (int k=0; k< this->always_sky.size(); k++) cerr << this->always_sky[k]<<endl; 

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
	vector<id> sky_union_candidates;
    for (auto it_partition=partitions.begin(); it_partition!=partitions.end(); it_partition++){
    	vector<id> Skyline;
    	Skyline=subspaceSkylineSize_TREE(full_Space, it_partition->second);
    	sky_union_candidates.insert(sky_union_candidates.end(),
    		Skyline.begin(), Skyline.end());
    }

    // compute never_sky and candidates
    vector<id> all_ids;
    for (int i=0; i<cfg->dataset_size; i++){
    	all_ids.push_back(i);
    }

    // never_sky = all_ids minus never_sky
    std::vector<int> never_sky(cfg->dataset_size);  //filled with zeros
  	std::vector<int>::iterator it_never;
  	std::sort (sky_union_candidates.begin(),sky_union_candidates.end());   
  	it_never=std::set_difference (all_ids.begin(),all_ids.end(), 
  		sky_union_candidates.begin(),sky_union_candidates.end(), never_sky.begin()); // diff between all ids and dom objects
  	never_sky.resize(it_never-never_sky.begin());  
    this->never_sky.swap(never_sky);
    cout <<"Never Skyline size: "<< this->never_sky.size()<<endl;

    //  candidates = sky_union_candidates minus sky
    std::vector<int> candidates(cfg->dataset_size);  //filled with zeros
  	std::vector<int>::iterator it_candidates;
  	std::sort (this->always_sky.begin(),this->always_sky.end()); 
  	it_candidates=std::set_difference (sky_union_candidates.begin(),sky_union_candidates.end(), 
  		this->always_sky.begin(),this->always_sky.end(), candidates.begin()); // diff between all ids and dom objects
  	candidates.resize(it_candidates-candidates.begin());  
    this->candidates.swap(candidates);
    cout << "Candidate set size: "<<this->candidates.size()<<endl;

 	cout << "Totally, there is "<<this->candidates.size() + this->never_sky.size() + this->always_sky.size()<<" objects"<<endl; 

 	// string const nomFichier1("always_sky_2.csv");
  //   ofstream monFlux1(nomFichier1.c_str());

  //   for (int i = 0; i < this->always_sky.size() ; i++)
  //   {
  //       for (int j = 0; j <= cfg->statDim_size ; j++){
  //           monFlux1 << this->to_dataset[this->always_sky[i]][j]<<";";  
  //       } 
  //       monFlux1 << this->po_dataset[this->always_sky[i]][0]<< endl;
  //   }

}

void dySky::compute_views(Config* cfg, vector<vector<Order>> preference_orders){
	cout << "dySky::compute_views" <<endl;

	vector<id> aSky_union_candidates;
	aSky_union_candidates.insert(aSky_union_candidates.end(), this->always_sky.begin(), this->always_sky.end());
	aSky_union_candidates.insert(aSky_union_candidates.end(), this->candidates.begin(), this->candidates.end());
	
	vector<Point> aSky_union_candidates_tuples;
	for (int i=0; i<aSky_union_candidates.size(); i++){
		Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
		for (int j=0;j<=cfg->statDim_size;j++){
			p[j]=to_dataset[aSky_union_candidates[i]][j];
		}
		for (int j=0;j<cfg->dyDim_size;j++){
			p[j+(cfg->statDim_size+1)]=po_dataset[aSky_union_candidates[i]][j];
		}
		//p[cfg->statDim_size+1]=po_dataset[aSky_union_candidates[i]][0];
		aSky_union_candidates_tuples.push_back(p);
	}
	
	print_dataset(aSky_union_candidates_tuples, " compute_view2 data", cfg->statDim_size+cfg->dyDim_size);

	vector<Order> orders;
	cout << "-->dySky::compute_view_recursively"<<endl;
	compute_view_recursively(cfg, 0, aSky_union_candidates_tuples, orders, this->sky_view, preference_orders);

	// int views_total_storage=0;
	// for (auto it=this->sky_view.begin();it!=this->sky_view.end();it++){
	// 	views_total_storage+=(it->second)->ids.size();
	// }
	// cout << "Views total storage: " << views_total_storage <<endl;

}

void dySky::compute_view_recursively(Config* cfg, int niveau, vector<Point> &dataset, vector<Order> orders_stack, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders){
	
	
	for(int i=0;i<preference_orders[niveau].size();i++){
		int best_value=preference_orders[niveau][i].first;
		int worst_value=preference_orders[niveau][i].second;
				
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
		orders_stack.push_back(Order(best_value,worst_value));
			
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//si c'est pas la dernière dimension 				
			if (niveau<cfg->dyDim_size -1){
				sky_view[Order(best_value,worst_value)]=new order_tree;;
				dySky::compute_view_recursively(cfg,  niveau+1, branch_dataset, orders_stack, sky_view[Order(best_value,worst_value)]->order_child, preference_orders);
			}

			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//si c'est la dernière dimension
			if (niveau==cfg->dyDim_size-1){

				//*************************************
				// calculer le skyline
				int All = (1<<(cfg->statDim_size+cfg->dyDim_size))-1;
			vector<Space> full_Space;
			listeAttributsPresents(All, cfg->statDim_size+cfg->dyDim_size, full_Space);
			vector<id> view_sky=subspaceSkylineSize_TREE(full_Space, branch_dataset);

			//**************************************
			//remove always skyline ids 
			std::vector<id> view_candidates(view_sky.size());  
		  	std::vector<id>::iterator it_view;
		  	std::sort (view_sky.begin(),view_sky.end());
		  	it_view=std::set_difference (view_sky.begin(),view_sky.end(), 
		  		this->always_sky.begin(),this->always_sky.end(), view_candidates.begin()); 
		  	view_candidates.resize(it_view-view_candidates.begin());  

		  	//***************************************
		  	//add other candidates
		  	for (int i=0; i< this->candidates.size();i++){
		  		bool to_add=false;
		  		for(int d=0;d<cfg->dyDim_size;d++){
		  			if( po_dataset[candidates[i]][d]!=orders_stack[d].first && po_dataset[candidates[i]][d]!=orders_stack[d].second){
		  			to_add=true;
		  			}
		  		}
		  		if(to_add){
		  			view_candidates.push_back(candidates[i]);
		  		}
		  	}

		  	//******************************************
		  	// sort ids
		  	std::sort (view_candidates.begin(),view_candidates.end());
		  	sky_view[Order(best_value,worst_value)]=new order_tree;
		  	sky_view[Order(best_value,worst_value)]->ids=view_candidates;

		}
		orders_stack.pop_back();
	}	
}



vector<id> dySky::compute_skyline(Config* cfg, vector<vector<Order>> preference){
	cout << "dySky::compute_skyline" <<endl;
	// Sky(T) = Always_sky Union ( Intersec(Views concerned ) ) 
	
	for (int i=0;i<preference.size();i++){
		for (int j=0;j<preference[i].size();j++){
			cout << preference[i][j].first<< ":"<< preference[i][j].second <<endl;
		}
		cout<<endl;
	}

	vector<id> result;

	for (int i=0;i<preference.size();i++){
		order_tree *ot=sky_view[preference[i][0]];
		int j=1;
		while (j<cfg->dyDim_size){
			ot=ot->order_child[preference[i][j]];
			j++;
		}
		if(i==0){
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

	result.insert(result.end(), always_sky.begin(), always_sky.end());


	return result;
}
