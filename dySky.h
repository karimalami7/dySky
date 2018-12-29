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
		
	vector<Point> to_dataset;
	vector<int> po_dataset;
	vector<id> always_sky;
	vector<id> never_sky;
	vector<id> candidates;
	unordered_map<Order,vector<id>, pairhash> skyline_view;
 
	public:

	void generate_to_data(Config* cfg); 
	void generate_po_data(Config* cfg);
	void compute_always_skyline(Config* cfg);
	int compute_candidates(Config* cfg);
	void compute_views(Config* cfg);
	void compute_view(Config* cfg, Order o);
	void compute_skyline(Config* cfg, vector<Order>);
	void print_dataset (Config* cfg);
	
};

void dySky::generate_to_data(Config* cfg){
	loadData("INDE","",cfg->dataset_size, cfg->statDim_size, cfg->statDim_val, this->to_dataset);

	
} 

void dySky::generate_po_data(Config* cfg){
	for (int i=0; i<cfg->dataset_size; i++) {
		this->po_dataset.push_back(rand()%cfg->dyDim_val);
	}
}

void dySky::compute_always_skyline(Config* cfg){
	
	int All = (1<<cfg->statDim_size)-1;
	vector<Space> full_Space;
	listeAttributsPresents(All, cfg->statDim_size, full_Space);
    this->always_sky=subspaceSkylineSize_TREE(full_Space, this->to_dataset);
    cout << "Always Skyline size: "<<this->always_sky.size()<<endl;
}

void dySky::print_dataset (Config* cfg){

    string const nomFichier1("INDE-"+to_string(cfg->statDim_size)+"-"+to_string(cfg->dataset_size)+".csv");
    ofstream monFlux1(nomFichier1.c_str());

    for (int i = 0; i < this->to_dataset.size() ; i++)
    {
        for (int j = 1; j <= cfg->statDim_size-1 ; j++){
            monFlux1 << this->to_dataset[i][j]<<";";  
        } 
        monFlux1 << this->to_dataset[i][cfg->statDim_size]<< endl;
    }
}

int dySky::compute_candidates(Config* cfg){

	// cluster dataset depending on the value in the po dimenion
	vector<vector<Point> > datasets(cfg->dyDim_val);
	for (int i=0; i<cfg->dataset_size; i++){
		datasets[po_dataset[i]].push_back(to_dataset[i]);
	}

	// compute skyline for every subdatasets -> always skyline + candidates
	int All = (1<<cfg->statDim_size)-1;
	vector<Space> full_Space;
	listeAttributsPresents(All, cfg->statDim_size, full_Space);
	vector<id> sky_union_candidates;
    vector<vector<id> >Skyline(cfg->dyDim_val);
    for (int i=0; i<cfg->dyDim_val; i++){
    	Skyline[i]=subspaceSkylineSize_TREE(full_Space, datasets[i]);
    	sky_union_candidates.insert(sky_union_candidates.end(),
    		Skyline[i].begin(), Skyline[i].end());
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

}

void dySky::compute_views(Config* cfg){
	for (int i=0; i<cfg->dyDim_val;i++){
		for (int j = 0; j <cfg->dyDim_val; ++j)
		{
			if (i!=j){
				//cout << "order " << i <<" "<< j<< endl; 
				compute_view(cfg,Order(i,j));
			}
		}
	}
	int views_total_storage=0;
	for (auto it=this->skyline_view.begin();it!=this->skyline_view.end();it++){
		views_total_storage+=(it->second).size();
	}
	cout << "Views total storage: " << views_total_storage <<endl;
}

void dySky::compute_view(Config* cfg, Order o){
	vector<id> sky_union_candidates;
	sky_union_candidates.insert(sky_union_candidates.end(), this->always_sky.begin(), this->always_sky.end());
	sky_union_candidates.insert(sky_union_candidates.end(), this->candidates.begin(), this->candidates.end());
	vector<Point> temp_dataset;
	
	for (int i=0; i<sky_union_candidates.size(); i++){
		if (this->po_dataset[sky_union_candidates[i]]==o.first){
			Point p=(int*)malloc((cfg->statDim_size+2)*sizeof(int));
			for (int j=0;j<=cfg->statDim_size;j++){
				p[j]=to_dataset[sky_union_candidates[i]][j];
			}
			p[cfg->statDim_size+1]=0;
			temp_dataset.push_back(p);
		}else if (this->po_dataset[sky_union_candidates[i]]==o.second){
			Point p=(int*)malloc((cfg->statDim_size+2)*sizeof(int));
			for (int j=0;j<=cfg->statDim_size;j++){
				p[j]=to_dataset[sky_union_candidates[i]][j];
			}
			p[cfg->statDim_size+1]=1;
			temp_dataset.push_back(p);
		}
	}

	// compute skyline set for this view
	
	int All = (1<<(cfg->statDim_size+1))-1;
	vector<Space> full_Space;
	listeAttributsPresents(All, cfg->statDim_size+1, full_Space);
	vector<id> view_sky=subspaceSkylineSize_TREE(full_Space, temp_dataset);
	std::vector<int> view_candidates(view_sky.size());  //filled with zeros
  	std::vector<int>::iterator it_view;
  	std::sort (view_sky.begin(),view_sky.end());
  	it_view=std::set_difference (view_sky.begin(),view_sky.end(), 
  		this->always_sky.begin(),this->always_sky.end(), view_candidates.begin()); // diff between all ids and dom objects
  	view_candidates.resize(it_view-view_candidates.begin());  
    this->skyline_view.insert(pair<Order, vector<id> >(o,view_candidates));

}



