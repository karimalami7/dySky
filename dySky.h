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
	vector<id> compute_skyline(Config* cfg, vector<Order> preference);
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

    for (int i = 0; i < this->to_dataset.size() ; i++)
    {
        for (int j = 1; j <= cfg->statDim_size ; j++){
            monFlux1 << this->to_dataset[i][j]<<";";  
        } 
        monFlux1 << this->po_dataset[i]<< endl;
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

 	string const nomFichier1("always_sky_2.csv");
    ofstream monFlux1(nomFichier1.c_str());

    for (int i = 0; i < this->always_sky.size() ; i++)
    {
        for (int j = 0; j <= cfg->statDim_size ; j++){
            monFlux1 << this->to_dataset[this->always_sky[i]][j]<<";";  
        } 
        monFlux1 << this->po_dataset[this->always_sky[i]]<< endl;
    }

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
	// cerr << "sky_union_candidates"<<endl;
	// for (auto it = sky_union_candidates.begin();it!=sky_union_candidates.end(); it++){	
	// 	cerr << (*it)<<endl;
	// }
	vector<Point> temp_dataset;
	
	for (int i=0; i<sky_union_candidates.size(); i++){
		if (this->po_dataset[sky_union_candidates[i]]==o.first){
			//cerr<<sky_union_candidates[i]<<"first"<<endl;
			Point p=(int*)malloc((cfg->statDim_size+2)*sizeof(int));
			for (int j=0;j<=cfg->statDim_size;j++){
				p[j]=to_dataset[sky_union_candidates[i]][j];
			}
			p[cfg->statDim_size+1]=0;
			temp_dataset.push_back(p);
		}else if (this->po_dataset[sky_union_candidates[i]]==o.second){
			//cerr<<sky_union_candidates[i]<<"second"<<endl;
			Point p=(int*)malloc((cfg->statDim_size+2)*sizeof(int));
			for (int j=0;j<=cfg->statDim_size;j++){
				p[j]=to_dataset[sky_union_candidates[i]][j];
			}
			p[cfg->statDim_size+1]=1;
			temp_dataset.push_back(p);
		}
	}
	//cerr << "temp_dataset"<<endl;
	//for (int i=0 ; i< temp_dataset.size(); i++) cerr <<temp_dataset[i][cfg->statDim_size+1]<<endl;
	// compute skyline set for this view
	
	int All = (1<<(cfg->statDim_size+1))-1;
	vector<Space> full_Space;
	listeAttributsPresents(All, cfg->statDim_size+1, full_Space);
	vector<id> view_sky=subspaceSkylineSize_TREE(full_Space, temp_dataset);
	//remove always skyline ids 
	std::vector<id> view_candidates(view_sky.size());  
  	std::vector<id>::iterator it_view;
  	std::sort (view_sky.begin(),view_sky.end());
  	it_view=std::set_difference (view_sky.begin(),view_sky.end(), 
  		this->always_sky.begin(),this->always_sky.end(), view_candidates.begin()); 
  	view_candidates.resize(it_view-view_candidates.begin());  
  	//add other candidates
  	for (auto it=this->candidates.begin(); it!=this->candidates.end(); it++){
  		if (this->po_dataset[(*it)]!=o.first && this->po_dataset[(*it)]!=o.second){
  			view_candidates.push_back(*it);
  		}
  	}
  	std::sort (view_candidates.begin(),view_candidates.end());
    this->skyline_view.insert(pair<Order, vector<id> >(o,view_candidates));
}

vector<id> dySky::compute_skyline(Config* cfg, vector<Order> preference){
	// compute union of all concerned views
	vector<id> result;
	if (preference.size()==1){
		result=this->skyline_view[preference[0]];
	}
	else{

		result=this->skyline_view[preference[0]];
		for (int i=1; i<preference.size(); i++){
			std::vector<id> v(result.size());   
  			std::vector<id>::iterator it;
  			it=std::set_intersection(result.begin(), result.end(), this->skyline_view[preference[i]].begin(),
  				this->skyline_view[preference[i]].end(), v.begin());                        
  			v.resize(it-v.begin());   
  			result=v;         
		}
	}
	return result;
}
