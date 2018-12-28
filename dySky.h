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
	map<std::set<Order>, vector<int> > view;
	vector<id> dominated_objects;

	public:

	void generate_to_data(Config* cfg); 
	void generate_po_data(Config* cfg);
	void compute_skyline(Config* cfg);
	void compute_view(set<Order>);
	void print_dataset (Config* cfg);
	int compute_dominated_objects(Config* cfg);
};

void dySky::generate_to_data(Config* cfg){
	loadData("INDE","",cfg->dataset_size, cfg->statDim_size, cfg->statDim_val, this->to_dataset);
} 

void dySky::generate_po_data(Config* cfg){
	for (int i=0; i<cfg->dataset_size; i++) {
		this->po_dataset.push_back(rand()%cfg->dyDim_val);
	}
}

void dySky::compute_skyline(Config* cfg){
	
	int All = (1<<cfg->statDim_size)-1;
	vector<Space> full_Space;
	listeAttributsPresents(All, cfg->statDim_size, full_Space);
    std::vector<id> Skyline;
    Skyline=subspaceSkylineSize_TREE(full_Space, this->to_dataset);
    cout << "minimum Skyline size: "<<Skyline.size()<<endl;
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

int dySky::compute_dominated_objects(Config* cfg){

	// cluster dataset depending on the value in the po dimenion
	vector<Point> datasets[cfg->dyDim_val];
	for (int i=0; i<cfg->dataset_size; i++){
		datasets[po_dataset[i]].push_back(to_dataset[i]);
	}

	// compute skyline for every subdatasets
	int All = (1<<cfg->statDim_size)-1;
	vector<Space> full_Space;
	listeAttributsPresents(All, cfg->statDim_size, full_Space);
    vector<id> Skyline[cfg->dyDim_val];
    for (int i=0; i<cfg->dyDim_val; i++){
    	Skyline[i]=subspaceSkylineSize_TREE(full_Space, datasets[i]);
    	this->dominated_objects.insert(this->dominated_objects.end(),
    		Skyline[i].begin(), Skyline[i].end());// it is actually skyline objects here
    }
    vector<id> all_ids;
    for (int i=0; i<cfg->dataset_size; i++){
    	all_ids.push_back(i);
    }
    //all_ids minus dominated_objects
    std::vector<int> minus_vector(cfg->dataset_size);  //filled with zeros
  	std::vector<int>::iterator it;
      
  	std::sort (this->dominated_objects.begin(),this->dominated_objects.end());   
  	it=std::set_difference (all_ids.begin(),all_ids.end(), 
  		this->dominated_objects.begin(),this->dominated_objects.end(), minus_vector.begin()); // diff between all ids and dom objects
  	minus_vector.resize(it-minus_vector.begin());  
    this->dominated_objects.swap(minus_vector);
    cout << "There is "<< dominated_objects.size()<<" dominated objects"<<endl;

 	

	int count_test=0;
	for (int i=0; i<cfg->dyDim_val; i++){
		cout << "size of Skyline of cluster " << i << ": "<<Skyline[i].size()<<endl;
		count_test+=Skyline[i].size();
	}
	cout <<endl;
	cout <<"sum (maximum Skyline size): "<<count_test<<endl;
}