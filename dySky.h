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
	map<std::set<Order>, vector<int>> view;

	public:

	void generate_to_data(Config* cfg); 
	void generate_po_data(Config* cfg);
	void compute_skyline(Config* cfg);
	void compute_view(set<Order>);
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

void dySky::compute_skyline(Config* cfg){
	
	int All = (1<<cfg->statDim_size)-1;
	vector<Space> full_Space;
	listeAttributsPresents(All, cfg->statDim_size, full_Space);
    std::vector<Point> Skyline;
    Skyline=subspaceSkylineSize_TREE(full_Space, this->to_dataset);
    cout << "Skyline size: "<<Skyline.size()<<endl;
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