/*
 * tos.h
 *
 *  Created on: Feb 2, 2019
 *      Author: karim
 *
 * Functionalities:  
 * 		+ 
 */

#include "common.h"
#include "config.h"
#include "generator/generateur.h"
#include "BSkyTree/bskytree.h"
using namespace std;

class Tos{

public:

	vector<Point> to_dataset;
	vector<int> po_dataset;
	map<chain, vector<id>> cache;

	vector<id> compute_skyline(vector<Graph<int>> chains, Config *cfg);

	void chain_graph_to_vec_representation(vector<Graph<int>> chains, Config *cfg, vector<chain> &chains_vec);

};

vector<id> Tos::compute_skyline(vector<Graph<int>> chains, Config *cfg){

	cout << "Tos::compute_skyline"<<endl;
	
	vector<chain> chains_vec;
	
	this->chain_graph_to_vec_representation(chains, cfg, chains_vec);

	//***************************************************
	// union the skyline of the chains
	vector<id> result;
	if (chains_vec.size()==1){
		result=this->cache[chains_vec[0]];
	}
	else{

		result=this->cache[chains_vec[0]];
		for (int i=1; i<chains_vec.size(); i++){
			sort(result.begin(),result.end());
			sort(cache[chains_vec[i]].begin(), cache[chains_vec[i]].end());
			std::vector<id> v(result.size()+cache[chains_vec[i]].size());   
  			std::vector<id>::iterator it;
  			it=std::set_union(result.begin(), result.end(), this->cache[chains_vec[i]].begin(),
  				this->cache[chains_vec[i]].end(), v.begin());                        
  			v.resize(it-v.begin());   
  			result=v;         
		}
	}
	return result;

}









void Tos::chain_graph_to_vec_representation(vector<Graph<int>> chains, Config *cfg, vector<chain> &chains_vec){
		
	cout << "Tos::chain_graph_to_vec_representation"<<endl;

	//**************************************************
	// transforms chains from graph representation to vector representation
	for (Graph<int> g : chains){
		g.print_edges();
		chain c;
		int num_childs=cfg->dyDim_val;
		vector<int> vertices=g.vertices;
		vector<int> vertices_added;
		while(num_childs>0){
			for (auto it : g.out_edges){
				if (it.second.size()==num_childs){
					c.push_back(it.first);
					vertices_added.push_back(it.first);
				}
			}
			num_childs--;
		}	
		sort(vertices.begin(), vertices.end());
		sort(vertices_added.begin(), vertices_added.end());
		std::vector<int> v(vertices.size()); 
		std::vector<int>::iterator it;
		it=std::set_difference (vertices.begin(), vertices.end(), vertices_added.begin(), vertices_added.end(), v.begin());
		v.resize(it-v.begin());
		c.push_back(v[0]);
		for (int s : c) cout << s << " ";
		cout <<endl;
		chains_vec.push_back(c);
	}
	//***************************************************
}