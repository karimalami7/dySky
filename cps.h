/*
 * cps.h
 *
 *  Created on: January 13, 2019
 *      Author: karim
 *
 * Functionalities:  
 * 		+ Decompose a preference into chains
 * 		+ Modify the dataset
 *		+ Compute skyline by BSkyTree
 */

#include "common.h"
#include "config.h"
#include "generator/generateur.h"
#include "BSkyTree/bskytree.h"
using namespace std;

class Cps{

public: 

	vector<Point> to_dataset;
	vector<int> po_dataset;
	vector<Graph<int>> chains;
	vector<id> skyline_result;

	void decompose_preference(Graph<int> p, Config *cfg);
	int compute_skyline(Config *cfg);
};

void Cps::decompose_preference(Graph<int> p, Config *cfg){

	Graph<int> transitive_preference;

	// compute transitive

	transitive_preference.compute_transitive_closure(p);
	cout << "version transitive du graphe"<<endl;
	transitive_preference.print_edges();
	// find incomparable pairs

	vector<pair<id,id>> incomparable_pairs;
	cout << "1 / Find incomparable pairs"<<endl;
	// for all vertices i
	for (int i=0; i<transitive_preference.vertices.size(); i++){
		int value1=transitive_preference.vertices[i];
		// for all other verices j
		for (int j=i+1;j<transitive_preference.vertices.size();j++){
			int value2=transitive_preference.vertices[j];
			// if this vertex has out edges
			if (((transitive_preference.out_edges.find(value1)!=transitive_preference.out_edges.end() && 
			transitive_preference.out_edges[value1].find(value2)==transitive_preference.out_edges[value1].end())
			||transitive_preference.out_edges.find(value1)==transitive_preference.out_edges.end())
			&&((transitive_preference.out_edges.find(value2)!=transitive_preference.out_edges.end() && 
			transitive_preference.out_edges[value2].find(value1)==transitive_preference.out_edges[value2].end())
			||transitive_preference.out_edges.find(value2)==transitive_preference.out_edges.end())){
				
				incomparable_pairs.push_back(pair<id,id>(value1,value2));
				incomparable_pairs.push_back(pair<id,id>(value2,value1));
				cout <<value1<<value2<<endl;
				cout <<value2<<value1<<endl;	
			}
					
		}
	}

	// construct consistency graph
	cout <<"2 / Construct Consistency graph"<<endl;
	Graph<int> consistency_graph;
	consistency_graph.vertices=vector<int>(incomparable_pairs.size());
	for (int i=0; i<incomparable_pairs.size(); i++){
		cout <<"== Processing pair: "<< incomparable_pairs[i].first << " " << incomparable_pairs[i].second<<endl;
		Graph<int> new_preference=p;
		// add the missing order to the preference
		new_preference.out_edges[incomparable_pairs[i].first].insert(incomparable_pairs[i].second);
		// compute transitive order for the new preference
		Graph<int> new_transitive_preference;
		new_transitive_preference.compute_transitive_closure(new_preference);
		// cout << "version transitive du graphe"<<endl;
		// new_transitive_preference.print_edges();
		// find induced missing order in the transitive
		for (int j=0; j<incomparable_pairs.size(); j++){
			if(i!=j){
				auto it_source=new_transitive_preference.out_edges.find(incomparable_pairs[j].first);
				if (it_source!=new_transitive_preference.out_edges.end()){
					if (it_source->second.find(incomparable_pairs[j].second)!=it_source->second.end()){
						cout << "++++++++ induced pair: "<<incomparable_pairs[j].first <<" "<<incomparable_pairs[j].second<<endl;
						consistency_graph.add_outedges(i,{j});
					}
				}
			}
		}	
	}
	// compute non induced pairs
	vector<pair<int,int>> non_induced_pairs;
	for(int i=0; i<incomparable_pairs.size();i++){
		bool induced=false;
		for (auto it_graph=consistency_graph.out_edges.begin();it_graph!=consistency_graph.out_edges.end();it_graph++){
			if (it_graph->second.find(i)!=it_graph->second.end()){
				induced=true;break;
			}
		}
		if (!induced){
			non_induced_pairs.push_back(incomparable_pairs[i]);
		}
	}
	// print non induced pairs
	cout <<endl<<"Non induced pairs: "<<endl;
	for (auto it_vector=non_induced_pairs.begin();it_vector!=non_induced_pairs.end();it_vector++){
		cout << it_vector->first<<" "<< it_vector->second <<endl;
	}
	
	// compute incompatibility graph
	cout <<endl<<"3 / Construct Incompatibility graph"<<endl;
	Graph<int> incompatibility_graph;
	incompatibility_graph.vertices=vector<int>(non_induced_pairs.size());
	for (int i=0; i<non_induced_pairs.size(); i++){
		cout <<"== Processing pair: "<< non_induced_pairs[i].first << " " << non_induced_pairs[i].second<<endl;
		Graph<int> new_preference=p;
		// add the missing order to the preference
		new_preference.out_edges[non_induced_pairs[i].first].insert(non_induced_pairs[i].second);
		// compute transitive order for the new preference
		Graph<int> new_transitive_preference;
		new_transitive_preference.compute_transitive_closure(new_preference);
		// cout << "version transitive du graphe"<<endl;
		// new_transitive_preference.print_edges();
		for (int j=0; j<non_induced_pairs.size(); j++){
			if(i!=j){
				auto it_source=new_transitive_preference.out_edges.find(non_induced_pairs[j].second);
				if (it_source!=new_transitive_preference.out_edges.end()){
					if (it_source->second.find(non_induced_pairs[j].first)!=it_source->second.end()){
						cout << "++++++++ incompatible pair: "<<non_induced_pairs[j].first <<" "<<non_induced_pairs[j].second<<endl;
						incompatibility_graph.add_outedges(i,{j});
					}
				}
			}
		}
	}	
	// print incompatibility graph
	// cout << "incompatibility graph"<<endl;
	// incompatibility_graph.print_edges();

	// color incompatibility graph
	cout <<endl<<"4 / Color Incompatibility graph"<<endl;
	incompatibility_graph.greedyColoring();

	for(int i=0;i<incompatibility_graph.vertices.size();i++){
		cout << "pair: " <<non_induced_pairs[i].first<<" "<< non_induced_pairs[i].second<<" --> color: "
			<<incompatibility_graph.vertex_color[i]<<endl;
	}

	//create chains
	int number_colors=1;
	for (int i=0;i<incompatibility_graph.vertex_color.size();i++){
		if (number_colors<incompatibility_graph.vertex_color[i]+1){
			number_colors=incompatibility_graph.vertex_color[i]+1;
		}
	}
	cout << "number colors: "<<number_colors<<endl;

	for (int color=0;color<number_colors;color++){
		cout<< "chain for color: "<<color<<endl;
		Graph<int> pre_chain;
		pre_chain=p;
		for (int j=0;j<incompatibility_graph.vertex_color.size();j++){
			// add vertices with the same color
			if (color==incompatibility_graph.vertex_color[j]){
				pre_chain.add_outedges(non_induced_pairs[j].first,{non_induced_pairs[j].second});
			}
		}
		Graph<int> chain;
		chain.compute_transitive_closure(pre_chain);
		chain.print_edges();
		chains.push_back(chain);
	}

}

int Cps::compute_skyline(Config *cfg){
	// encoding values into <number_colors> tuple
	vector<vector<int>> values_encoding;
	for (int i=0;i<cfg->dyDim_val;i++){
		vector<int> encoding; 
		for (int j=0;j<chains.size();j++){
			auto it_graph=chains[j].out_edges.find(i);
			if (it_graph!=chains[j].out_edges.end()){
				encoding.push_back(cfg->dyDim_val-it_graph->second.size()-1);
			}
			else{
				encoding.push_back(cfg->dyDim_val-1);
			}
		}
		values_encoding.push_back(encoding);
	}

	for (int i=0;i<cfg->dyDim_val;i++){
		cout << "value: " << i << " encoding: "; 
		for (int j=0;j<chains.size();j++){
			cout << values_encoding[i][j] <<" ";
		}	
		cout<<endl;
	}

	int All = (1<<cfg->statDim_size+chains.size())-1;
	vector<Space> full_Space;
	listeAttributsPresents(All, cfg->statDim_size+chains.size(), full_Space);
	vector<Point> temp_dataset;
	for (int i=0; i< cfg->dataset_size; i++){
		Point p=(int*)malloc((cfg->statDim_size+chains.size()+1)*sizeof(int));
			for (int j=0;j<=cfg->statDim_size;j++){
				p[j]=to_dataset[i][j];
			}
			for (int k=0;k<chains.size();k++){
				p[cfg->statDim_size+k+1]=values_encoding[po_dataset[i]][k];
			}
			temp_dataset.push_back(p);
	}

    skyline_result=subspaceSkylineSize_TREE(full_Space, temp_dataset);
    cout << "Skyline size: "<<skyline_result.size()<<endl;
    return skyline_result.size();
}


