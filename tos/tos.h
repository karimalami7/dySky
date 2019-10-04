/*
 * tos.h
 *
 *  Created on: Feb 2, 2019
 *      Author: karim
 *
 * Functionalities:  
 * 		+ 
 */

#include "../common/common.h"
#include "../common/config.h"
#include "../generator/generateur.h"
#include "../BSkyTree/bskytree.h"
using namespace std;

class Tos{

public:

	vector<Point> to_dataset;
	vector<vector<int>> po_dataset;
	map<chain, chain_tree*> sky_view;
	vector<vector<chain>> chains_vec_cross;
	vector<vector<chain>> paths;

	Tos(Config *cfg);
	void compute_views(Config *cfg, uint64_t *storage);
	void compute_view_recursively(Config *cfg, int niveau, vector<Preference> preference_stack, map<chain, chain_tree*> &view_node, uint64_t *storage);
	void compute_view_1d(Config *cfg, map<chain, chain_tree*> &view_node, uint64_t *storage);
	vector<id> compute_skyline(Config *cfg);
	void define_paths(vector<Preference> p, Config *cfg);
	void chain_graph_to_vec_representation(vector<vector<Graph<int>>> chains, Config *cfg);
	void chain_cross(Config* cfg, vector<chain> v, vector<vector<chain>> chains_vec, int niv);
};

Tos::Tos(Config *cfg){
	this->paths=vector<vector<chain>>(cfg->dyDim_size);
}

void Tos::compute_view_1d(Config *cfg, map<chain, chain_tree*> &view_node, uint64_t *storage){
	
	chain tos_values(cfg->dyDim_val);
	for (int i=0; i<cfg->dyDim_val; i++){tos_values[i]=i;}
	vector<chain> all_permutation;
	do{
		all_permutation.push_back(tos_values);
	}
	while ( std::next_permutation(tos_values.begin(),tos_values.begin()+cfg->dyDim_val) );

	#pragma omp parallel for schedule(dynamic)
	for (int i=0;i<all_permutation.size();i++){
  	  	Preference p_to;
		p_to.add_vertices(cfg->dyDim_val);
		for (id source=0;source<cfg->dyDim_val-1;source++){
			unordered_set<id> v_dest;
			v_dest.insert(all_permutation[i][source+1]);
			p_to.add_edges(all_permutation[i][source],v_dest);
		}
		Cps *cps_tos = new Cps(cfg);
		cps_tos->decompose_preference(p_to,cfg,0);	
		cps_tos->encoding(cfg);
		cps_tos->compute_skyline(cfg, false);
		#pragma omp critical
		{
			view_node[all_permutation[i]]=new chain_tree;
		}
		view_node[all_permutation[i]]->ids.swap(cps_tos->skyline_result);
		*storage=*storage+view_node[all_permutation[i]]->ids.size();
		delete cps_tos;
	}
}

void Tos::compute_views(Config *cfg, uint64_t *storage){
	//cout << "Tos::compute_views"<<endl;
  	
  	if(cfg->dyDim_size==1){
  		this->compute_view_1d(cfg, this->sky_view, storage);
  	}
  	else{
		chain tos_values(cfg->dyDim_val);
		for (int i=0; i<cfg->dyDim_val; i++){tos_values[i]=i;}
		vector<chain> all_permutation;
		do{
			all_permutation.push_back(tos_values);
		}
		while ( std::next_permutation(tos_values.begin(),tos_values.begin()+cfg->dyDim_val) );


		//*************************************

		vector<vector<Preference>> preference_stack(all_permutation.size());
		
		#pragma omp parallel for schedule(dynamic)
		for (int j=0;j<all_permutation.size();j++) {
			//************************************
			// memory
			struct sysinfo sys_info;
			uint64_t info_ram;
			if (!(sysinfo(&sys_info) == -1)) {
				info_ram=sys_info.totalram - sys_info.freeram;
				info_ram = (info_ram * sys_info.mem_unit)/1024;
				printf("memory %d\n", info_ram);
			}
			//************************************

	  	  	Preference p_to;
			p_to.add_vertices(cfg->dyDim_val);
			for (id source=0;source<cfg->dyDim_val-1;source++){
				unordered_set<id> v_dest;
				v_dest.insert(all_permutation[j][source+1]);
				p_to.add_edges(all_permutation[j][source],v_dest);
			}
			preference_stack[j].push_back(p_to);	
			#pragma omp critical
			{	
				this->sky_view[all_permutation[j]]=new chain_tree;
			}		
  			this->compute_view_recursively(cfg, 1, preference_stack[j], this->sky_view[all_permutation[j]]->chain_child, storage);

  		}
  	}
  	
}

void Tos::compute_view_recursively(Config *cfg, int niveau, vector<Preference> preference_stack, map<chain, chain_tree*> &view_node, uint64_t *storage){
	//cout << "Tos::compute_view_recursively, niveau "<< niveau<<endl;

	chain tos_values(cfg->dyDim_val);
	for (int i=0; i<cfg->dyDim_val; i++){tos_values[i]=i;}
	vector<chain> all_permutation;
	do{
		all_permutation.push_back(tos_values);
	}
	while ( std::next_permutation(tos_values.begin(),tos_values.begin()+cfg->dyDim_val) );

	//************************************
	// memory
	// struct sysinfo sys_info;
	// uint64_t info_ram;
	// if (!(sysinfo(&sys_info) == -1)) {
	// 	info_ram=sys_info.totalram - sys_info.freeram;
	// 	info_ram = (info_ram * sys_info.mem_unit)/1024;
	// 	printf("memory %d\n", info_ram);
	// }

	for (int j=0;j<all_permutation.size();j++) {
  	  	Preference p_to;
		p_to.add_vertices(cfg->dyDim_val);
		for (id source=0;source<cfg->dyDim_val-1;source++){
			unordered_set<id> v_dest;
			v_dest.insert(all_permutation[j][source+1]);
			p_to.add_edges(all_permutation[j][source],v_dest);
		}
		preference_stack.push_back(p_to);
		if (niveau<cfg->dyDim_size -1){
			view_node[all_permutation[j]]=new chain_tree;
			this->compute_view_recursively(cfg,  niveau+1, preference_stack, view_node[all_permutation[j]]->chain_child, storage);
		}
		else if (niveau==cfg->dyDim_size-1){
			Cps *cps_tos = new Cps(cfg);
			for (int i=0;i<cfg->dyDim_size;i++){
				cps_tos->decompose_preference(preference_stack[i],cfg,i);	
			}
			cps_tos->encoding(cfg);
			cps_tos->compute_skyline(cfg, false);
			view_node[all_permutation[j]]=new chain_tree;
		  	view_node[all_permutation[j]]->ids.swap(cps_tos->skyline_result);
		  	*storage=*storage+view_node[all_permutation[j]]->ids.size();
			delete cps_tos;
		}	
		preference_stack.pop_back();
  	} 
}

void Tos::define_paths(vector<Preference> p, Config *cfg){


	for (int i=0; i<p.size();i++){
		vector<Order> all_orders;
		p[i].print_edges();
		for (auto it=p[i].out_edges.begin(); it!=p[i].out_edges.end();it++){
			for (auto it2=it->second.begin();it2!=it->second.end();it2++){
				all_orders.push_back(Order(it->first,*it2));
			}
		}
	

		// generer toutes les permutations qui contiennent ces orders
		vector<int> tos_values(cfg->dyDim_val);
	  	for (int l=0; l<cfg->dyDim_val; l++){tos_values[l]=l;}

	  	do{
	  		bool path_to_include=true;
  			// for (auto p : tos_values){
  			// 	cout << p << " : ";
  			// }
  			// cout <<endl;

	  		for(auto a=all_orders.begin(); a!=all_orders.end() &&  path_to_include ;a++){
	  			// check if this order is included in path
	  			for (int j=0;j<cfg->dyDim_val && path_to_include ;j++){
	  				// cout << "j "<<j<<endl;
	  				// cout << "j "<<j<< " : "<< tos_values[j]<<" : " << a->first<<endl;
	  				if (tos_values[j]==a->first){
	  					// cout << "j "<<j<< " : "<< tos_values[j]<<" : " << a->first<<endl;
	  					for (int k=0;k<j && path_to_include;k++){
	  						// cout << "k "<<k<<endl;
	  						if(tos_values[k]==a->second) {
	  							// cout << "k "<<k<< " : "<< tos_values[j]<<" : " << a->second<<endl;
	  							// cout << "hs"<<endl;
	  							path_to_include=false;
	  						}
	  					}

	  				}
	  			}
	  		}
	  		// cout <<"ici"<<endl;
	  		if(path_to_include){
	  			this->paths[i].push_back(tos_values);
	  			// for (auto p : tos_values){
	  			// 	cout << p << " : ";
	  			// }
	  			// cout <<endl;
	  		}

	  	}while ( std::next_permutation(tos_values.begin(),tos_values.begin()+cfg->dyDim_val) );

	}
}

vector<id> Tos::compute_skyline(Config *cfg){

	// cout << "Tos::compute_skyline"<<endl;
	
	// this->chain_graph_to_vec_representation(chains, cfg);

	vector<chain> v;
	this->chains_vec_cross.clear();
	Tos::chain_cross(cfg, v, this->paths, 0);

	// cout << "chains_vec_cross size: "<<this->chains_vec_cross.size()<<endl;
	// for (int i=0;i<this->chains_vec_cross.size();i++){
	// 	cout << "croisement " <<i<<endl;
	// 	for (int j=0;j<this->chains_vec_cross[i].size();j++){
	// 		for(auto value: this->chains_vec_cross[i][j]){
	// 			cout << value <<" ";
	// 		}
	// 		cout<<endl;	
	// 	}
	// 	cout<<endl;
	// }
	//***************************************************
	// union the skyline of the chains
	vector<id> result;

	for (int i=0;i<this->chains_vec_cross.size();i++){
		chain_tree *ct=this->sky_view[this->chains_vec_cross[i][0]];
		int j=1;
		while (j<cfg->dyDim_size){

			ct=ct->chain_child[this->chains_vec_cross[i][j]];
			j++;
		}
		// cout <<"ids"<<endl;
		// for (int i=0; i<ct->ids.size();i++){
		// 	cout << ct->ids[i] <<endl;
		// }
		//cout <<endl;

		if(i==0){
			result=ct->ids;

		}
		else{	
			std::vector<id> v(result.size()+ct->ids.size());   
		  	std::vector<id>::iterator it;
		  	it=std::set_union(result.begin(), result.end(), ct->ids.begin(), ct->ids.end(), v.begin());                        
		  	v.resize(it-v.begin());   
		  	result=v; 			
		}
		//cout <<"result"<<endl;
		for (int i=0; i<result.size();i++){
			//cout << result[i] <<endl;
		}
	}

	// cout << "print ids at the end" <<endl; 
	// for (auto tuple : result){
	// 	cout << tuple << endl;
	// }
	
	return result;
}









void Tos::chain_graph_to_vec_representation(vector<vector<Graph<int>>> chains, Config *cfg){
		
	//cout << "Tos::chain_graph_to_vec_representation"<<endl;

	// for (int i=0;i<chains.size();i++){
	// 	for (int j=0;j<chains[i].size();j++){
	// 		for(auto value: chains[i][j]){
	// 			cout << value <<" ";
	// 		}
	// 		cout<<endl;	
	// 	}
	// 	cout<<endl;
	// }


	//**************************************************
	// transforms chains from graph representation to vector representation
	vector<vector<chain>> chains_vec(cfg->dyDim_size);
	for (int i=0;i<cfg->dyDim_size;i++){
		//cout << "Dimension: "<<i<<endl;
		for (Graph<int> g : chains[i]){
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
			//for (int s : c) cout << s << " ";
			//cout <<endl;
			chains_vec[i].push_back(c);
		}		
	}
	//***************************************************

	// cout << "chains_vec size: " <<chains_vec.size()<<endl;
	// for (auto ch : chains_vec){
	// 	cout << ch.size() << " : " ;
	// }
	// cout <<endl;

	vector<chain> v;
	this->chains_vec_cross.clear();
	Tos::chain_cross(cfg, v, chains_vec, 0);

}


void Tos::chain_cross(Config* cfg, vector<chain> v, vector<vector<chain>> chains_vec, int niv){

	for (int i=0; i<chains_vec[niv].size(); i++){
		v.push_back(chains_vec[niv][i]);
		if (niv<cfg->dyDim_size-1){
			chain_cross(cfg, v, chains_vec, niv+1);
			}
		if (niv==cfg->dyDim_size-1){

			this->chains_vec_cross.push_back(v);
		}
		v.pop_back();
	}
}