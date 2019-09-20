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

#include "../common/common.h"
#include "../common/config.h"
#include "../generator/generateur.h"
#include "../BSkyTree/bskytree.h"
#include <omp.h>
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
	void views_selection(Config* cfg, uint64_t max_storage);
	void compute_views(Config* cfg, vector<vector<Order>> preference_orders, uint64_t *storage);
	void compute_view_1d(Config* cfg, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, uint64_t *storage);
	void compute_view_recursively_md(Config* cfg, int niveau, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, uint64_t *storage);
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




void dySky::compute_view_1d(Config* cfg, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, uint64_t *storage){
	
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

void dySky::compute_views(Config* cfg, vector<vector<Order>> preference_orders, uint64_t *storage){
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

		// this->partition_ids=vector<map<int, vector<int>>>(cfg->dyDim_size);
		// for (int i=0; i<this->candidates.size(); i++){
		// 	for (int j=0; j<cfg->dyDim_size; j++){
		// 		partition_ids[j][this->po_dataset[this->candidates[i]][j]].push_back(this->candidates[i]);
		// 	}	
		// }

		preference_orders[0].push_back(Order(-1,-1));
		for(int s=0;s<preference_orders[0].size();s++){
			this->sky_view[Order(preference_orders[0][s].first,preference_orders[0][s].second)]=new order_tree;
		}
		int niveau=0;
		#pragma omp parallel for schedule(dynamic)
		for(int s=0;s<preference_orders[0].size();s++){
			int best_value=preference_orders[0][s].first;
			int worst_value=preference_orders[0][s].second;
			//cout <<"niveau: "<<niveau<< " --> "<< best_value<<" : "<<worst_value<<endl;		
			
			if (best_value!=-1){
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
				compute_view_recursively_md(cfg, niveau+1, branch_dataset, this->sky_view[Order(best_value,worst_value)]->order_child, preference_orders, storage);
				
				/////
				//destroy Point pointers
				for ( auto p : branch_dataset)
				delete p;				
			}
			else{
				//here, we parition by values in dimension 0 and we call the recursive function Order (-1,-1)
		 		
				//cout <<"niveau: "<<niveau<< " --> "<< "-1"<<" : "<<"-1"<<endl;
				//auto start_time=omp_get_wtime();
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

				/////
				//destroy Point pointers
				for ( auto partition : partitions)
				for ( auto p : partition.second)
				delete p;
				//	
				//cout <<"time for (-1,-1): "<<omp_get_wtime()-start_time<<endl;
				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				// add other nsky
				// vector<id> nsky_for_all_orders_global;
				// compute_other_nsky(cfg, nsky_for_all_orders_global, this->candidates, orders_stack[0]);
			}
		}
	}
}

void dySky::compute_view_recursively_md(Config* cfg, int niveau, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, uint64_t *storage){
	
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
			/////
			//destroy Point pointers
		for (auto partition : partitions)
		for ( auto p : partition.second)
		delete p;
			//	
  	 	
  	 	// cout << "print ids of notSky: " <<endl; 
		// for (auto tuple : notSky){
		// 	cout << tuple << endl;
		// }
		
		*storage=*storage+notSky.size();
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


	bool* cSkyline=new bool[cfg->dataset_size];
	for (int i=0;i<cfg->dataset_size;i++) cSkyline[i]=true;

	for (int i=0;i<preference_cross.size();i++){ // loop on orders
		order_tree *ot=sky_view[preference_cross[i][0]];
		//cout << preference_cross[i][0].first<< ":"<< preference_cross[i][0].second <<endl;
		int j=1;
		while (j<cfg->dyDim_size){ // enter here when multi po dimension
			ot=ot->order_child[preference_cross[i][j]];
			//cout << preference_cross[i][j].first<< ":"<< preference_cross[i][j].second <<endl;
			j++;
		}
		
		for ( int id_tuple : ot->ids) cSkyline[id_tuple]=false;		
	}
	// remove never_sky 
	
	for (int id_tuple : this->never_sky) cSkyline[id_tuple]=false;			

	// compute positive skyline

	vector<id> skyline_result=vector<id>(cfg->dataset_size);	
	int index=0;

	for (int i=0; i<cfg->dataset_size; i++){
		if (cSkyline[i]==true){
			skyline_result[index]=i;
			index++;
		}
		
	}
	delete[] cSkyline;
	skyline_result.resize(index);

	return skyline_result;
}



