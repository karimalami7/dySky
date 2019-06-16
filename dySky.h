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
	//void compute_always_skyline(Config* cfg);
	int compute_candidates(Config* cfg);
	
	void print_dataset (Config* cfg);
	void print_dataset (vector<Point> data, string name, int d);
	dySky(Config *cfg);
	
	void compute_views(Config* cfg, vector<vector<Order>> preference_orders);
	void compute_view_1d(Config* cfg, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders);
	void compute_view_recursively_md(Config* cfg, int niveau, vector<Point> &dataset, vector<id> &sky_for_all_orders_global,vector<Order> orders_stack, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, vector<id> remaining_candidates);
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

// void dySky::compute_always_skyline(Config* cfg){
// 	//cout << "dySky::compute_always_skyline"<<endl;

// 	int All = (1<<cfg->statDim_size)-1;
// 	vector<Space> full_Space;
// 	listeAttributsPresents(All, cfg->statDim_size, full_Space);
//     this->always_sky=subspaceSkylineSize_TREE(full_Space, to_dataset);
//     //cout << "Always Skyline size: "<<this->always_sky.size()<<endl;
//     // cerr << "::" <<endl;
//     // for (int k=0; k< this->always_sky.size(); k++) cerr << this->always_sky[k]<<endl; 
// 	// int All = (1<<cfg->statDim_size+1)-1;
// 	// vector<Space> full_Space;
// 	// listeAttributsPresents(All, cfg->statDim_size+1, full_Space);
// 	// vector<Point> temp_dataset;
// 	// for (int i=0; i< cfg->dataset_size; i++){
// 	// 	Point p=(int*)malloc((cfg->statDim_size+2)*sizeof(int));
// 	// 		for (int j=0;j<=cfg->statDim_size;j++){
// 	// 			p[j]=to_dataset[i][j];
// 	// 		}
// 	// 		p[cfg->statDim_size+1]=po_dataset[i];
// 	// 		temp_dataset.push_back(p);
// 	// }
// 	// // gooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooood      results
//  //    this->always_sky=subspaceSkylineSize_TREE(full_Space, temp_dataset);
//  //    //cout << "Always Skyline size: "<<this->always_sky.size()<<endl;
//  //    cerr << "::" <<endl;
//  //    for (int k=0; k< this->always_sky.size(); k++) cerr << this->always_sky[k]<<endl; 

// }

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
	//cout << "dySky::compute_candidates"<<endl;

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

    for (auto cand: this->candidates){
    	cout <<cand << " : ";
    } 
    cout <<endl;

    // compute never_sky
    // vector<id> all_ids;
    // for (int i=0; i<cfg->dataset_size; i++){
    // 	all_ids.push_back(i);
    // }

    // never_sky = all_ids minus sky_union_candidates
   //  std::vector<int> never_sky(cfg->dataset_size);  //filled with zeros
  	// std::vector<int>::iterator it_never;
  	// it_never=std::set_difference (all_ids.begin(),all_ids.end(), 
  	// 	this->candidates.begin(),this->candidates.end(), never_sky.begin()); // diff between all ids and dom objects
  	// never_sky.resize(it_never-never_sky.begin());  
   //  this->never_sky.swap(never_sky);
    //cout <<"Never Skyline size: "<< this->never_sky.size()<<endl;

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


	


	vector<Point> candidates_tuples;
	for (int i=0; i<this->candidates.size(); i++){
		Point p=(int*)malloc((cfg->statDim_size+cfg->dyDim_size+1)*sizeof(int));
		for (int j=0;j<=cfg->statDim_size;j++){
			p[j]=to_dataset[this->candidates[i]][j];
		}
		for (int j=0;j<cfg->dyDim_size;j++){
			p[j+(cfg->statDim_size+1)]=po_dataset[this->candidates[i]][j];
		}
		//p[cfg->statDim_size+1]=po_dataset[aSky_union_candidates[i]][0];
		candidates_tuples.push_back(p);
	}

	//print_dataset(aSky_union_candidates_tuples, " compute_view2 data", cfg->statDim_size+cfg->dyDim_size);

	if (cfg->dyDim_size==1){
		compute_view_1d(cfg, candidates_tuples, this->sky_view, preference_orders);
	}
	else{
		vector<Point> to_candidates_tuples(this->candidates.size());
		vector<vector<int>> po_candidates_tuples(this->candidates.size());
		for (int i=0; i<this->candidates.size(); i++){
			to_candidates_tuples[i]=to_dataset[this->candidates[i]];
			po_candidates_tuples[i]=po_dataset[this->candidates[i]];
		}
		vector<Order> orders_stack;
		//cout << "-->dySky::compute_view_recursively"<<endl;


		compute_view_recursively_md(cfg, 0, candidates_tuples, this->candidates,orders_stack, this->sky_view, preference_orders, this->candidates);
	}
	// int views_total_storage=0;
	// for (auto it=this->sky_view.begin();it!=this->sky_view.end();it++){
	// 	views_total_storage+=(it->second)->ids.size();
	// }
	// //cout << "Views total storage: " << views_total_storage <<endl;

}

void dySky::compute_view_1d(Config* cfg, vector<Point> &dataset, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders){
	
	for(int i=0;i<preference_orders[0].size();i++){
		int best_value=preference_orders[0][i].first;
		int worst_value=preference_orders[0][i].second;
		cout << best_value << ":" << worst_value <<endl;

		vector<Point> branch_dataset;
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// partitionner les donnees par rapport à cette dimension
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

		//print_dataset(branch_dataset, " compute_view2 data", cfg->statDim_size+cfg->dyDim_size);

		//*************************************
		// calculer le skyline
		int All = (1<<(cfg->statDim_size+cfg->dyDim_size))-1;
		vector<Space> full_Space;
		listeAttributsPresents(All, cfg->statDim_size+cfg->dyDim_size, full_Space);
		vector<id> sky=subspaceSkylineSize_TREE(full_Space, branch_dataset);

		//**************************************
		//remove always skyline ids 
		unordered_set<id> view;
		view.insert(sky.begin(),sky.end());  
	  	// std::vector<id>::iterator it_view;
	  	// std::sort (view_sky.begin(),view_sky.end());
	  	// it_view=std::set_difference (view_sky.begin(),view_sky.end(), 
	  	// 	this->always_sky.begin(),this->always_sky.end(), view_candidates.begin()); 
	  	// view_candidates.resize(it_view-view_candidates.begin());  
	  	cout << "view before add candidates"<<endl;
	  	for (int m : view){
			cout <<m<<":";
		}
		cout <<endl;
	  	//***************************************
	  	//add other candidates
	  	for (int i=0; i< this->candidates.size();i++){
	  		if( po_dataset[this->candidates[i]][0]!=best_value && po_dataset[this->candidates[i]][0]!=worst_value){
	  			view.insert(this->candidates[i]);
	  		}
	  	}

	  	cout << "view after add candidates"<<endl;
	  	for (int m : view){
			cout <<m<<":";
		}
		cout <<endl;
	  	//******************************************
	  	// sort ids
	  	vector<id> view_vector;
	  	view_vector.insert(view_vector.begin(),view.begin(),view.end());
	  	std::sort (view_vector.begin(),view_vector.end());
	  	sky_view[Order(best_value,worst_value)]=new order_tree;
	  	sky_view[Order(best_value,worst_value)]->ids=view_vector;
	}
}

void dySky::compute_view_recursively_md(Config* cfg, int niveau, vector<Point> &dataset, vector<id> &sky_for_all_orders_global, vector<Order> orders_stack, unordered_map<Order, order_tree*, pairhash> &sky_view, vector<vector<Order>> preference_orders, vector<id> remaining_candidates){
	
	
	for(int s=0;s<preference_orders[niveau].size();s++){
		int best_value=preference_orders[niveau][s].first;
		int worst_value=preference_orders[niveau][s].second;
		cout <<"niveau: "<<niveau<< " --> "<< best_value<<" : "<<worst_value<<endl;		


		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// FOR ALL DIMENSIONS
		unordered_map<int,Point> to_map;
		unordered_map<int,vector<int>> po_map;
		vector<id> directly_sky;
		// for (int i=0; i<this->po_dataset.size();i++){
		// 	cout << i <<" : " <<this->po_dataset[i][0]<<" : "<< this->po_dataset[i][1]<<endl;
		// }
		for (int i=0; i<this->candidates.size(); i++){
			//cout <<this->candidates[i]<<" : "<<po_dataset[this->candidates[i]][niveau]<< " : "<<best_value <<" : "<<worst_value<<endl;
			if (this->po_dataset[this->candidates[i]][niveau]==best_value){
				//cout << this->candidates[i]<<endl;
				Point p=(int*)malloc((cfg->statDim_size+2)*sizeof(int));
				vector<int> t;
				for (int j=0;j<=cfg->statDim_size;j++) p[j]=this->to_dataset[this->candidates[i]][j];
				p[cfg->statDim_size+1]=0;
				for (int j=0;j<cfg->dyDim_size;j++) if (j!=niveau) t.push_back(this->po_dataset[this->candidates[i]][j]);
				to_map.insert(pair<int,Point>(p[0],p) );
				po_map.insert(pair<int,vector<int>>(p[0],t));
				//candidate_partition.push_back(to_partial_dataset[i][0]);
			}else if (this->po_dataset[this->candidates[i]][niveau]==worst_value){
				//cout << this->candidates[i]<<endl;
				Point p=(int*)malloc((cfg->statDim_size+2)*sizeof(int));
				vector<int> t;
				for (int j=0;j<=cfg->statDim_size;j++) p[j]=this->to_dataset[this->candidates[i]][j];
				p[cfg->statDim_size+1]=1;
				for (int j=0;j<cfg->dyDim_size;j++) if (j!=niveau) t.push_back(this->po_dataset[this->candidates[i]][j]);
				to_map.insert(pair<int,Point>(p[0],p) );
				po_map.insert(pair<int,vector<int>>(p[0],t));
				//candidate_partition.push_back(to_partial_dataset[i][0]);
			}
			else{
				directly_sky.push_back(this->candidates[i]);
			}
		}
		map<vector<int>, vector<Point>> partitions;
		for (auto it=to_map.begin(); it!=to_map.end(); it++){
			partitions[po_map[it->first]].push_back(it->second);
		}
		cout << "nombre partition"<<partitions.size()<<endl;
		// compute skyline for every partition -> always skyline + candidates
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

	    sort(directly_sky.begin(),directly_sky.end());
	  	cout <<"directly_sky"<<endl;
	  	for (auto e : directly_sky){
	  		cout <<e<<" : ";
	  	}
	  	cout <<endl;
	    sort(sky_for_all_orders.begin(),sky_for_all_orders.end());
		//UNION
		std::vector<id> v(this->to_dataset.size());   
	  	std::vector<id>::iterator it;
	  	it=std::set_union(sky_for_all_orders.begin(), sky_for_all_orders.end(), directly_sky.begin(), directly_sky.end(), v.begin());                        
	  	v.resize(it-v.begin());   
	  	sky_for_all_orders=v; 

	  	cout <<"sky_for_all_orders"<<endl;
	  	for (auto e : sky_for_all_orders){
	  		cout <<e<<" : ";
	  	}
	  	cout <<endl;

	  	// //INTERSECTION
	  	std::vector<id> l(this->to_dataset.size());   
	  	std::vector<id>::iterator it2;
	  	it2=std::set_intersection(sky_for_all_orders.begin(), sky_for_all_orders.end(), sky_for_all_orders_global.begin(), sky_for_all_orders_global.end(), l.begin());                        
	  	l.resize(it2-l.begin());   
	  	sky_for_all_orders=l; 

	  	cout <<"sky_for_all_orders after global"<<endl;
	  	for (auto e : sky_for_all_orders){
	  		cout <<e<<" : ";
	  	}
	  	cout <<endl;

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
			orders_stack.push_back(Order(best_value,worst_value));
			sky_view[Order(best_value,worst_value)]=new order_tree;;
			this->compute_view_recursively_md(cfg,  niveau+1, branch_dataset, sky_for_all_orders, orders_stack, sky_view[Order(best_value,worst_value)]->order_child, preference_orders, remaining_candidates);
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
			vector<id> branch_dataset_ids(branch_dataset.size());
			for (int b=0;b<branch_dataset.size();b++){
				branch_dataset_ids[b]=branch_dataset[b][0];
			}
			sort(branch_dataset_ids.begin(),branch_dataset_ids.end());
		  	cout <<"branch_dataset ids : "<<endl;
		  	for (auto e : branch_dataset_ids){
		  		cout <<e<<" : ";
		  	}
		  	cout <<endl;

			orders_stack.push_back(Order(best_value,worst_value));
			//*************************************
			// calculer le skyline
			int All = (1<<(cfg->statDim_size+cfg->dyDim_size))-1;
			vector<Space> full_Space;
			listeAttributsPresents(All, cfg->statDim_size+cfg->dyDim_size, full_Space);
			vector<id> view=subspaceSkylineSize_TREE(full_Space, branch_dataset);
		  	cout <<"view limited"<<endl;
		  	for (auto e : view){
		  		cout <<e<<" : ";
		  	}
		  	cout <<endl;
			//


			sort(view.begin(),view.end());
			//MINUS
			std::vector<id> notSky(branch_dataset_ids.size());   
		  	std::vector<id>::iterator it4;
		  	it4=std::set_difference(branch_dataset_ids.begin(), branch_dataset_ids.end(), view.begin(), view.end(), notSky.begin());                        
		  	notSky.resize(it4-notSky.begin()); 
		  	cout <<"notsky"<<endl;
		  	for (auto e : notSky){
		  		cout <<e<<" : ";
		  	}
		  	cout <<endl;  
		  	//

		  	//***************************************
		  	//add other candidates
		  	// for (int i=0; i< sky_for_all_orders_global.size();i++){
		  	// 	bool to_add=false;
		  	// 	for(int d=0;d<cfg->dyDim_size;d++){
		  	// 		if( po_dataset[sky_for_all_orders_global[i]][d]!=orders_stack[d].first && po_dataset[sky_for_all_orders_global[i]][d]!=orders_stack[d].second){
		  	// 		to_add=true;
		  	// 		}
		  	// 	}
		  	// 	if(to_add){
		  	// 		view.push_back(sky_for_all_orders_global[i]);
		  	// 	}
		  	// }

		  	//******************************************
		  	// sort ids
		  	std::sort (view.begin(),view.end());
			//UNION
			std::vector<id> h(this->to_dataset.size());   
		  	std::vector<id>::iterator it3;
		  	it3=std::set_union(view.begin(), view.end(), sky_for_all_orders.begin(), sky_for_all_orders.end(), h.begin());                        
		  	h.resize(it3-h.begin());   
		  	view=h; 
		  	cout <<"view union sky for all orders"<<endl;
		  	for (auto e : view){
		  		cout <<e<<" : ";
		  	}
		  	cout <<endl;

		  	//MINUS
			std::vector<id> g(view.size());   
		  	std::vector<id>::iterator it5;
		  	it5=std::set_difference(view.begin(), view.end(), notSky.begin(), notSky.end(), g.begin());                        
		  	g.resize(it5-g.begin());   
		  	view=g; 
		  	cout <<"view minus not sky for this order"<<endl;
		  	for (auto e : view){
		  		cout <<e<<" : ";
		  	}
		  	cout <<endl;

		  	//
		  	//
		  	sky_view[Order(best_value,worst_value)]=new order_tree;
		  	sky_view[Order(best_value,worst_value)]->ids=view;

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

	for (int i=0;i<preference.size();i++){ // loop on orders
		order_tree *ot=sky_view[preference[i][0]];
		int j=1;
		while (j<cfg->dyDim_size){ // enter here when multi po dimension
			ot=ot->order_child[preference[i][j]];
			j++;
		}
		cout << "print ids for one preference" <<endl; 
		for (auto tuple : ot->ids){
			cout << tuple << endl;
		}
		cout <<endl;
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



	//result.insert(result.end(), always_sky.begin(), always_sky.end());

	cout << "print ids at the end" <<endl; 
	for (auto tuple : result){
		cout << tuple << endl;
	}

	return result;
}
