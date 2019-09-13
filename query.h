/*
 * query.h
 *
 *  Created on: Avril 24, 2019
 *      Author: karim
 */

class Query {

public:

	vector<Preference> preference;
	vector<vector<Order>> preference_orders;
	vector<vector<Order>> preference_orders_cross;
	
	void generate_preference(Config* cfg);
	void graph_to_orderPairs(Config* cfg);
	void cross_orders_over_dimensions(Config* cfg);
	void recur_cross(Config* cfg, vector<Order> v, int niv);
};

void Query::generate_preference(Config* cfg){

	// generate a dag using the values of the po domain

	// 1. choose randomely two values vi vj

	// 2. check if the preference + (vi,vj) is still a DAG

		// 2.1 if yes, store it

		// 2.2 if no, throw it

	// repeat 1 and 2 until the preference reaches the wanted properties

	preference=vector<Preference>(cfg->dyDim_size);

	for (int i=0; i<cfg->dyDim_size; i++ ){
	
		int iteration=0;
		int max_iteration=rand()%(2*cfg->dyDim_val);
		//int max_iteration=1;

		this->preference[i].add_vertices(cfg->dyDim_val);

		while (iteration < max_iteration || this->preference[i].out_edges.size()==0) {
			int v1=rand()%cfg->dyDim_val;
			int v2=rand()%cfg->dyDim_val;
			while (v1==v2){
				v2=rand()%cfg->dyDim_val;
			}  
			if(this->preference[i].is_DAG(Order(v1,v2))){
				this->preference[i].add_outedge(v1,v2);
			}
			iteration++;
		}
		this->preference[i].compute_transitive_closure(this->preference[i]);
	}
}

void Query::graph_to_orderPairs(Config* cfg){

	this->preference_orders=vector<vector<Order>>(cfg->dyDim_size);
	for (int i=0; i<cfg->dyDim_size; i++){
		for (auto it=this->preference[i].out_edges.begin(); it!=this->preference[i].out_edges.end();it++){
			for (auto it2=it->second.begin(); it2!=it->second.end(); it2++){
				this->preference_orders[i].push_back(Order((it->first),(*it2)));
			}
		}		
	}
}

void Query::recur_cross(Config* cfg, vector<Order> v, int niv){

	for (int i=0; i<this->preference_orders[niv].size()+1; i++){
		if (i<this->preference_orders[niv].size()) {
			v.push_back(this->preference_orders[niv][i]);
		}else{
			v.push_back(Order(-1,-1));
		}
		
		if (niv<cfg->dyDim_size-1){
			recur_cross(cfg, v, niv+1);
			}
		if (niv==cfg->dyDim_size-1){
			this->preference_orders_cross.push_back(v);
		}
		v.pop_back();
	}
}

void Query::cross_orders_over_dimensions(Config* cfg){
	
	vector<Order> v;
	recur_cross(cfg, v, 0);

}

