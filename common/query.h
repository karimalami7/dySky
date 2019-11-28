/*
 * query.h
 *
 *  Created on: Avril 24, 2019
 *      Author: karim
 */

class Query {

public:

	vector<Preference> preference; // for every dimension, we have a DAG
	vector<vector<Order>> preference_orders; // for every dimension, we have a set of "Order
	vector<vector<chain>> preference_chains;
	vector<vector<Order>> preference_orders_cross; // we cross orders between dimensions
	
	Query();
	Query(Config* cfg);
	void graph_to_orderPairs(Config* cfg);
	void cross_orders_over_dimensions(Config* cfg);
	void recur_cross(Config* cfg, vector<Order> v, int niv);
};
Query::Query(){
	
}
Query::Query(Config* cfg){

	preference= vector<Preference>(cfg->dyDim_size);
	preference_chains= vector<vector<chain>>(cfg->dyDim_size);
	
	for (int i=0; i<cfg->dyDim_size; i++ ){
		this->preference[i].generate_preference(cfg);
		// cout <<"preference: "<<endl;
		// this->preference[i].print_edges();
		// cout << endl;
		this->preference[i].compute_transitive_reduction(this->preference[i]);
		// cout <<"preference with transitive reduction: "<<endl;
		// this->preference[i].print_edges();
		// cout << endl;
		preference_chains[i]=this->preference[i].paths(cfg);
		this->preference[i].compute_transitive_closure(this->preference[i]);
		// cout <<"preference with transitive closure: "<<endl;
		// this->preference[i].print_edges();
		// cout << endl;
	}
	this->graph_to_orderPairs(cfg);
	this->cross_orders_over_dimensions(cfg);

}

void Query::graph_to_orderPairs(Config* cfg){

	this->preference_orders=vector<vector<Order>>(cfg->dyDim_size);
	for (int i=0; i<cfg->dyDim_size; i++){
		set<int> ordered_values;
		for (auto it=this->preference[i].out_edges.begin(); it!=this->preference[i].out_edges.end();it++){
			ordered_values.insert(it->first);
			for (auto it2=it->second.begin(); it2!=it->second.end(); it2++){
				this->preference_orders[i].push_back(Order((it->first),(*it2)));
				ordered_values.insert(*it2);
			}
		}
		vector<int>	all_values(cfg->dyDim_val);
		for (int i=0;i<cfg->dyDim_val;i++) all_values[i]=i;

		vector<int> non_ordered_values(all_values.size());   
	  	vector<int>::iterator it;
	  	it=std::set_difference(all_values.begin(), all_values.end(), ordered_values.begin(), ordered_values.end(), non_ordered_values.begin());                        
	  	non_ordered_values.resize(it-non_ordered_values.begin());	

		for (int value : non_ordered_values){ 
			this->preference_orders[i].push_back(Order(value,value));
		}	
	}

	// print number of sequences
	int num_seq=1;
	for (int i=0; i<cfg->dyDim_size; i++){
		cout << "number of orders: "<< this->preference_orders[i].size()<<endl;
		num_seq*=this->preference_orders[i].size();
	}
	cerr << "number of sequences of this query: "<< num_seq<<endl;

}

void Query::recur_cross(Config* cfg, vector<Order> v, int niv){

	for (int i=0; i<this->preference_orders[niv].size(); i++){
		if (i<this->preference_orders[niv].size()) {
			v.push_back(this->preference_orders[niv][i]);
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
	
	if (cfg->dyDim_size==1){
		for (auto order : this->preference_orders[0]){
			vector<Order> v;
			v.push_back(order);
			this->preference_orders_cross.push_back(v);
		}
	}
	else{
		vector<Order> v;
		recur_cross(cfg, v, 0);		
	}


}

