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
		int max_iteration=rand()%(5*cfg->dyDim_val);
		
		this->preference[i].add_vertices(cfg->dyDim_val);

		//while (iteration < max_iteration || this->preference[i].out_edges.size()==0) {
		while (iteration < max_iteration) {
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

