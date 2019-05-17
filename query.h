/*
 * query.h
 *
 *  Created on: Avril 24, 2019
 *      Author: karim
 */

class Query {

public:
	vector<Preference> preference;

	void generate_preference(Config* cfg);
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
	
		int size_preference=0;

		this->preference[i].add_vertices(cfg->dyDim_val);

		while (size_preference < cfg->dyDim_val*2) {
			int v1=rand()%cfg->dyDim_val;
			int v2=rand()%cfg->dyDim_val;
			while (v1==v2){
				v2=rand()%cfg->dyDim_val;
			}  
			if(this->preference[i].is_DAG(Order(v1,v2))){
				this->preference[i].add_outedge(v1,v2);
			}
			size_preference++;
		}		
	}


}