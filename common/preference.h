/*
 * preference.h
 *
 *  Created on: January 5, 2019
 *      Author: karim
 */

class Preference : public Graph<int>{

	public:

	void add_vertices(int max_id);
	void add_edges(id vertex_src, unordered_set<id> vertices_dest);
	void generate_preference(Config* cfg);
	vector<int> find_heads(Config* cfg);
	void find_paths();
};

void Preference::add_vertices(int max_id){
	for (int i=0; i<max_id; i++){
		this->vertices.push_back(i);
	}
}

void Preference::add_edges(id vertex_src, unordered_set<id> vertices_dest){
		this->out_edges.insert(pair<id,unordered_set<id>>(vertex_src,vertices_dest));
}

void Preference::generate_preference(Config* cfg){

	// generate a dag using the values of the po domain

	// 1. choose randomely two values vi vj

	// 2. check if the preference + (vi,vj) is still a DAG

		// 2.1 if yes, store it

		// 2.2 if no, throw it

	// repeat 1 and 2 until the preference reaches the wanted properties
	
	int iteration=0;
	int max_iteration=rand()%(5*cfg->dyDim_val);
	
	this->add_vertices(cfg->dyDim_val);

	//while (iteration < max_iteration || this->preference[i].out_edges.size()==0) {
	while (iteration < max_iteration) {
		int v1=rand()%cfg->dyDim_val;
		int v2=rand()%cfg->dyDim_val;
		while (v1==v2){
			v2=rand()%cfg->dyDim_val;
		}  
		if(this->is_DAG(Order(v1,v2))){
			this->add_outedge(v1,v2);
		}
		iteration++;
	}
	this->compute_transitive_closure(*this);
}

vector<int> Preference::find_heads(Config* cfg){
	// find the vertices that have no incoming edge

	unordered_set<int> set_vertices_with_incoming_edge;

	for (auto vertex_edges: this->out_edges){
		for (auto vertex: vertex_edges.second){
			set_vertices_with_incoming_edge.insert(vertex);
		}
	}

	vector<int> vertices_with_incoming_edge;
	vertices_with_incoming_edge.insert(vertices_with_incoming_edge.begin(), set_vertices_with_incoming_edge.begin(), set_vertices_with_incoming_edge.end());
	sort(vertices_with_incoming_edge.begin(), vertices_with_incoming_edge.end());
	
	std::vector<int> v(this->vertices.size());                      
	std::vector<int>::iterator it;
	it=std::set_difference (this->vertices.begin(), this->vertices.end(), vertices_with_incoming_edge.begin(), vertices_with_incoming_edge.end(), v.begin());
	v.resize(it-v.begin());

	// for (auto vertex: vertices_with_incoming_edge){
	// 	cout << "vertex: "<< vertex<< endl;
	// }

	return vertices_with_incoming_edge;
}