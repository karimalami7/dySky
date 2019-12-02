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
	vector<vector<int>> paths(Config* cfg);
	void aux_path(int value, vector<int> path, vector<vector<int>> &paths);
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
	
	this->add_vertices(cfg->dyDim_val);

	//random preference setting

	double density=0;
	int max_number_edges=(cfg->dyDim_val*(cfg->dyDim_val-1))/2;

	int iteration=0;
	//int max_iteration=rand()%(10*cfg->dyDim_val);
	int number_edges;
	while (density < 0.8) {
		int v1=rand()%cfg->dyDim_val;
		int v2=rand()%cfg->dyDim_val;
		while (v1==v2){
			v2=rand()%cfg->dyDim_val;
		}  
		if(this->is_DAG(Order(v1,v2), &number_edges)){
			this->add_outedge(v1,v2);
			density=(float)number_edges/(float)max_number_edges;
			cout << "density" << density <<endl;
		}
		
	}

	//manual preference setting

	// this->add_outedge(0,1);
	// this->add_outedge(1,2);
	
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
	
	std::vector<int> heads(this->vertices.size());                      
	std::vector<int>::iterator it;
	it=std::set_difference (this->vertices.begin(), this->vertices.end(), vertices_with_incoming_edge.begin(), vertices_with_incoming_edge.end(), heads.begin());
	heads.resize(it-heads.begin());

	cout <<"heads of the preference: "<<endl;
	for (auto vertex: heads){
		cout << "vertex: "<< vertex<< endl;
	}
	cout <<endl;

	return heads;
}

void Preference::aux_path(int value, vector<int> path, vector<vector<int>> &paths){
	if (this->out_edges.find(value)==this->out_edges.end()){
		paths.push_back(path);
	}
	else {
		for (auto val: out_edges[value]){
			path.push_back(val);
			aux_path(val, path, paths);	
			path.pop_back();
		}
	}
}

vector<vector<int>> Preference::paths(Config* cfg){

	vector<int> heads=this->find_heads(cfg);
	vector<vector<int>> paths;
	int niveau=0;
	for (auto head: heads){
		vector<int> path;
		path.push_back(head);
		aux_path(head, path, paths);
	}
	cout <<"paths of the preference: "<<endl;
	for (auto path: paths){
		cout << "path: ";
		for (auto val: path) {
			cout << val << " ";
		}
		cout <<endl;
	}
	cout <<endl;
	return paths;
}

