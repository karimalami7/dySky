class Preference : public Graph<int>{
public:
void add_vertices(int max_id);
void add_edges(id vertex_src, unordered_set<id> vertices_dest);
bool is_DAG(Order new_edge);

};

void Preference::add_vertices(int max_id){
	for (int i=0; i<max_id; i++){
		this->vertices.push_back(i);
	}
}

void Preference::add_edges(id vertex_src, unordered_set<id> vertices_dest){
		this->out_edges.insert(pair<id,unordered_set<id>>(vertex_src,vertices_dest));
}

bool Preference::is_DAG(Order new_edge){

	Preference p_trans;
	p_trans.compute_transitive_closure(*this);

	// detect cycle
	bool cycle_exists=false;

	auto it_src=p_trans.out_edges.find(new_edge.second);
	if (it_src!=p_trans.out_edges.end()){
		auto it_dest=it_src->second.find(new_edge.first);
		if(it_dest!=it_src->second.end()){
			cycle_exists=true;
		}
	}

	return !cycle_exists;

}
