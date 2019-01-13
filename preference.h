class Preference : public Graph<int>{
public:
void add_vertices(int max_id);
void add_edges(id vertex_src, vector<id> vertices_dest);

};

void Preference::add_vertices(int max_id){
	for (int i=0; i<max_id; i++){
		this->vertices.push_back(i);
	}
}

void Preference::add_edges(id vertex_src, vector<id> vertices_dest){
		this->out_edges.insert(pair<id,vector<id>>(vertex_src,vertices_dest));
}

