class Preference : public Graph<int>{
public:
void add_vertices(int max_id);
void print_vertices();
};

void Preference::add_vertices(int max_id){
	for (int i=0; i<max_id; i++){
		this->vertices.push_back(i);
	}
}

void Preference::print_vertices(){
	for (int i=0; i<this->vertices.size(); i++){
		cout <<this->vertices[i]<<endl;
	}
	cout << endl;
}
