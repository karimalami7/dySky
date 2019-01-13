/*
 * graph.h
 *
 *  Created on: December 28, 2018
 *      Author: karim
 */

template <typename T>
class Graph
{
public:
	vector<T> vertices;
	unordered_map<T,unordered_set<T> > out_edges;

	void add_vertices(vector<T> vertices);
	void add_outedges(T vertex_source, unordered_set<T> vertices_destination);
	unordered_map<T,unordered_set<T> > get_edges();
	void print_vertices();
	void print_edges();
	void compute_transitive_closure(Graph<T> p);
};

template <typename T>
void Graph<T>::add_vertices(vector<T> vertices){
	this->vertices=vertices;
}

template <typename T>
void Graph<T>::add_outedges(T vertex_source, unordered_set<T> vertices_destination){
	auto it = this->out_edges.find(vertex_source);
	if(it!=this->out_edges.end()){
		for (auto it2=vertices_destination.begin();it2!=vertices_destination.end();it2++){
			it->second.insert(*it2);
		}
	}
	else{
		this->out_edges.insert(pair<T,unordered_set<T> >(vertex_source,vertices_destination));
	}
}
template <typename T>
unordered_map<T,unordered_set<T> > Graph<T>::get_edges(){
	return this->out_edges;
}
template <typename T>
void Graph<T>::print_vertices(){
	cout <<"vertices: ";
	for (int i=0; i<this->vertices.size(); i++){
		cout<<this->vertices[i]<<" ";
	}
	cout << endl;
}
template <typename T>
void Graph<T>::print_edges(){
	cout <<"edges: "<<endl;
	for (auto it=this->out_edges.begin(); it!=this->out_edges.end(); it++){
		
		cout<<"source: "<<it->first<< endl;

		for (auto it2=it->second.begin(); it2!=it->second.end(); it2++){
			cout<<"-->: "<<(*it2)<< endl;
		}
	}
	cout << endl;
}

void recursive_add(unordered_map<int,unordered_set<int> > &preference_in, id v_src, unordered_set<int> &vertices_dest){
	
	for (auto it=preference_in[v_src].begin();it!=preference_in[v_src].end();it++){
		vertices_dest.insert((*it));
		if (preference_in.find((*it))!=preference_in.end()){
			recursive_add(preference_in,(*it),vertices_dest);
		}
	}

}


template <typename T>
void Graph<T>::compute_transitive_closure(Graph<T> p){

	//unordered_map<int,unordered_set<int> > preference_in=p.get_edges();
	this->vertices=p.vertices;
	for (auto it=p.out_edges.begin();it!=p.out_edges.end();it++){
		unordered_set<int> vertices_dest;
		recursive_add(p.out_edges,it->first,vertices_dest);
		this->add_outedges(it->first,vertices_dest);
	}
	cout << "version transtive du graphe"<<endl;
	this->print_edges();
}