/*
 * graph.h
 *
 *  Created on: December 28, 2018
 *      Author: karim
 */

template <typename T>
class Graph
{
protected:
	vector<T> vertices;
	unordered_map<T,vector<T> > out_edges;

public:
	void add_vertices(vector<T> vertices);
	void add_outedges(T vertex_source, vector<T> vertices_destination);
	unordered_map<T,vector<T> > get_edges();
	void print_vertices();
	void print_edges();
};

template <typename T>
void Graph<T>::add_vertices(vector<T> vertices){
	this->vertices=vertices;
}

template <typename T>
void Graph<T>::add_outedges(T vertex_source, vector<T> vertices_destination){
	auto it = this->out_edges.find(vertex_source);
	if(it!=this->out_edges.end()){
		for (int i=0;i<vertices_destination.size();i++){
			it->second.push_back(vertices_destination[i]);
		}
	}
	else{
		this->out_edges=pair<T,vector<T> >(vertex_source,vertices_destination);
	}
}
template <typename T>
unordered_map<T,vector<T> > Graph<T>::get_edges(){
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