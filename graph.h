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