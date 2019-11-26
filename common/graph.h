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
	vector<int> vertex_color;
	
	// Setters 
	void set_vertices(vector<T> vertices);
	void set_outedges(T vertex_source, unordered_set<T> vertices_destination);
	// Getters
	vector<T> get_vertices();
	unordered_map<T,unordered_set<T> > get_edges();
	
	// Methods	
	void add_outedge(T vertex_source, T vertex_destination);
	void print_vertices();
	void print_edges();
	void compute_transitive_closure(Graph<T> p);
	void greedyColoring(); 
	bool is_subgraph(Graph<T> p); 
	bool is_DAG(Order new_edge);

};


// Setters 
template <typename T>
void Graph<T>::set_vertices(vector<T> vertices){
	this->vertices=vertices;
}

template <typename T>
void Graph<T>::set_outedges(T vertex_source, unordered_set<T> vertices_destination){
	auto it = this->out_edges.find(vertex_source);
	// if vertex_source has already an out_edge
	if(it!=this->out_edges.end()){
		for (auto it2=vertices_destination.begin();it2!=vertices_destination.end();it2++){
			it->second.insert(*it2);
		}
	}
	else{
		this->out_edges.insert(pair<T,unordered_set<T> >(vertex_source,vertices_destination));
	}
}
// Getters
template <typename T>
vector<T> Graph<T>::get_vertices(){
	return this->vertices;
}

template <typename T>
unordered_map<T,unordered_set<T> > Graph<T>::get_edges(){
	return this->out_edges;
}

// Methods
template <typename T>
void Graph<T>::add_outedge(T vertex_source, T vertex_destination){
	auto it = this->out_edges.find(vertex_source);
	if(it!=this->out_edges.end()){
		it->second.insert(vertex_destination);
	}
	else{
		unordered_set<T> vertices_destination;
		vertices_destination.insert(vertex_destination);
		this->out_edges.insert(pair<T,unordered_set<T> >(vertex_source,vertices_destination));
	}
}

template <typename T>
void Graph<T>::print_vertices(){
	//cout <<"vertices: ";
	for (int i=0; i<this->vertices.size(); i++){
		//cout<<this->vertices[i]<<" ";
	}
	//cout << endl;
}
template <typename T>
void Graph<T>::print_edges(){
	//cout <<"edges: "<<endl;
	for (auto it=this->out_edges.begin(); it!=this->out_edges.end(); it++){
		
		cout<<"source: "<<it->first<< endl;

		for (auto it2=it->second.begin(); it2!=it->second.end(); it2++){
			cout<<"-->: "<<(*it2)<< endl;
		}
	}
	//cout << endl;
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
		this->set_outedges(it->first,vertices_dest);
	}
}

// Assigns colors (starting from 0) to all vertices and prints 
// the assignment of colors 
template <typename T>
void Graph<T>::greedyColoring() 
{ 	
	int V;    // No. of vertices 
    list<int> *adj;    // A dynamic array of adjacency lists 

    V=this->vertices.size();
    adj = new list<int>[V];

    for (auto it_umap=this->out_edges.begin();it_umap!=this->out_edges.end();it_umap++){
    	for (auto it_uset=it_umap->second.begin();it_uset!=it_umap->second.end();it_uset++){
    		adj[it_umap->first].push_back(*it_uset);
    	}	
    }

    int result[V]; 
  
    // Assign the first color to first vertex 
    result[0]  = 0; 
  
    // Initialize remaining V-1 vertices as unassigned 
    for (int u = 1; u < V; u++) 
        result[u] = -1;  // no color is assigned to u 
  
    // A temporary array to store the available colors. True 
    // value of available[cr] would mean that the color cr is 
    // assigned to one of its adjacent vertices 
    bool available[V]; 
    for (int cr = 0; cr < V; cr++) 
        available[cr] = false; 
  
    // Assign colors to remaining V-1 vertices 
    for (int u = 1; u < V; u++) 
    { 
        // Process all adjacent vertices and flag their colors 
        // as unavailable 
        list<int>::iterator i; 
        for (i = adj[u].begin(); i != adj[u].end(); ++i) 
            if (result[*i] != -1) 
                available[result[*i]] = true; 
  
        // Find the first available color 
        int cr; 
        for (cr = 0; cr < V; cr++) 
            if (available[cr] == false) 
                break; 
  
        result[u] = cr; // Assign the found color 
  
        // Reset the values back to false for the next iteration 
        for (i = adj[u].begin(); i != adj[u].end(); ++i) 
            if (result[*i] != -1) 
                available[result[*i]] = false; 
    } 
  
    // print the result 
    for (int u = 0; u < V; u++){ 
        ////cout << "Vertex " << u << " --->  Color "<< result[u] << endl; 
        this->vertex_color.push_back(result[u]);
    }
} 


template <typename T>
bool Graph<T>::is_subgraph(Graph<T> p) 
{	
	for(auto it_outedges=this->out_edges.begin();it_outedges!=this->out_edges.end();it_outedges++){
		auto it_source_edge=p.out_edges.find(it_outedges->first);
		if (it_source_edge==this->out_edges.end()){
			return false;
		}
		else {
			for (auto it_dests = it_outedges->second.begin(); it_dests != it_outedges->second.end(); it_dests++){
				auto dest=it_source_edge->second.find(*it_dests);
				if (dest==it_source_edge->second.end()){
					return false;
				}
			}
		}
	}
	return true;
}

template <typename T>
bool Graph<T>::is_DAG(Order new_edge){

	Graph<T> p_trans;
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