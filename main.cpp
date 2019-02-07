/*
 * main.cpp
 *
 *  Created on: December 19, 2018
 *      Author: karim
 */


#include "common.h"
#include "config.h"
#include "dySky.h"
#include "graph.h"
#include "preference.h"
#include "cps.h"
#include "arg.h"
#include "tos.h"
#include <omp.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////

void printUsage() {
	printf("\ndySky - computes Skyline over datasets with dynamic partial order\n\n");
	printf("USAGE: ./dySky [Option Argument]\n");
	printf( " -n: size of the dataset\n" );
	printf( " -s: number of static dimensions\n" );
	printf( " -k: distinct values of static dimensions\n" );
	printf( " -d: number of dynamic dimensions\n" );
	printf( " -m: distinct values of dynamic dimensions\n" );		
	printf( " -t: run with num_threads, e.g., 20 (default \"4\")\n" );
	printf("Example: ");
	printf("./dySky -n 100000 -s 8 -k 100 -d 2 -m 10 -t 8 \n\n");
}


int main(int argc, char** argv) {

	Config *cfg = new Config;
	cfg->dataset_size=100000;
	cfg->statDim_size=8;
	cfg->statDim_val=20;
	cfg->dyDim_size=1;
	cfg->dyDim_val=4;
	cfg->num_threads=8;
	cfg->verbose=false;

	int c = 0;
	opterr = 0;

	while ((c = getopt(argc, argv, "n:s:k:d:m:t")) != -1) {
		switch (c) {
		case 'n':
			cfg->dataset_size = atoi(optarg);
			break;
		case 's':
			cfg->statDim_size = atoi(optarg);
			break;
		case 'k':
			cfg->statDim_val = atoi(optarg);
			break;
		case 'd':
			cfg->dyDim_size = atoi(optarg);
			break;
		case 'm':
			cfg->dyDim_val = atoi(optarg);
			break;
		case 't':
			cfg->num_threads = atoi(optarg);
			break;
		default:
			if (isprint(optopt))
				fprintf( stderr, "Unknown option `-%c'.\n", optopt);
			printUsage();
			return 1;
		}
	}




	// dySky
	cerr << "=====dySky=====" <<endl;
	// start for pre processing
	cerr << "---preprocessing---"<<endl;
	dySky dysky;
	// generate total order data
	dysky.generate_to_data(cfg);
	// generate partial order data
	dysky.generate_po_data(cfg);
	// compute skyline by BSkyTree on static order dimensions (Always Sky)
	double start_time=omp_get_wtime();
	double start_time2=omp_get_wtime();
	dysky.compute_always_skyline(cfg);
	cerr<<"--> Time for compute_always_skyline: "<< omp_get_wtime()-start_time2 << endl;
	// compute candidates 
	start_time2=omp_get_wtime();
	dysky.compute_candidates(cfg);
	cerr<<"--> Time for compute_candidates: "<< omp_get_wtime()-start_time2 << endl;
	// compute views
	start_time2=omp_get_wtime();
	dysky.compute_views(cfg);
	cerr<<"--> Time for compute_views: "<< omp_get_wtime()-start_time2 << endl;
	dysky.print_dataset(cfg);

	// end of pre processing
	
	cerr<<"--> Time for all dySky: "<< omp_get_wtime()-start_time << endl;


	//ARG  

	// cache queries

	cerr << "=====Arg=====" <<endl;
	start_time=omp_get_wtime();
	// store skylines
	cerr << "---Compute skylines to cache---" <<endl;
	// Arg arg;
	// arg.cached_preferences.push_back(p);
	// Cps *cps_arg = new Cps();
	// cps_arg->decompose_preference(p,cfg);
	// cps_arg->to_dataset=dysky.to_dataset;
	// cps_arg->po_dataset=dysky.po_dataset;
	// cps_arg->compute_skyline(cfg);
	// arg.cached_skylines.push_back(cps_arg->skyline_result);
	cerr<<"--> Time for all ARG: "<< omp_get_wtime()-start_time << endl;


	// TOS

	// compute all total orders

	cerr << "=====TOS=====" <<endl;
	start_time=omp_get_wtime();
	cerr << "---Compute all TO skylines---" <<endl;

  	Tos tos;
  	vector<int> tos_values;
  	for (int i=0; i<cfg->dyDim_val; i++){tos_values.push_back(i);}
  	do {
  	  	Preference p_to;
		p_to.add_vertices(cfg->dyDim_val);
		for (id source=0;source<cfg->dyDim_val-1;source++){
			unordered_set<id> v_dest;
			v_dest.insert(tos_values[source+1]);
			p_to.add_edges(tos_values[source],v_dest);
		}
		Cps *cps_tos = new Cps();
		cps_tos->decompose_preference(p_to,cfg);
		cps_tos->to_dataset=dysky.to_dataset;
		cps_tos->po_dataset=dysky.po_dataset;
		cps_tos->compute_skyline(cfg);
		tos.cached_preferences.push_back(p_to);
		tos.cached_skylines.push_back(cps_tos->skyline_result);
  	} while ( std::next_permutation(tos_values.begin(),tos_values.begin()+cfg->dyDim_val) );
  	cerr<<"--> Time for all TOS: "<< omp_get_wtime()-start_time << endl;

  	// skyline query answering

  	// input preference
	cerr << "Input preference: "<<endl;
	Preference p;
	p.add_vertices(cfg->dyDim_val);
	p.print_vertices();
	// for (id source=0;source<cfg->dyDim_val-1;source++){
	// 	unordered_set<id> v_dest;
	// 	v_dest.insert(source+1);
	// 	p.add_edges(source,v_dest);
	// }
	p.add_edges(0,{1,2,3});
	p.print_edges();
	Preference p_trans;
	p_trans.compute_transitive_closure(p);

	// skyline query answering by dySky
	cerr << "---query answering by dySky---"<<endl;

	vector<Order> preference_orders;
	unordered_map<id,unordered_set<id> > out_edges=p_trans.get_edges();
	for (auto it=out_edges.begin(); it!=out_edges.end();it++){
		for (auto it2=it->second.begin(); it2!=it->second.end(); it2++){
			preference_orders.push_back(Order((it->first),(*it2)));
		}
	}
	start_time2=omp_get_wtime();
	cerr << "size result: "<< dysky.compute_skyline(cfg, preference_orders).size()+dysky.always_sky.size()<<endl;
	cerr<<"--> Time for computing skyline by dySky : "<< omp_get_wtime()-start_time2 << endl;

  	// compute an issued query by ARG
	cerr << "---compute an issued query by ARG---" <<endl;
	p.is_subgraph(p);

	// compute an issued query by TOS
	cerr << "---compute an issued query by TOS---" <<endl;

  	// skyline query answering by CPS
	// CPS
	cerr << "=====CPS=====" <<endl;
	// start for preference decompositon
	start_time=omp_get_wtime();
	cerr << "---preference decompositon---"<<endl;
	Cps cps;
	start_time2=omp_get_wtime();
	cps.decompose_preference(p,cfg);
	cerr<<"--> Time for preference decompositon: "<< omp_get_wtime()-start_time2 << endl;
	//start for query answering
	cerr << "---query answering---"<<endl;
	start_time2=omp_get_wtime();
	cps.to_dataset=dysky.to_dataset;
	cps.po_dataset=dysky.po_dataset;
	cps.compute_skyline(cfg);
	cerr<<"--> Time for query answering: "<< omp_get_wtime()-start_time2 << endl;
	cerr<<"--> Time for all CPS: "<< omp_get_wtime()-start_time << endl;


}

