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
	cfg->dyDim_val=10;
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

	// start for pre processing

	dySky dysky;
	// generate total order data
	dysky.generate_to_data(cfg);
	// generate partial order data
	dysky.generate_po_data(cfg);
	// compute skyline by BSkyTree on static order dimensions (Always Sky)
	double start_time=omp_get_wtime();
	dysky.compute_always_skyline(cfg);
	cout<<"--> Time for compute_always_skyline: "<< omp_get_wtime()-start_time << endl;
	// compute candidates 
	start_time=omp_get_wtime();
	dysky.compute_candidates(cfg);
	cout<<"--> Time for compute_candidates: "<< omp_get_wtime()-start_time << endl;
	// compute views
	start_time=omp_get_wtime();
	dysky.compute_views(cfg);
	cout<<"--> Time for compute_views: "<< omp_get_wtime()-start_time << endl;
	dysky.print_dataset(cfg);

	// end of pre processing

	// start of skyline query answering

	// input preference
	Preference p;
	p.add_vertices(cfg->dyDim_val);
	p.print_vertices();
	for (id source=0;source<cfg->dyDim_val-1;source++){
		vector<id> v_dest;
		for (int dest=source+1;dest<cfg->dyDim_val;dest++){
			v_dest.push_back(dest);
		}
		p.add_edges(source,v_dest);
	}
	p.print_edges();

	// skyline query answering by dySky

	vector<Order> preference_orders;
	unordered_map<id,vector<id> > out_edges=p.get_edges();
	for (auto it=out_edges.begin(); it!=out_edges.end();it++){
		for (auto it2=it->second.begin(); it2!=it->second.end(); it2++){
			preference_orders.push_back(Order((it->first),(*it2)));
		}
	}
	start_time=omp_get_wtime();
	cout << "size result: "<< dysky.compute_skyline(cfg, preference_orders).size()<<endl;
	cout<<"--> Time for compute_skyline: "<< omp_get_wtime()-start_time << endl;

	// skyline query answering by CPS

	Cps cps;


}
