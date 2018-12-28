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
	cfg->statDim_val=100;
	cfg->dyDim_size=1;
	cfg->dyDim_val=8;
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

	// if (cfg->dataset_size.empty() || cfg->statDim_size.empty()) {
	// 	printUsage();
	// 	return 1;
	// }

	cout << "all ok"<<endl;

	dySky dysky;
	// generate total order data
	dysky.generate_to_data(cfg);
	// generate partial order data
	dysky.generate_po_data(cfg);
	// input preference
	Preference p;
	p.add_vertices(cfg->dyDim_val);
	p.print_vertices();
	// transform partial order dimension to k total order dimension
	//dysky.transform_po-to()
	// compute skyline by BSkyTree
	dysky.compute_skyline(cfg);

}
