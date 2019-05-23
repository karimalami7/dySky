/*
 * main.cpp
 *
 *  Created on: December 19, 2018
 *      Author: karim
 */


#include "common.h"
#include "config.h"
#include "graph.h"
#include "preference.h"
#include "query.h"
#include "cps.h"
#include "dySky.h"
//#include "arg.h"
//#include "tos.h"
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
	printf("Example: ");
	printf("./dySky -n 100000 -s 8 -k 100 -d 2 -m 10 \n\n");
}


int main(int argc, char** argv) {

	Config *cfg = new Config;
	cfg->dataset_size=1000;
	cfg->statDim_size=4;
	cfg->statDim_val=100;
	cfg->dyDim_size=1;
	cfg->dyDim_val=3;
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
		default:
			if (isprint(optopt))
				fprintf( stderr, "Unknown option `-%c'.\n", optopt);
			printUsage();
			return 1;
		}
	}
	time_t val=144630157;
    srand (val);

  	//////////////////////////////////////////////////////////////////////////////
  	// Preprocessing
  	//

	cerr << "---PREPROCESSING---"<<endl<<endl;

	//**********************************************
	// dySky
	cerr << "=====dySky=====" <<endl;
	// start for pre processing
	
	dySky dysky(cfg);
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
	dysky.compute_views2(cfg);
	cerr<<"--> Time for compute_views: "<< omp_get_wtime()-start_time2 << endl;
	// start_time2=omp_get_wtime();
	// dysky.compute_views(cfg);
	// cerr<<"--> Time for compute_views: "<< omp_get_wtime()-start_time2 << endl;
	dysky.print_dataset(cfg);
	cerr<<"--> Time for all dySky: "<< omp_get_wtime()-start_time << endl;
	//***********************************************
	
	cerr<<endl;

	//***********************************************
	//ARG: compute skyline view wrt some partial order

	// cerr << "=====Arg=====" <<endl;
	// start_time=omp_get_wtime();
	// Arg arg;
	// arg.to_dataset=dysky.to_dataset;
	// arg.po_dataset=dysky.po_dataset;
	// arg.compute_views(cfg);

	// cerr<<"--> Time for all ARG: "<< omp_get_wtime()-start_time << endl;
	// //************************************************

	// cerr<<endl;
	
	// //************************************************
	// // TOS: compute a skyline view wrt all possible total orders

	// cerr << "=====TOS=====" <<endl;
	// start_time=omp_get_wtime();

 //  	Tos tos;
 //  	chain tos_values(cfg->dyDim_val);
 //  	for (int i=0; i<cfg->dyDim_val; i++){tos_values[i]=i;}
 //  	do {
 //  	  	Preference p_to;
	// 	p_to.add_vertices(cfg->dyDim_val);
	// 	for (id source=0;source<cfg->dyDim_val-1;source++){
	// 		unordered_set<id> v_dest;
	// 		v_dest.insert(tos_values[source+1]);
	// 		p_to.add_edges(tos_values[source],v_dest);
	// 	}
	// 	Cps *cps_tos = new Cps();
	// 	cps_tos->decompose_preference(p_to,cfg);
	// 	cps_tos->to_dataset=dysky.to_dataset;
	// 	cps_tos->po_dataset=dysky.po_dataset;
	// 	cps_tos->compute_skyline(cfg);
	// 	tos.cache[tos_values]=cps_tos->skyline_result;

 //  	} while ( std::next_permutation(tos_values.begin(),tos_values.begin()+cfg->dyDim_val) );
 //  	cerr<<"--> Time for all TOS: "<< omp_get_wtime()-start_time << endl;
  	//***********************************************
  	
  	cerr<<endl;
  	cout<<endl;

  	/////////////////////////////////////////////////////////////////////////////
  	// skyline query answering
  	//


  	cerr << "---QUERY ANSWERING---"<<endl<<endl;
  	// input preference
	//cerr << "Input preference: "<<endl;
	Query q;
	q.generate_preference(cfg);
	q.graph_to_orderPairs(cfg);
	q.cross_orders_over_dimensions(cfg);
	cout << "query preferences: "<<endl;
	for (int i=0; i<cfg->dyDim_size; i++){
		q.preference[i].print_edges();
	}	

	// skyline query answering by dySky using materialized 
	cerr << "=====dySky: materialized views=====" <<endl;
	cout << "=====dySky: materialized views=====" <<endl;
	start_time2=omp_get_wtime();
	int size_result; 
	size_result=dysky.compute_skyline2(cfg, q.preference_orders_cross).size();
	cerr << "--> Result size: "<< size_result<<endl;
	cerr << "--> Time: "<< omp_get_wtime()-start_time2 << endl;
	// start_time2=omp_get_wtime();
	// size_result=dysky.compute_skyline(cfg, q.preference_orders, q).size();
	// cerr << "--> Result size: "<< size_result<<endl;
	// cerr << "--> Time: "<< omp_get_wtime()-start_time2 << endl;
	// cerr <<endl;

	cerr << "=====dySky: virtual views=====" <<endl;
	cout << "=====dySky: virtual views=====" <<endl;
	// delete materialized views and compute only the views need for the issued query
	// dysky.skyline_view=vector<unordered_map<Order,vector<id>, pairhash>>(cfg->dyDim_size);
	// start_time=omp_get_wtime();
	// dysky.compute_views(cfg, preference_orders);
	// size_result=dysky.compute_skyline(cfg, preference_orders, q).size();
	// cerr << "--> Result size: "<< size_result<<endl;
	// cerr << "--> Time: "<< omp_get_wtime()-start_time << endl;

	// cerr <<endl;

  	// skyline query answering by CPS
	cerr << "=====CPS=====" <<endl;
	cout << "=====CPS=====" <<endl;
	// start for preference decompositon
	start_time=omp_get_wtime();
	//cerr << "---preference decompositon---"<<endl;
	Cps cps(cfg);
	start_time2=omp_get_wtime();
	for (int i=0; i<q.preference.size();i++){
		cps.decompose_preference(q.preference[i],cfg,i);		
	}
	cps.encoding(cfg);
	cerr<<"--> Time for preference decompositon: "<< omp_get_wtime()-start_time2 << endl;
	//start for query answering
	//cerr << "---query answering---"<<endl;
	start_time2=omp_get_wtime();
	cps.to_dataset=dysky.to_dataset;
	cps.po_dataset=dysky.po_dataset;
	for (int i=0;i<cfg->dyDim_size;i++){
		cps.compute_skyline_perDimension(cfg, i);	
	}
	size_result= cps.compute_skyline(cfg);
	cerr<< "--> Result size: "<<size_result<<endl;
	cerr<< "--> Time for query answering: "<< omp_get_wtime()-start_time2 << endl;
	cerr<< "--> Time for all CPS: "<< omp_get_wtime()-start_time << endl;

	cerr <<endl;

  	// compute an issued query by ARG
	// cerr << "=====Arg=====" <<endl;
	// start_time=omp_get_wtime();
	// arg.compute_skyline(cfg,q);
	// cerr<< "--> Result size: "<<arg.skyline_result.size() <<endl;
	// cerr << "--> Time: "<< omp_get_wtime()-start_time << endl;	
	// cerr <<endl;

	// // compute an issued query by TOS
	// cerr << "=====TOS=====" <<endl;
	// start_time=omp_get_wtime();
	// Cps cps_for_tos;
	// cps_for_tos.decompose_preference(q.preference,cfg);
	// cerr << "--> number of decomposed chains: " << cps_for_tos.chains.size() <<endl;
	// cerr<< "--> Result size: "<<tos.compute_skyline(cps_for_tos.chains, cfg).size() <<endl;
	// cerr << "--> Time: "<< omp_get_wtime()-start_time << endl;

	//**************************************************************************************

}

