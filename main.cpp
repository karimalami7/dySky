/*
 * main.cpp
 *
 *  Created on: December 19, 2018
 *      Author: karim
 */

#include "config.h"
#include "common.h"
#include "graph.h"
#include "preference.h"
#include "query.h"
#include "cps.h"
#include "dySky.h"
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
	printf( " -q: number of queries to answer\n" );	
	printf("Example: ");
	printf("./dySky -n 100000 -s 8 -k 100 -d 2 -m 10 -q 2 \n\n");
}


int main(int argc, char** argv) {

	Config *cfg = new Config;
	cfg->dataset_size=1000;
	cfg->statDim_size=4;
	cfg->statDim_val=100;
	cfg->dyDim_size=1;
	cfg->dyDim_val=3;
	cfg->workload_size=2;
	cfg->verbose=false;

	int c = 0;
	opterr = 0;

	while ((c = getopt(argc, argv, "n:s:k:d:m:q:")) != -1) {
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
		case 'q':
			cfg->workload_size = atoi(optarg);
			break;
		default:
			if (isprint(optopt))
				fprintf( stderr, "Unknown option `-%c'.\n", optopt);
			printUsage();
			return 1;
		}
	}
	time_t val=14430157;
    srand (val);

    printConfig(cfg);

    string const fileName("logs");
    ofstream myFile(fileName);

  	//////////////////////////////////////////////////////////////////////////////
  	// Preprocessing
  	//

	cerr << "---PREPROCESSING---"<<endl<<endl;
	//cout << "---PREPROCESSING---"<<endl<<endl;

	//**********************************************
	// dySky
	cerr << "=====dySky=====" <<endl;
	// start for pre processing
	
	dySky dysky_m(cfg);
	// generate total order data
	dysky_m.generate_to_data(cfg);
	// generate partial order data
	dysky_m.generate_po_data(cfg);
	
	vector<vector<Order>> all_orders(cfg->dyDim_size);
	generate_all_orders(cfg,all_orders);
		
	double start_time=omp_get_wtime();
	double start_time2=omp_get_wtime();
	dysky_m.compute_always_skyline(cfg);
	cerr<<"--> Time for compute_always_skyline: "<< omp_get_wtime()-start_time2 << endl;
	// compute candidates 
	start_time2=omp_get_wtime();
	dysky_m.compute_candidates(cfg);
	cerr<<"--> Time for compute_candidates: "<< omp_get_wtime()-start_time2 << endl;
	// compute views
	start_time2=omp_get_wtime();
	dysky_m.compute_views(cfg, all_orders);
	cerr<<"--> Time for compute_views: "<< omp_get_wtime()-start_time2 << endl;
	dysky_m.print_dataset(cfg);
	cerr<<"--> Time for all dySky: "<< omp_get_wtime()-start_time << endl;
	//***********************************************
	cout << "dysky completed" <<endl;
	cerr<<endl;

	//***********************************************
	//ARG: compute skyline view wrt some partial order

	cerr << "=====Arg=====" <<endl;
	start_time=omp_get_wtime();
	Arg arg;
	arg.to_dataset=dysky_m.to_dataset;
	arg.po_dataset=dysky_m.po_dataset;
	arg.compute_views(cfg);

	cerr<<"--> Time for all ARG: "<< omp_get_wtime()-start_time << endl;
	//************************************************

	cerr<<endl;
	
	//************************************************
	// TOS: compute a skyline view wrt all possible total orders

	cerr << "=====TOS=====" <<endl;
	start_time=omp_get_wtime();

  	Tos tos;
  	tos.to_dataset=dysky_m.to_dataset;
  	tos.po_dataset=dysky_m.po_dataset;
  	tos.compute_views(cfg);

  	cerr<<"--> Time for all TOS: "<< omp_get_wtime()-start_time << endl;
  	//***********************************************
  	
  	cerr<<endl;
  	//cout<<endl;

  	/////////////////////////////////////////////////////////////////////////////
  	// skyline query answering
  	//


  	cerr << "---QUERY ANSWERING---"<<endl<<endl;
  	//cout << "---QUERY ANSWERING---"<<endl<<endl;
  	// input preference
	//cerr << "Input preference: "<<endl;
	vector<Query> workload(cfg->workload_size);
	for (int q=0; q<cfg->workload_size; q++){

		cerr <<"Query n° "<<q<<endl<<endl;
		//cout <<"Query n° "<<q<<endl<<endl;

		vector <int> results; // to compare results of all methods

		workload[q].generate_preference(cfg);
		workload[q].graph_to_orderPairs(cfg);
		workload[q].cross_orders_over_dimensions(cfg);
		//cout << "query preferences: "<<endl;
		// for (int i=0; i<cfg->dyDim_size; i++){
		// 	workload[q].preference[i].print_edges();
		// }	

		// skyline query answering by dySky using materialized 
		cerr << "=====dySky: materialized views=====" <<endl;
		//cout << "=====dySky: materialized views=====" <<endl;
		start_time2=omp_get_wtime();
		int size_result; 
		size_result=dysky_m.compute_skyline(cfg, workload[q].preference_orders_cross).size();
		results.push_back(size_result);
		cerr << "--> Result size: "<< size_result<<endl;
		cerr << "--> Time: "<< omp_get_wtime()-start_time2 << endl;

		cerr << "=====dySky: virtual views=====" <<endl;
		//cout << "=====dySky: virtual views=====" <<endl;
		// delete materialized views and compute only the views need for the issued query
		dySky dysky_v(cfg);
		dysky_v.to_dataset=dysky_m.to_dataset;
		dysky_v.po_dataset=dysky_m.po_dataset;
		start_time=omp_get_wtime();
		dysky_v.compute_always_skyline(cfg);
		dysky_v.compute_candidates(cfg);
		dysky_v.compute_views(cfg, workload[q].preference_orders);
		size_result=dysky_v.compute_skyline(cfg, workload[q].preference_orders_cross).size();
		results.push_back(size_result);
		cerr << "--> Result size: "<< size_result<<endl;
		cerr << "--> Time: "<< omp_get_wtime()-start_time << endl;

		cerr <<endl;

	  	// skyline query answering by CPS
		cerr << "=====CPS=====" <<endl;
		//cout << "=====CPS=====" <<endl;
		// start for preference decompositon
		start_time=omp_get_wtime();
		//cerr << "---preference decompositon---"<<endl;
		Cps cps(cfg);
		start_time2=omp_get_wtime();
		for (int i=0; i<workload[q].preference.size();i++){
			cps.decompose_preference(workload[q].preference[i],cfg,i);		
		}
		cps.encoding(cfg);
		cerr<<"--> Time for preference decompositon: "<< omp_get_wtime()-start_time2 << endl;
		//start for query answering
		//cerr << "---query answering---"<<endl;
		start_time2=omp_get_wtime();
		cps.to_dataset=dysky_m.to_dataset;
		cps.po_dataset=dysky_m.po_dataset;
		// for (int i=0;i<cfg->dyDim_size;i++){
		// 	cps.compute_skyline_perDimension(cfg, i);	
		// }
		size_result= cps.compute_skyline(cfg);
		results.push_back(size_result);
		cerr<< "--> Result size: "<<size_result<<endl;
		cerr<< "--> Time for query answering: "<< omp_get_wtime()-start_time2 << endl;
		cerr<< "--> Time for all CPS: "<< omp_get_wtime()-start_time << endl;

		cerr <<endl;

	  	// compute an issued query by ARG
		cerr << "=====Arg=====" <<endl;
		start_time=omp_get_wtime();
		arg.compute_skyline(cfg,workload[q]);
		size_result=arg.skyline_result.size();
		results.push_back(size_result);
		cerr<< "--> Result size: "<< size_result<<endl;
		cerr << "--> Time: "<< omp_get_wtime()-start_time << endl;	
		cerr <<endl;

		// compute an issued query by TOS
		cerr << "=====TOS=====" <<endl;
		//cout << "=====TOS=====" <<endl;
		start_time=omp_get_wtime();
		Cps cps_for_tos(cfg);
		for (int i=0;i<cfg->dyDim_size;i++){
			cps_for_tos.decompose_preference(workload[q].preference[i],cfg,i);
		}
		//cerr << "--> number of decomposed chains: " << cps_for_tos.chains.size() <<endl;
		size_result=tos.compute_skyline(cps_for_tos.chains, cfg).size();
		results.push_back(size_result);
		cerr<< "--> Result size: "<< size_result <<endl;
		cerr << "--> Time: "<< omp_get_wtime()-start_time << endl;

		//**************************************************************************************

		if ( adjacent_find( results.begin(), results.end(), not_equal_to<int>() ) != results.end() )
		{
			cout << "Query n° "<<q<<endl;
			for (int i=0; i<cfg->dyDim_size; i++){
				workload[q].preference[i].print_edges();
			}	

			cout << "Results: "<<endl;
			cout << "dysky_m: "<< results[0]<<endl;
			cout << "dysky_v: "<< results[1]<<endl;
			cout << "cps: "<< results[2]<<endl;
			cout << "arg: "<< results[3]<<endl;
			cout << "tos: "<< results[4]<<endl;
		}
		cout<<endl;
		cerr<<endl;
	}
}

