/*
 * main.cpp
 *
 *  Created on: December 19, 2018
 *      Author: karim
 */

#include "common/config.h"
#include "common/common.h"
#include "common/graph.h"
#include "common/preference.h"
#include "common/query.h"
#include "CPS/cps.h"
#include "dySky/dySky.h"
#include "dySky/dySky_h.h"
#include "dySky/dySky_v.h"
#include "dySky/dySky_v_chains.h"
#include "Ref/arg.h"
#include "OST/tos.h"
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
	time_t val=144301578;
    srand (val);

    printConfig(cfg);

    string const fileName("logs");
    ofstream myFile(fileName);
    uint64_t storage;
    bool selectedMethod[]={
    	false, //dysky_m
    	true, //dysky_v
    	true, //cps
    	false, //tos
    	false, //arg
    	false, //dysky_h
    	true, //dysky_v_chains
    };
  	//////////////////////////////////////////////////////////////////////////////
  	// Preprocessing
  	//

	cerr << "---PREPROCESSING---"<<endl<<endl;
	cout << "---PREPROCESSING---"<<endl<<endl;

	// start for pre processing
	cerr << "=====data and queries generation=====" <<endl;
	dySky dysky_m(cfg);
	// generate total order data
	dysky_m.generate_to_data(cfg);
	// generate partial order data
	dysky_m.generate_po_data(cfg);
	//
	dysky_m.print_dataset(cfg);
	//
	//********************************************
	// generate workload
	int compteur_query=0;
	vector<Query> workload(cfg->workload_size);
	for (int i=0; i<cfg->workload_size; i++){
		cerr << "compteur query: "<<i<<endl;
		cout << "compteur query: "<<i<<endl;
		Query q(cfg);
		workload[i]=q;
	}
	vector<Query> workload2(1000);
	//********************************

	//**********************************
	// recursive parallelization
	omp_set_nested(cfg->dyDim_size); 
	// 

	cfg->to_dataset=dysky_m.to_dataset;
	cfg->po_dataset=dysky_m.po_dataset;
	//
	vector<vector<Order>> all_orders(cfg->dyDim_size);
	generate_all_orders(cfg,all_orders);
	double start_time;
	double start_time2;

	//**********************************************
	// dySky_m
	if (selectedMethod[0]==true){
		cerr << "=====dySky=====" <<endl;	
		cout << "=====dySky=====" <<endl;	
		start_time=omp_get_wtime();
		// compute candidates 
		start_time2=omp_get_wtime();
		dysky_m.compute_candidates(cfg);
		cerr<<"--> Time for compute_candidates: "<< omp_get_wtime()-start_time2 << endl;
		// compute views
		start_time2=omp_get_wtime();
		storage=0;
		dysky_m.compute_views(cfg, all_orders, &storage);
		cerr<<"--> Total storage: "<< storage << endl;
		cerr<<"--> Time for compute_views: "<< omp_get_wtime()-start_time2 << endl;
		cerr<<"--> Time for all dySky: "<< omp_get_wtime()-start_time << endl;
	}
	//***********************************************

	cerr<<endl;
	cout<<endl;

	//************************************************
	// TOS: compute a skyline view wrt all possible total orders
	Tos tos(cfg);
	if(selectedMethod[3]){
		cerr << "=====TOS=====" <<endl;
		cout << "=====TOS=====" <<endl;
		start_time=omp_get_wtime();
	  	storage=0;
	  	tos.compute_views(cfg, &storage);
	  	cerr<<"--> Total storage: "<< storage << endl;
	  	cerr<<"--> Time for all TOS: "<< omp_get_wtime()-start_time << endl;
  	}
	//************************************************

	cerr<<endl;
	cout<<endl;

  	//***********************************************
	//ARG: compute skyline view wrt some partial order
	Arg arg;
	if(selectedMethod[4]==true){
		for (int i=0; i<workload2.size(); ++i){
			//cerr << compteur_query<<endl;
			Query q(cfg);
			workload2[i]=q;
			//compteur_query++;
		}
		cerr << "=====Arg=====" <<endl;
		cout << "=====Arg=====" <<endl;
		start_time=omp_get_wtime();
		storage=0;
		arg.compute_views(cfg, &storage, workload2);
		cerr<<"--> Total storage: "<< storage << endl;
		cerr<<"--> Time for all ARG: "<< omp_get_wtime()-start_time << endl;
	}	
  	//***********************************************

	//**********************************************
	// dysky_h
	dySky_h dysky_h(cfg);
	if (selectedMethod[5]==true){
		cerr << "=====dysky_h=====" <<endl;	
		cout << "=====dysky_h=====" <<endl;	
		dysky_h.to_dataset=dysky_m.to_dataset;
		dysky_h.po_dataset=dysky_m.po_dataset;
		start_time=omp_get_wtime();
		// compute candidates 
		start_time2=omp_get_wtime();
		dysky_h.compute_candidates(cfg);
		cerr<<"--> Time for compute_candidates: "<< omp_get_wtime()-start_time2 << endl;
		// compute views
		start_time2=omp_get_wtime();
		storage=0;
		dysky_h.compute_views(cfg, all_orders, &storage);
		dysky_h.views_selection(cfg, 16, workload);
		cerr<<"--> Total storage: "<< storage << endl;
		cerr<<"--> Time for compute_views: "<< omp_get_wtime()-start_time2 << endl;
		cerr<<"--> Time for all dySky: "<< omp_get_wtime()-start_time << endl;
	}
	//***********************************************

  	cerr<<endl;
  	cout<<endl;

  	/////////////////////////////////////////////////////////////////////////////
  	// skyline query answering
  	//

	string const nomFichier1("logs/performance-"+to_string(cfg->dataset_size)+"-"+to_string(cfg->statDim_size)+"-"+to_string(cfg->dyDim_size)+"-"+to_string(cfg->dyDim_val));
    ofstream monFlux1(nomFichier1.c_str());

  	// cerr << "---QUERY ANSWERING---"<<endl<<endl;
  	//cout << "---QUERY ANSWERING---"<<endl<<endl;

	//**************************************************************************************
	// print query answering performance

    monFlux1 << "Running Configuration:" <<endl;
    monFlux1 << "Dataset size: "<<cfg->dataset_size<<endl;
    monFlux1 << "Nb. Stat. Dimension: "<<cfg->statDim_size<<endl;
    monFlux1 << "Nb. Values Stat. Dim.: "<<cfg->statDim_val<<endl;
    monFlux1 << "Nb. Dyn. Dimension: : "<<cfg->dyDim_size<<endl;
    monFlux1 << "Nb. Values Dyn. Dim.: "<<cfg->dyDim_val<<endl;
    monFlux1 << "Workload size: "<<cfg->workload_size<<endl<<endl;	

	if (selectedMethod[0]) monFlux1 << "dysky_m"<< " : ";
	if (selectedMethod[1]) monFlux1 << "dysky_v"<< " : ";
	if (selectedMethod[2]) monFlux1 << "cps"<< " : ";
	if (selectedMethod[3]) monFlux1 << "tos"<< " : ";
	if (selectedMethod[4]) monFlux1 << "arg"<< " : ";
	if (selectedMethod[5]) monFlux1 << "dysky_h"<< " : ";
	if (selectedMethod[6]) monFlux1 << "dysky_v_chains"<< " : ";

	map<string, double> processing_time;
	processing_time["dysky_m"]=0;
	processing_time["dysky_v"]=0;
	processing_time["cps"]=0;
	processing_time["tos"]=0;
	processing_time["arg"]=0;
	processing_time["dysky_h"]=0;
	processing_time["dysky_v_chains"]=0;

	monFlux1 << endl;
	//*****************************************************************************************

	
  	// input preference
	//cerr << "Input preference: "<<endl;

	for (int q=0; q<cfg->workload_size; q++){

		cerr <<"Query n° "<<q<<endl<<endl;
		cout <<"Query n° "<<q<<endl<<endl;

		map<string, int> results; // to compare results of all methods
		int size_result; 

		// workload[q].generate_preference(cfg);
		// workload[q].graph_to_orderPairs(cfg);
		// workload[q].cross_orders_over_dimensions(cfg);
		////cout << "query preferences: "<<endl;
		// for (int i=0; i<cfg->dyDim_size; i++){
		// 	workload[q].preference[i].print_edges();
		// 	//cout <<endl;
		// }	

		// dySky_m: skyline query answering by dySky using materialized views
		if(selectedMethod[0]==true){
			cerr << "=====dySky: materialized views=====" <<endl;
			//cout << "=====dySky: materialized views=====" <<endl;
			start_time=omp_get_wtime();
			results["dysky_m"]=dysky_m.compute_skyline(cfg, workload[q].preference_orders_cross).size();
			processing_time["dysky_m"]=processing_time["dysky_m"]+(omp_get_wtime()-start_time);
			cerr << "--> Result size: "<< results["dysky_m"]<<endl;
			cerr << "--> Time: "<< processing_time["dysky_m"] << endl;
		}

		cerr <<endl;

		// dySky_v: skyline query answering by dySky using virtual views
		if(selectedMethod[1]==true){
			cerr << "=====dySky: on the fly with orders=====" <<endl;
			cout << "=====dySky: on the fly with orders=====" <<endl;
			dySky_v dysky_v(cfg);
			//dySky_v_chains dysky_v(cfg);
			dysky_v.to_dataset=dysky_m.to_dataset;
			dysky_v.po_dataset=dysky_m.po_dataset;
			start_time=omp_get_wtime();
			start_time2=omp_get_wtime();			
			dysky_v.compute_candidates(cfg);
			cerr<<"--> Time for compute_candidates: "<< omp_get_wtime()-start_time2 << endl;		
			start_time2=omp_get_wtime();
			results["dysky_v"]=dysky_v.compute_skyline(cfg, workload[q].preference_orders).size();
			//results["dysky_v"]=dysky_v.compute_skyline(cfg, workload[q].preference_chains).size();
			cerr<<"--> Time for compute_skyline: "<< omp_get_wtime()-start_time2 << endl;
			processing_time["dysky_v"]=processing_time["dysky_v"]+(omp_get_wtime()-start_time);
			cerr << "--> Result size: "<< results["dysky_v"]<<endl;
			cerr << "--> Time: "<< processing_time["dysky_v"] << endl;
		}

		// dySky_v: skyline query answering by dySky using virtual views
		if(selectedMethod[6]==true){
			cerr << "=====dySky: on the fly with chains=====" <<endl;
			cout << "=====dySky: on the fly with chains=====" <<endl;
			dySky_v_chains dysky_v(cfg);
			dysky_v.to_dataset=dysky_m.to_dataset;
			dysky_v.po_dataset=dysky_m.po_dataset;
			start_time=omp_get_wtime();
			start_time2=omp_get_wtime();			
			dysky_v.compute_candidates(cfg);
			cerr<<"--> Time for compute_candidates: "<< omp_get_wtime()-start_time2 << endl;		
			start_time2=omp_get_wtime();
			results["dySky_v_chains"]=dysky_v.compute_skyline(cfg, workload[q].preference_chains).size();
			cerr<<"--> Time for compute_skyline: "<< omp_get_wtime()-start_time2 << endl;
			processing_time["dySky_v_chains"]=processing_time["dySky_v_chains"]+(omp_get_wtime()-start_time);
			cerr << "--> Result size: "<< results["dySky_v_chains"]<<endl;
			cerr << "--> Time: "<< processing_time["dySky_v_chains"] << endl;
		}
		

	  	// skyline query answering by CPS
	  	if(selectedMethod[2]==true){
			cerr << "=====CPS=====" <<endl;
			cout << "=====CPS=====" <<endl;
			// start for preference decompositon
			//cerr << "---preference decompositon---"<<endl;
			Cps cps(cfg);
			start_time=omp_get_wtime();
			start_time2=omp_get_wtime();
			for (int i=0; i<workload[q].preference.size();i++){
				cps.decompose_preference(workload[q].preference[i],cfg,i);		
			}
			cps.encoding(cfg);
			cerr<<"--> Time for preference decompositon: "<< omp_get_wtime()-start_time2 << endl;
			start_time2=omp_get_wtime();
			results["cps"]=cps.compute_skyline(cfg, false);
			cerr<< "--> Result size: "<<results["cps"]<<endl;
			cerr<< "--> Time for query answering: "<< omp_get_wtime()-start_time2 << endl;
			processing_time["cps"]=processing_time["cps"]+(omp_get_wtime()-start_time);
			cerr<< "--> Time for all CPS: "<< processing_time["cps"] << endl;
		}

		cerr <<endl;

		// skyline query answering by TOS
		if(selectedMethod[3]==true){
			cerr << "=====TOS=====" <<endl;
			cout << "=====TOS=====" <<endl;
			start_time=omp_get_wtime();
			tos.paths=vector<vector<chain>>(cfg->dyDim_size);
			tos.define_paths(workload[q].preference,cfg);
			results["tos"]=tos.compute_skyline(cfg).size();
			processing_time["tos"]=processing_time["tos"]+(omp_get_wtime()-start_time);
			cerr<< "--> Result size: "<< results["tos"] <<endl;
			cerr << "--> Time: "<< processing_time["tos"] << endl;			
		}

	  	// skyline query answering by ARG
	  	if(selectedMethod[4]==true){
			cerr << "=====Arg=====" <<endl;
			cout << "=====Arg=====" <<endl;
			start_time=omp_get_wtime();
			arg.compute_skyline(cfg,workload[q]);
			results["arg"]=arg.skyline_result.size();
			cerr<< "--> Result size: "<< results["arg"] <<endl;
			processing_time["arg"]=processing_time["arg"]+(omp_get_wtime()-start_time);
			cerr << "--> Time: "<< processing_time["arg"] << endl;	
	  	}

		cerr <<endl;

		// dysky_h: skyline query answering by dySky using partially materialized views
		if(selectedMethod[5]==true){
			cerr << "=====dySky: hybrid=====" <<endl;
			cout << "=====dySky: hybrid=====" <<endl;
			start_time=omp_get_wtime();
			// compute hybrid skyline
			results["dysky_h"]=dysky_h.hybrid_compute_skyline(cfg, workload[q].preference_orders_cross).size();
			processing_time["dysky_h"]=processing_time["dysky_h"]+(omp_get_wtime()-start_time);
			cerr << "--> Result size: "<< results["dysky_h"]<<endl;
			cerr << "--> Time: "<< processing_time["dysky_h"] << endl;
		}

		cerr <<endl;
		cerr <<endl;
		cerr <<endl;


		//**************************************************************************************
		// check and print if results are different
		vector <int> v;
		for (auto it : results) v.push_back(it.second);
		if ( adjacent_find( v.begin(), v.end(), not_equal_to<int>() ) != v.end() )
		{
			cout << "Different answer for Query n° "<<q<<endl;
			for (int i=0; i<cfg->dyDim_size; i++){
				cout << "dim "<<i<<endl;
				workload[q].preference[i].print_edges();
			}	

			cout << "Results: "<<endl;
			if (selectedMethod[0]) cout << "dysky_m: "<< results["dysky_m"]<<endl;
			if (selectedMethod[1]) cout << "dysky_v: "<< results["dysky_v"]<<endl;
			if (selectedMethod[2]) cout << "cps: "<< results["cps"]<<endl;
			if (selectedMethod[3]) cout << "tos: "<< results["tos"]<<endl;
			if (selectedMethod[4]) cout << "arg: "<< results["arg"]<<endl;
			if (selectedMethod[5]) cout << "dysky_h: "<< results["dysky_h"]<<endl;
			if (selectedMethod[6]) cout << "dysky_v_chains: "<< results["dySky_v_chains"]<<endl;
		}
		else {
			cout <<"same results"<<endl;
		}
		cout<<endl;
		cerr<<endl;
	}

	if (selectedMethod[0]) monFlux1 << processing_time["dysky_m"]/cfg->workload_size<< " : ";
	if (selectedMethod[1]) monFlux1 << processing_time["dysky_v"]/cfg->workload_size<< " : ";
	if (selectedMethod[2]) monFlux1 << processing_time["cps"]/cfg->workload_size<< " : ";
	if (selectedMethod[3]) monFlux1 << processing_time["tos"]/cfg->workload_size<< " : ";
	if (selectedMethod[4]) monFlux1 << processing_time["arg"]/cfg->workload_size<< " : ";
	if (selectedMethod[5]) monFlux1 << processing_time["dysky_h"]/cfg->workload_size<< " : ";
	if (selectedMethod[6]) monFlux1 << processing_time["dySky_v_chains"]/cfg->workload_size<< " : ";

	monFlux1 << endl;
}

