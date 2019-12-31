/*
 * config.h
 *
 *  Created on: December 19, 2018
 *      Author: karim
 */

#ifndef CONFIG_H_
#define CONFIG_H_
#include <string>
#include <vector>

using namespace std;

typedef struct Config {
	int dataset_size;
	int statDim_size;
	int statDim_val;
	int dyDim_size;
	int dyDim_val;
	int workload_size;
	bool verbose;
	vector<int*> to_dataset;
	vector<vector<int>> po_dataset;
	char * distrib;
} Config;


#endif /* CONFIG_H_ */
