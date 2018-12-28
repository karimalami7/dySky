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

typedef struct Config {
	int dataset_size;
	int statDim_size;
	int statDim_val;
	int dyDim_size;
	int dyDim_val;
	int num_threads;
	bool verbose;
} Config;


#endif /* CONFIG_H_ */
