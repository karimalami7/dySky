#/bin/bash
###############################
#
#	dySky - computes Skyline over datasets with dynamic partial order
#	USAGE: ./dySky [Option Argument]
#	-n: size of the dataset
#	-s: number of static dimensions
#	-k: distinct values of static dimensions
#	-d: number of dynamic dimensions
#	-m: distinct values of dynamic dimensions
#	-q: number of queries to answer
#	Example:
#	./dySky -n 100000 -s 8 -k 100 -d 2 -m 10 -q 2 \n\n

rm cout
rm cerr 

for dataset_size in 10000 100000 1000000  
do
	for statDim in 6 8 10
	do
		for dyDim in 1 2
		do
			for value in 4 6 8
			do
				./dySky -n $dataset_size -s $statDim -k 100 -d $dyDim -m $value -q 5 > cout
				mv cout ./logs/cout-$dataset_size-$statDim-$dyDim-$value
			done
		done
	done
done



# ./dySky -n 100000 -s 8 -k 100 -d 2 -m 10 -q 2