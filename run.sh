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

# rm cout
# rm cerr 

for t in 1
do
	for dataset_size in 100000
	do
		k=$(($dataset_size/100))
		for distrib in "ANTI"
		do
			for statDim in 6
			do
				for dyDim in 2
				do
					for value in 10
					do
						for queries in 3
						do
							echo "-n $dataset_size -r $distrib -s $statDim -k $k -d $dyDim -m $value -q $queries"
							OMP_NUM_THREADS=$t ./exec_dySky -r $distrib -n $dataset_size -s $statDim -k $k -d $dyDim -m $value -q $queries #> cout 2> cerr
							mv cout ./logs/cout-$dataset_size-$distrib-$statDim-$dyDim-$value
							mv cerr ./logs/cerr-$dataset_size-$distrib-$statDim-$dyDim-$value
						done
					done
				done
			done
		done
	done
done


# ./dySky -n 100000 -s 8 -k 100 -d 2 -m 10 -q 2
