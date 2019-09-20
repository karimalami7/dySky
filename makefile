CC=g++-7
CFLAGS=-O3 -fopenmp -finline -std=c++11  
MODULES=main.cpp generator/*.cpp BSkyTree/bskytree.cpp common/common.cpp

main: 
	
	$(CC) $(CFLAGS) $(MODULES) -o exec_dySky

clean: 

	rm dysky