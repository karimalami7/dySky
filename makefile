CC=g++-7
CFLAGS=-O3 -fopenmp -finline -std=c++11  
MODULES=main.cpp generator/*.cpp BSkyTree/bskytree.cpp common.cpp

main: 
	
	$(CC) $(CFLAGS) $(MODULES) -o dySky

clean: 

	rm dySky