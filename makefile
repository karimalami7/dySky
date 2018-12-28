CC=g++
CFLAGS=-O3 -std=c++11  
MODULES=main.cpp generator/*.cpp BSkyTree/bskytree.cpp common.cpp

main: 
	
	$(CC) $(CFLAGS) $(MODULES) -o dySky

clean: 

	rm dySky