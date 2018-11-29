all: kmeans

kmeans: temp/main.o temp/Menu.o
	g++ -g -std=c++11 temp/main.o temp/Menu.o
	mv a.out kmeans

temp/Menu.o: src/Menu.cpp include/Menu.h
	g++ -c -g -std=c++11 src/Menu.cpp
	mv Menu.o temp/Menu.o

temp/main.o: src/main.cpp include/Menu.h
	g++ -c -g -std=c++11 src/main.cpp
	mv main.o temp/main.o

clean:
	rm -f temp/*
	rm -f kmeans