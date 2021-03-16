CXXFLAGS = -O3 -fPIC
CXX = g++
#CXXFLAGS = -g
default: FastMatSparse.o

DMS.o: DMS.cpp
	$(CXX) $(CXXFLAGS) -I. -c DMS.cpp -o DMS.o


DegMatSparse.o: DegMatSparse.cpp
	$(CXX) $(CXXFLAGS) -I. -c DegMatSparse.cpp -o DegMatSparse.o

FastMatSparse.o: FastMatSparse.cpp
	$(CXX) $(CXXFLAGS) -I. -c FastMatSparse.cpp -o FastMatSparse.o

MapMatSparse.o: MapMatSparse.cpp
	$(CXX) $(CXXFLAGS) -I. -c MapMatSparse.cpp -o MapMatSparse.o

main.o: main.cpp DegMatSparse.o MapMatSparse.o
	$(CXX) $(CXXFLAGS) -I. -c main.cpp -o main.o

read.o: read.cpp FastMatSparse.o DegMatSparse.o MapMatSparse.o isearch.o DMS.o
	$(CXX) $(CXXFLAGS) -I. -c read.cpp -o read.o

isearch.o: isearch.cpp
	$(CXX) $(CXXFLAGS) -I. -c isearch.cpp -o isearch.o

matgen.o: matgen.cpp DegMatSparse.o MapMatSparse.o
	$(CXX) $(CXXFLAGS) -I. -c matgen.cpp -o matgen.o

app: main.o
	$(CXX) $(CXXFLAGS) -I. isearch.o DegMatSparse.o MapMatSparse.o main.o -o app

read: read.o
	$(CXX) $(CXXFLAGS) -I. isearch.o FastMatSparse.o DegMatSparse.o MapMatSparse.o read.o -o read

matgen: matgen.o
	$(CXX) $(CXXFLAGS) -I. isearch.o DegMatSparse.o MapMatSparse.o matgen.o -o matgen

clean:
	rm -f *o read app
