all:
	g++ h5extract.cpp -o h5extract -std=c++11 -lhdf5_cpp -lhdf5
