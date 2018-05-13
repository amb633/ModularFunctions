#include "utility_functions.hpp"

void printMatrix( vector<vector<double>>* matrix )
{
	// assuming square matrix
	int rank = (*matrix).size();
	for ( int i = 0 ; i < rank ; i++ ){
		for ( int j = 0 ; j < rank ; j++ ){
			cout << ((*matrix)[i])[j] << "	";
		}
		cout << endl;
	}

}

void printVector( vector<double>* vector_1 )
{
	int rank = (*vector_1).size();
	for ( int i = 0 ; i < rank ; i++ ){
		cout << (*vector_1)[i] << "   ";
	}
	cout << endl;
}


void identityMatrix( int rank , vector<vector<double>>* I )
{
	// function to create identity matrix of size rank

	for ( int c = 0 ; c < rank ; c++ ){
		vector<double> row;
		for ( int r = 0 ; r < rank ; r++ ){
			if ( r == c ) row.push_back(1.0);
			else row.push_back(0.0);
		}
		(*I).push_back(row);
	}

}

void zeroMatrix( int rank , vector<vector<double>>* z )
{
	// function to create zero matrix of size rank
	for ( int c = 0 ; c < rank ; c++ ){
		vector<double> row;
		for ( int r = 0 ; r < rank ; r++ ){
			row.push_back(0.0);
		}
		(*z).push_back(row);
	}
}