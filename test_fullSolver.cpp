#include "fullSolver.hpp"
#include "utility_functions.hpp"


void test_direct_solver()
{
	// function to test the full direct solver

	vector<vector<double>> A;
	A.push_back( { 4 , 4 , 2 , 3 });
	A.push_back( { 2 , 3 , 3 , 1 });
	A.push_back( { 4 , 5 , 3 , 2 });
	A.push_back( { 1 , 3 , 2 , 4 });

	vector<double> b = { -1 , 1 , 4 , -3 };

	vector<double> x;
	fullSolver( &A , &b , &x );

	printVector( &x );
}