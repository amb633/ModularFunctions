#include "fullSolver.hpp"
#include "utility_functions.hpp"


void test_direct_solver()
{
	// function to test the full direct solver

    vector<vector<double>> A = {{-4.0, 1.0, 0.0, 0.0, 1.0}, {4.0, -4.0, 1.0, 0.0, 0.0}, { 0.0, 1.0, -4.0, 1.0, 0.0}, {0.0, 0.0, 1.0, -4.0, 1.0}, {1.0, 0.0, 0.0, 1.0, -4.0}};
    vector<double> b = {1.0, 0.0, 0.0, 0.0 ,0.0};

	vector<double> x;
	fullSolver( &A , &b , &x );

	printVector( &x );
}

void test_jacobi_iter_solver(){
    
    //function to test the iterative jacobi solver
    vector<vector<double>> A = {{-4.0, 1.0, 0.0, 0.0, 1.0}, {4.0, -4.0, 1.0, 0.0, 0.0}, { 0.0, 1.0, -4.0, 1.0, 0.0}, {0.0, 0.0, 1.0, -4.0, 1.0}, {1.0, 0.0, 0.0, 1.0, -4.0}};
    vector<double> b = {1.0, 0.0, 0.0, 0.0 ,0.0};
    
    
    vector<double> DInv;
    get_diagInv(&A, &DInv);
    vector<double> DInvB; // initializing the size of the vector
    vectorProduct(&DInv, &b, &DInvB); // calculating the D^-1 * b vector
    vector<double> x = DInvB;
    
    jacobiSolver( &A, &b, &x, 1.0e-7 );
    
    printVector( &x );
    
}

void test_SOR_iter_solver(){
    vector<vector<double>> A = {{-4.0, 1.0, 0.0, 0.0, 1.0}, {4.0, -4.0, 1.0, 0.0, 0.0}, { 0.0, 1.0, -4.0, 1.0, 0.0}, {0.0, 0.0, 1.0, -4.0, 1.0}, {1.0, 0.0, 0.0, 1.0, -4.0}};
    vector<double> b = {1.0, 0.0, 0.0, 0.0 ,0.0};
    
    vector<double> DInv;
    get_diagInv(&A, &DInv);
    vector<double> DInvB; // initializing the size of the vector
    vectorProduct(&DInv, &b, &DInvB); // calculating the D^-1 * b vector
    vector<double> x = DInvB;
    
    SORSolver(&A, &b, &x, 1.0e-7, 0.5);
    
    printVector(&x);
    
}
