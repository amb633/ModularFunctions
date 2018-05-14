#pragma once
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

void printMatrix( vector<vector<double>>* matrix );
void printVector( vector<double>* vector_1 );

void identityMatrix( int rank , vector<vector<double>>* I );
void zeroMatrix( int rank , vector<vector<double>>* z );
void zeroVector( int rank , vector<double>* v );
void copyMatrix( vector<vector<double>>* original , vector<vector<double>>* copy );

void matrixProduct( vector<vector<double>>* matrix_1 , vector<vector<double>>* matrix_2 , vector<vector<double>>* result );
void vectorProduct( vector<vector<double>>* matrix_1 , vector<double>* vector_1 , vector<double>* result );
void vectorProduct( vector<double>* diagMat , vector<double>* vec, vector<double>* result );

void addVectors( vector<double>* vector_1 , vector<double>* vector_2 , vector<double>* sum );
void vectorSum(vector<double>* A, vector<double>* B, vector<double>* result); //inplace sum (result can be vector A or B);
void multiplyVectors( vector<double>* vector_1 , vector<double>* vector_2 , vector<double>* product );

void scaleVector( double scalar , vector<double>* vector_1 , vector<double>* result );
void shiftVector( double scalar , vector<double>* vector_1 , vector<double>* result );
void expVector( vector<double>* input , vector<double>* result);

void vectorNorm( vector<double>* vector_1 , double& result );
void vectorNorm( vector<double>* vector_1 , vector<double>* vector_2 , double& result );

void get_random( double min , double max , double& random_number );

    
// returns element at i,j in matrix AF / element i in vector VF
double retrieveElement( vector<vector<double>>* AF, int rowInd, int colInd);
double retrieveElement( vector<double>* VF , int rowInd );

// print the matrix AF / vector VF
//void print_full_mat( vector< vector<double>>* AF );
//void print_full_vec( vector<double>* VF );

// set element at i,j in AF to newValue
void changeElement( vector< vector<double>>* AF , int rowInd , int colInd , double newValue );

// create a copy of matrix AF and store it in CF
//int copyMatrix( vector< vector<double>>* CF , vector< vector<double>>* AF );

// multiply all (non-zero) elements in AF with scale
void scaleMatrix( vector< vector<double>>* Matrix , double scale );

// swaps row i,j in matrix AF
void rowPermute( vector< vector <double>>* AF , int i , int j );

// returns row_j = row_j + a * row_i
void rowScale(vector< vector<double>>* AF, int i, int j, double a );

// calculates the matrix product of AF and VF and stores it in result
//int productAx( vector<double>* result , vector< vector<double>>* AF , vector<double>* VF );

// calculate norm of v - AX
//int calculateNorm( double& norm , vector<double>* v , vector<double>* Ax );



