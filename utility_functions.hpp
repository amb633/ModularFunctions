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

void addVectors( vector<double>* vector_1 , vector<double>* vector_2 , vector<double>* sum );
void multiplyVectors( vector<double>* vector_1 , vector<double>* vector_2 , vector<double>* product );

void scaleVector( double scalar , vector<double>* vector_1 , vector<double>* result );
void shiftVector( double scalar , vector<double>* vector_1 , vector<double>* result );
void expVector( vector<double>* input , vector<double>* result);

void vectorNorm( vector<double>* vector_1 , double& result );
void vectorNorm( vector<double>* vector_1 , vector<double>* vector_2 , double& result );

void get_random( double min , double max , double& random_number );