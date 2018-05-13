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

void matrixProduct( vector<vector<double>>* matrix_1 , vector<vector<double>>* matrix_2 , vector<vector<double>>* result );
void vectorProduct( vector<vector<double>>* matrix_1 , vector<double>* vector_2 , vector<double>* result );

void addVectors( vector<double>* vector_1 , vector<double>* vector_2 , vector<double>* sum );
void multiplyVectors( vector<double>* vector_1 , vector<double>* vector_2 , vector<double>* product );
void scaleVector( double scalar , vector<double>* vector_1 , vector<double>* result );
void shiftVector( double scalar , vector<double>* vector_1 , vector<double>* result );

void vectorNorm( vector<double>* vector_1 );
void vectorNorm( vector<double>* vector_1 , vector<double>* vector_2 );

