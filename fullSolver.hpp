//
//  fullSolver2.hpp
//  PA3
//
//  Created by Ariana Bruno on 4/7/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#ifndef fullSolver_hpp
#define fullSolver_hpp

#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "utility_functions.hpp"

using namespace std;

void fullSolver( vector<vector<double>>* matrix , vector<double>* b , vector<double>* solution );
void conditionMatrix( vector<vector<double>>* matrix , vector<double>* b );

void calcAtomicVector( vector<vector<double>>* matrix, int row , vector<vector<double>>* m);

void backwardSubstitution( vector<vector<double>>* U , vector<double>* y , vector<double>* x );

void test_direct_solver();

void get_diagInv( vector<vector<double>>* input, vector<double>* diag );

// decompose AF matrix into diagonal elements (stored in DF) and non-diagonal elements (stored in LUF)
void decomposeMatrix( vector<double>* DF , vector< vector<double>>* LUF , vector< vector<double>>* AF );

// jacobi iteration for jacobi solver function
void jacobiIter( vector<double>* X , vector<double>* DF , vector< vector<double>>* LUF , vector<double>* BF );

//jacobi solver, input x as initial guess, output result in updated x vector, pass in tolerance to achieve convergence
void jacobiSolver( vector<vector<double>>* A, vector<double>* b, vector<double>* x, double tol );

void test_jacobi_iter_solver();

// the SOR iter for updating x vector
void SORIter( vector<vector<double>>* A, vector<double>* b, vector<double>* x, double omega);

// the SOR iter for updating x vector
void SORIter( vector<vector<double>>* A, vector<double>* b, vector<double>* x, double omega);

// the full SOR solver, passing in x with initial guesses and given the tolerance and omega
void SORSolver( vector<vector<double>>* A, vector<double>* b , vector<double>* x, double tol, double omega);

void test_SOR_iter_solver();


#endif /* fullSolver_hpp */
