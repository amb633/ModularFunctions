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

using namespace std;

void fullSolver( vector<vector<double>>* matrix , vector<double>* b , vector<double>* solution );
void conditionMatrix( vector<vector<double>>* matrix , vector<double>* b );

void calcAtomicVector( vector<vector<double>>* matrix, int row , vector<vector<double>>* m);

void backwardSubstitution( vector<vector<double>>* U , vector<double>* y , vector<double>* x );

void test_direct_solver();


#endif /* fullSolver_hpp */
