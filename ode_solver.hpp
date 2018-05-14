#pragma once
#include "utility_functions.hpp"

#define FORWARD_EULER 0
#define RK34 1
#define HEUN_ONE 2

void ode_exponential_function( double time , vector<double>* input , vector<double>* output );
void test_forward_euler();
void test_heun_oneStep();

namespace ode_solvers
{
void ODE_SOLVER( void (*function)(double, vector<double>* , vector<double>* ),
	double time , double march , vector<double>* input , vector<double>* slope , int method );

void forward_euler ( void (*function)( double, vector<double>* , vector<double>* ) ,
	double time , double march , vector<double>* input , vector<double>* slope );

void heun_oneStep ( void (*function) ( double , vector<double>* , vector<double>* ) ,
	double time , double march , vector<double>* input , vector<double>* slope );
}