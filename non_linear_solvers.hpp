#pragma once

#include "utility_functions.hpp"
#include "fullSolver.hpp"

void generic_nonlinear_function ( double time , vector<double>* input , vector<double>* output );

namespace secantMinimization
{
	void secant_delta( void (*function)(double , vector<double>* , vector<double>* ) , 
		double time , vector<double> initial_guess_1 , vector<double> initial_guess_2 , 
		vector<double>* parameter_solutions , vector<double>* performance_metrics , bool normalize = 0 );

	void recurrence_relation( double time , vector<double>* initial_guess_1 , vector<double>* initial_guess_2 , 
	vector<double>* output_1 , vector<double>* output_2 , vector<double>* result );

	void gradient_calculation( vector<double>* current_parameters , vector<double>* previous_parameters , 
	vector<double>* current_output , vector<double>* previous_output , vector<double>* gradients );

	void hessian_calculation( vector<double>* current_parameters , vector<double>* previous_parameters_1 , vector<double>* previous_parameters_2 , 
	vector<double>* current_output , vector<double>* previous_output_1 , vector<double>* previous_output_2 ,
	vector<vector<double>>* hessians );
}

void generic_nonlinear_function ( double time , vector<double>* input , vector<double>* output );
void exponential_function( double time , vector<double>* input , vector<double>* output );
void test_recurrence_relation();
void test_secant_gradient();
void test_secant_hessian();