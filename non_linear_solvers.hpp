#pragma once

#include "utility_functions.hpp"
#include "fullSolver.hpp"


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

namespace quasiNewtonMinimization
{
	void quasiNewton_delta( void (*function)( double , vector<double>* , vector<double>* ) ,
	double time , vector<double> initial_guess_1 , double pert , 
	vector<double>* parameter_solutions , vector<double>* performance_metrics , vector<double>* desired_result , bool normalize = 0);

	void gradient_calculation( vector<double>* current_parameters , double pert ,
	vector<double>* current_output , vector<double>*perturbed_output , vector<double>* gradients );

	void hessian_calculation( vector<double>* current_parameters , double pert , 
	vector<double>* current_output , vector<double>* perturbed_output , vector<vector<double>>* perturbed_output_cross ,
	vector<vector<double>>* hessians );

}

void tSearch( void (*function)( double , vector<double>* , vector<double>* ), vector<double>* delta, vector<double>* t_delta, vector<double>* x, double t_min, double t_max, vector<double>* true_result);

void generic_nonlinear_function ( double time , vector<double>* input , vector<double>* output );
void exponential_function( double time , vector<double>* input , vector<double>* output );
void test_recurrence_relation();
void test_secant_gradient();
void test_secant_hessian();

void test_quasiNewton_gradient();
void test_quasiNewton_hessian();
