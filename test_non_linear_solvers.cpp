#include "non_linear_solvers.hpp"

void generic_nonlinear_function ( double time , vector<double>* input , vector<double>* output )
{
	// evaluate output = 2*x^2 + 4*x*y + 2*y^2
	// x = input[0]
	// y = input[1]
	// time is redundant in this case
	double x , y , temp;
	x = (*input)[0];
	y = (*input)[1];
	temp = 2.0*x*x + 4.0*x*y + 2.0*y*y;
	(*output).push_back(temp);
	return;
}

void exponential_function( double time , vector<double>* input , vector<double>* output )
{
	double x , temp;
	x = (*input)[0];
	temp = exp( 100.0*x ) - 1.0;
	(*output).push_back(temp);
	return;
}

void test_recurrence_relation()
{
	void (*function)( double , vector<double>* , vector<double>* ) = generic_nonlinear_function;
	vector<double> guess_1 = { 1.0 , 1.0 };
	vector<double> guess_2 = { 2.0 , 3.0 };
	vector<double> output_1 , output_2;
	generic_nonlinear_function( 0.0 , &guess_1 , &output_1 );
	generic_nonlinear_function( 0.0 , &guess_2 , &output_2 );

	vector<double> guess_3 ;
	secantMinimization::recurrence_relation( 0.0 , &guess_1 , &guess_2 , &output_1 , &output_2 , &guess_3 );
	printVector( &guess_3 );
	return;
}

void test_secant_gradient()
{
	vector<double> guess_1 = { 1.0 , 1.0 };
	vector<double> guess_2 = { 2.0 , 3.0 };
	vector<double> output_1 , output_2;
	generic_nonlinear_function( 0.0 , &guess_1 , &output_1 );
	generic_nonlinear_function( 0.0 , &guess_2 , &output_2 );

	vector<double> gradients;
	secantMinimization::gradient_calculation( &guess_2 , &guess_1 , &output_2 , &output_1 , &gradients );
	printVector( &gradients );
	return;
}

void test_secant_hessian()
{
	vector<double> guess_1 = {1.0 , 1.0};
	vector<double> guess_2 = {2.0 , 3.0};
	vector<double> guess_3 = {4.0 , 2.0};
	vector<double> output_1 , output_2 , output_3;
	generic_nonlinear_function( 0.0 , &guess_1 , &output_1 );
	generic_nonlinear_function( 0.0 , &guess_2 , &output_2 );
	generic_nonlinear_function( 0.0 , &guess_3 , &output_3 );
	// cout << output_1[0] << "	" << output_2[0] <<"	"<< output_3[0] << endl;

	vector<vector<double>> hessians;
	secantMinimization::hessian_calculation( &guess_3 , &guess_2 , &guess_1 , &output_3 , &output_2 , &output_1 , &hessians );
	printMatrix( &hessians );
	return;
}