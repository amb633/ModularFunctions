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

void test_quasiNewton_gradient()
{
	vector<double> guess_1 = {1.0 , 1.0};
	double pert = 1e-3;
	vector<double> output_1 , perturbed_parameters , perturbed_output , gradients;
	generic_nonlinear_function( 0.0 , &guess_1 , &output_1 );

	for ( int i = 0 ; i < guess_1.size() ; i++ ){
		vector<double> param_dummy = guess_1;
		param_dummy[i] = (1.0 + pert)*param_dummy[i];
		perturbed_parameters.push_back( param_dummy[i] );
		generic_nonlinear_function( 0.0 , &param_dummy , &perturbed_output );
	}

	quasiNewtonMinimization::gradient_calculation( &guess_1 , pert , &output_1 , &perturbed_output , &gradients );
	printVector( &gradients );
	return;
}

void test_quasiNewton_hessian()
{
	vector<double> guess_1 = {1.0 , 1.0};
	double pert = 1e-3;
	vector<double> output_1 , perturbed_parameters , perturbed_output , gradients;
	generic_nonlinear_function( 0.0 , &guess_1 , &output_1 );

	for ( int i = 0 ; i < guess_1.size() ; i++ ){
		vector<double> param_dummy = guess_1;
		param_dummy[i] = (1.0 + pert)*param_dummy[i];
		perturbed_parameters.push_back( param_dummy[i] );
		generic_nonlinear_function( 0.0 , &param_dummy , &perturbed_output );
	}

	vector<vector<double>> perturbed_output_cross;
	vector<vector<double>> hessians;
	zeroMatrix( guess_1.size() , &perturbed_output_cross );

	for ( int i = 0 ; i < guess_1.size() ; i++ ){
		for ( int j = 0 ; j < guess_1.size() ; j++ ){
			vector<double> param_dummy = guess_1;
			if ( i == j ){
				param_dummy[i] = (1.0 + 2.0*pert)*param_dummy[i];
			}
			else {
				param_dummy[i] = (1.0 + pert)*param_dummy[i];
				param_dummy[j] = (1.0 + pert)*param_dummy[j];
			}
			vector<double> temp;
			generic_nonlinear_function( 0.0 , &param_dummy , &temp );
			perturbed_output_cross[i][j] = temp[0];
//            perturbed_output_cross[j][i] = temp[0];
		}
	}

	quasiNewtonMinimization::hessian_calculation( &guess_1 , pert , &output_1 , &perturbed_output , 
		&perturbed_output_cross , &hessians );
	printMatrix( &hessians );
	return;

}
