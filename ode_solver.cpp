#include "ode_solver.hpp"

void ode_exponential_function( double time , vector<double>* input , vector<double>* output )
{
	vector<double> temp;
	double exp_comp = 4.0*exp( 0.8*time );
	scaleVector( -0.5 , input , &temp );
	shiftVector( exp_comp , &temp , output );
}

void ground_truth( double time , vector<double>* input , vector<double>* output )
{
	double temp = 4.0*(exp( 0.8*time ) - exp( -0.5*time ))/1.3;
	double ans = temp + 2.0*exp( -0.5*time );
	(*output).push_back( ans );
}

void ode_solvers::forward_euler ( void (*function)( double, vector<double>* , vector<double>* ) ,
	double time , double march , vector<double>* input , vector<double>* slope )
{
	function( time , input , slope );
	return;
}

void ode_solvers::heun_oneStep( void (*function)( double , vector<double>* , vector<double>* ) ,
	double time , double march , vector<double>* input , vector<double>* slope )
{
	vector<double> f1 , f2 , eu , eu_update , predict , sum;
	
	// get the prediction using the forward euler method
	forward_euler( function , time , march , input , &eu );
	scaleVector( march , &eu , &eu_update );
	addVectors( input , &eu_update , &predict );
	
	function( time , input , &f1 );
	function( time + march , &predict , &f2 );

	addVectors( &f1 , &f2 , &sum );
	scaleVector( 0.5 , &sum , slope );
	return;
}

void ode_solvers::ODE_SOLVER( void (*function)(double, vector<double>* , vector<double>* ),
	double time , double march , vector<double>* input , vector<double>* new_values , int method )
{
	if ( time == 0.0 ){
		(*new_values) = (*input);
		return;
	}
	switch( method ){
		case FORWARD_EULER: {
			vector<double> phi , update;
			forward_euler( function , time , march , input , &phi );
			scaleVector( march , &phi , &update );
			addVectors( input , &update , new_values );
			break;
		}
		case HEUN_ONE: {
			vector<double> phi , update;
			heun_oneStep( function , time , march , input , &phi );
			scaleVector(march , &phi , &update);
			addVectors( input , &update , new_values );
			break;
		}
	}
	return;
}

void test_forward_euler()
{
	void (*ode_exp_fcn)( double , vector<double>* , vector<double>* ) = ode_exponential_function;
	vector<double> old_values = {2.0};
	double march = 1.0;
	for ( double time = 0.0 ; time < 5.0 ; time += march ){
		vector<double> phi , update , new_values , actual_values;
		ode_solvers::forward_euler( ode_exp_fcn , time , march , &old_values , &phi );
		scaleVector(march , &phi , &update);
		addVectors( &old_values , &update , &new_values );
		ground_truth( time , &old_values , &actual_values );
		printVector( &new_values );
		// cout << "	" ;
		// printVector( &actual_values );
		old_values.erase(old_values.begin(), old_values.end());
		old_values = new_values;
	}
	return;
}

void test_heun_oneStep()
{
	void (*ode_exp_fcn)( double , vector<double>* , vector<double>* ) = ode_exponential_function;
	vector<double> old_values = {2.0};
	double march = 1.0;
	for ( double time = 0.0 ; time < 5.0 ; time += march ){
		vector<double> phi , update , new_values , actual_values;
		ode_solvers::heun_oneStep( ode_exp_fcn , time , march , &old_values , &phi );
		scaleVector(march , &phi , &update);
		addVectors( &old_values , &update , &new_values );
		ground_truth( time , &old_values , &actual_values );
		printVector( &new_values );
		// cout << "	" ;
		// printVector( &actual_values );
		old_values.erase(old_values.begin(), old_values.end());
		old_values = new_values;
	}
}