#include "non_linear_solvers.hpp"


void secantMinimization::secant_delta( void (*function)( double , vector<double>* , vector<double>* ) ,
	double time , vector<double> initial_guess_1 , vector<double> initial_guess_2 , 
	vector<double>* parameter_solutions , vector<double>* performance_metrics , bool normalize )
{
	// evaluate function for first guess
	vector<double> output_1;
	function( time , &initial_guess_1 , &output_1 );

	// evaluate function for second guess
	vector<double> output_2;
	function( time , &initial_guess_2 , &output_2 );

	// use recurrence relation to find third guess
	vector<double> initial_guess_3;
	recurrence_relation( time , &initial_guess_1 , &initial_guess_2 , &output_1 , &output_2 , &initial_guess_3 );
	// evaluate function for third guess
	vector<double> output_3;
	function( time , &initial_guess_3 , &output_3 );

	// initialize vectors for first iteration
	vector<double> gradients;
	vector<vector<double>> hessians;
	vector<double> deltas;

	// calculate the gradient 
	gradient_calculation( &initial_guess_3 , &initial_guess_2 , &output_3 , &output_2 , &gradients );

	// calculate the hessian
	hessian_calculation( &initial_guess_3 , &initial_guess_2 , &initial_guess_1 , 
		&output_3 , &output_2 , &output_1 , &hessians );

	// call the solver function to find the deltas
	fullSolver( &hessians , &gradients , &deltas );

	// do a line search

	int counter = 0;
	double relative_residual = 1.0 ;
	double absolute_residual;

	while( relative_residual > 1e-9 ){
		// erase the oldest position , and move the positions by 1
		initial_guess_1.erase( initial_guess_1.begin(), initial_guess_1.end());
		initial_guess_1 = initial_guess_2;
		output_1[0] = output_2[0];

		initial_guess_2.erase(initial_guess_2.begin(), initial_guess_2.end());
		initial_guess_2 = initial_guess_3;
		output_2[0] = output_3[0];

		initial_guess_3.erase(initial_guess_3.begin(), initial_guess_3.end());
		output_3.erase(output_3.begin(),output_3.end());

		// update current parameters
		vector<double> deltas_n;
		scaleVector( -1.0 , &deltas , &deltas_n );
		addVectors( &initial_guess_2 , &deltas_n , &initial_guess_3 );
		function( time , &initial_guess_3 , &output_3 );

		// erase data for next iteration
		gradients.erase(gradients.begin(), gradients.end());
		hessians.erase(hessians.begin(), hessians.end());
		deltas.erase(deltas.begin(), deltas.end());

		// calculate next update iteration
		gradient_calculation( &initial_guess_3 , &initial_guess_2 , &output_3 , &output_2 , &gradients );
		hessian_calculation( &initial_guess_3 , &initial_guess_2 , &initial_guess_1 , 
			&output_3 , &output_2 , &output_1 , &hessians );

		fullSolver( &hessians , &gradients , &deltas );

		absolute_residual = 0.0;
		relative_residual = 0.0;
		// cout << "size of delta : " << deltas.size() << endl;
		for ( int i = 0 ; i < deltas.size() ; i++ ){
			absolute_residual += deltas[i]*deltas[i];
			relative_residual += (deltas[i]*deltas[i]) / (initial_guess_3[i] * initial_guess_3[i]);
		}

		counter++;
		if ( counter > 50000 ) break;
	}

	(*parameter_solutions) = initial_guess_3;
	(*performance_metrics).push_back(counter);
	(*performance_metrics).push_back(absolute_residual);
	(*performance_metrics).push_back(relative_residual);
	return;
}

void secantMinimization::recurrence_relation( double time , vector<double>* initial_guess_1 , vector<double>* initial_guess_2 , 
	vector<double>* output_1 , vector<double>* output_2 , vector<double>* result )
{

	// assuming the function has a single output
	// otherwise calculate the norm of the vectors (?)
	double op1 = (*output_1)[0];
	double op2 = (*output_2)[0];

	// cout << "output_1 = " << op1 << endl;
	// cout << "output_2 = " << op2 << endl;

	for ( int i = 0 ; i < (*initial_guess_1).size() ; i++ ){
		double temp = (((*initial_guess_1)[i])*op2 - ((*initial_guess_2)[i])*op1) / (op2 - op1);
		(*result).push_back(temp);
	}

	return;
}

void secantMinimization::gradient_calculation( vector<double>* current_parameters , vector<double>* previous_parameters , 
	vector<double>* current_output , vector<double>* previous_output , vector<double>* gradients )
{
	if ( (*current_parameters).size() != (*previous_parameters).size() ){
		cout << "error in parameter history in gradient calculation " << endl;
		return;
	}
	// assuming the function has a single output
	// otherwise calculate the norm of the vectors (?)
	double diff = (*current_output)[0] - (*previous_output)[0];

	for ( int i = 0 ; i < (*current_parameters).size() ; i++ ){
		double current = (*current_parameters)[i];
		double prev =  (*previous_parameters)[i];
		double grad = ( diff ) / ( current - prev );
		(*gradients).push_back(grad);
	}
	return;
}

void secantMinimization::hessian_calculation( vector<double>* current_parameters , vector<double>* previous_parameters_1 , vector<double>* previous_parameters_2 , 
	vector<double>* current_output , vector<double>* previous_output_1 , vector<double>* previous_output_2 ,
	vector<vector<double>>* hessians )
{
	if ( (*current_parameters).size() != (*previous_parameters_1).size() || (*current_parameters).size() != (*previous_parameters_2).size() ){
		cout << "error in parameter history in hessian calculation " << endl;
		return;
	}
	int rank = (*current_parameters).size();
	zeroMatrix( rank , hessians );

	double diff = ((*current_output)[0]) - 2.0*((*previous_output_1)[0]) + ((*previous_output_2)[0]);

	for ( int i = 0 ; i < rank ; i++ ){
		double current = (*current_parameters)[i];
		double prev_1 = (*previous_parameters_1)[i];
		double prev_2 = (*previous_parameters_2)[i];
		double hess = diff/ ((current - prev_1)*(prev_1 - prev_2));
		(*hessians)[i][i] = hess;
	}
	return;
}