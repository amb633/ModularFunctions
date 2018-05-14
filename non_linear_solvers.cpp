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

void quasiNewtonMinimization::quasiNewton_delta( void (*function)( double , vector<double>* , vector<double>* ) ,
	double time , vector<double> previous_parameters , double pert , 
	vector<double>* parameter_solutions , vector<double>* performance_metrics , vector<double>* desired_result , bool normalize )
{

	int counter = 0;
	double relative_residual = 1.0 ;
	double absolute_residual = 1.0;
	vector<double> updated_parameters;

	while( absolute_residual > 1e-9 ){

		//printVector( &previous_parameters );
		vector<double> previous_output , perturbed_output;
		vector<double> perturbed_parameters;
		vector<double> gradients , deltas;
		vector<vector<double>> perturbed_output_cross , hessians;
		zeroMatrix( previous_parameters.size() , &perturbed_output_cross );
		function( time , &previous_parameters , &previous_output );

		for ( int i = 0 ; i < previous_parameters.size() ; i++ ){
			vector<double> param_dummy = previous_parameters;
			param_dummy[i] = (1.0+pert)*param_dummy[i];
			perturbed_parameters.push_back(param_dummy[i]);
			function( time , &param_dummy , &perturbed_output );
		}

		for ( int i = 0 ; i < previous_parameters.size() ; i++ ){
			for ( int j = 0 ; j < previous_parameters.size() ; j++ ){
				vector<double> param_dummy = previous_parameters;
				if ( i == j ) param_dummy[i] = (1.0 + 2.0*pert)*param_dummy[i];
				else{
					param_dummy[i] = (1.0 + pert)*param_dummy[i];
					param_dummy[j] = (1.0 + pert)*param_dummy[j];
				}
				vector<double> temp;
				function( time , &param_dummy , &temp );
				perturbed_output_cross[i][j] = temp[0];
				perturbed_output_cross[j][i] = temp[0];
			}
		}

		gradient_calculation( &previous_parameters , pert , &previous_output , &perturbed_output , &gradients );
		hessian_calculation( &previous_parameters , pert , &previous_output , &perturbed_output , &perturbed_output_cross, &hessians);
		fullSolver( &hessians , &gradients , &deltas );

		vector<double> deltas_n;
		scaleVector( -1.0 , &deltas , &deltas_n );
        vector<double> scaled_output, scaled_deltas;
        tSearch(function, &deltas_n, &scaled_deltas, &previous_parameters, 0, 1, desired_result);
		addVectors( &previous_parameters , &scaled_deltas , &updated_parameters );

		absolute_residual = 0.0;
		relative_residual = 0.0;

//        for ( int i = 0 ; i < deltas.size() ; i++ ){
//            absolute_residual += scaled_deltas[i]*scaled_deltas[i];
//            relative_residual += (scaled_deltas[i]*scaled_deltas[i]) / (previous_parameters[i] * previous_parameters[i]);
//        }
        vectorNorm(&updated_parameters, &previous_parameters, relative_residual);
        
        vector<double> current_result;
        function( time , &updated_parameters , &current_result );
        vectorNorm(desired_result,&current_result, absolute_residual);
        
		previous_parameters.erase(previous_parameters.begin(), previous_parameters.end());
		previous_parameters = updated_parameters;
		updated_parameters.erase(updated_parameters.begin(), updated_parameters.end());
		counter++;
		if ( counter > 50000 ) break;
	}	
	(*parameter_solutions) = previous_parameters;
	(*performance_metrics).push_back(counter);
	(*performance_metrics).push_back(absolute_residual);
	(*performance_metrics).push_back(relative_residual);
	return;
}

void quasiNewtonMinimization::gradient_calculation( vector<double>* current_parameters , double pert ,
	vector<double>* current_output , vector<double>* perturbed_output , vector<double>* gradients )
{

	for ( int i = 0 ; i < (*current_parameters).size() ; i++ ){
		double num = (*perturbed_output)[i] - (*current_output)[0];
		double den = pert*((*current_parameters)[i]);
		double grad = num / den;
		(*gradients).push_back(grad);
	}
	return;
}

void quasiNewtonMinimization::hessian_calculation( vector<double>* current_parameters , double pert , 
	vector<double>* current_output , vector<double>* perturbed_output , vector<vector<double>>* perturbed_output_cross ,
	vector<vector<double>>* hessians )
{
	int rank = (*current_parameters).size();
	zeroMatrix( rank , hessians );
	for ( int i = 0 ; i < rank ; i++ ){
		for ( int j = 0 ; j < rank ; j++ ){
			double num = (*perturbed_output_cross)[i][j] - (*perturbed_output)[i] - (*perturbed_output)[i] + (*current_output)[0];
			double d1 = pert*((*current_parameters)[i]);
			double d2 = pert*((*current_parameters)[j]);
			double hess = num/ (d1 * d2);
			(*hessians)[i][j] = hess;
		}
	}
	return;
}

void tSearch( void (*function)( double , vector<double>* , vector<double>* ), vector<double>* delta, vector<double>* t_delta, vector<double>* x, double t_min, double t_max, vector<double>* true_result){
    
    vector<double> t_opt_result;
    vector<double> scaled_delta_mid, update_param_mid, result_mid;
    vector<double> scaled_delta_min, update_param_min, result_min;
    vector<double> scaled_delta_max, update_param_max, result_max;
    
    double t_mid = (t_max + t_min)/2;
    
    scaleVector(t_min, delta, &scaled_delta_min);
    scaleVector(t_max, delta, &scaled_delta_max);
    scaleVector(t_mid, delta, &scaled_delta_mid);
    
    addVectors(&scaled_delta_min, x, &update_param_min);
    addVectors(&scaled_delta_max, x, &update_param_max);
    addVectors(&scaled_delta_mid, x, &update_param_mid);
    
    double time = 0; //dummy not needed for this function;
    function( time, &update_param_min, &result_min);
    function( time, &update_param_max, &result_max);
    function( time, &update_param_mid, &result_mid);
    
    double min_result, max_result, mid_result, opt_result;
    vectorNorm(&result_min, true_result, min_result);
    vectorNorm(&result_max, true_result, max_result);
    vectorNorm(&result_mid, true_result, mid_result);
    (*t_delta) = scaled_delta_mid;
    opt_result = mid_result;
    if( mid_result > 1e-7 ){
        if( min_result < opt_result ){
            (*t_delta) = scaled_delta_min;
            tSearch(function, delta, t_delta, x, t_min, t_mid, true_result);
        } else if ( max_result < opt_result ){
            (*t_delta) = scaled_delta_max;
            tSearch(function, delta, t_delta, x, t_mid, t_max, true_result);
        }
    }
    
}
