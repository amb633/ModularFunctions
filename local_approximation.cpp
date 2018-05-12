#include "local_approximation.hpp"

/* Local Approximation Analysis */


void generic_polynomial_function ( double time , vector<double>* input , vector<double>* result )
{
	double answer = 7.0*pow(time,3.0) - 8.0*pow(time,2.0) - 3.0*time + 3.0;
	(*result).push_back( answer );
}

namespace local_approximation {

	void forwardEuler( void (*function)(double , vector<double>* input, vector<double>* ) ,
		double time , double march , vector<double>* input , vector<double>* gradient )
	{
		vector<double> result_1 , result_2;

		function( time , input , &result_1 );
		function( time + march , input , &result_2 );

		double feu = ( result_2[0] - result_1[0] ) / march;

		(*gradient).push_back( feu );
	}

	void backwardEuler( void (*function)( double , vector<double>* , vector<double>* ) ,
		double time , double march , vector<double>* input , vector<double>* gradient )
	{
		vector<double> result_1 , result_2;

		function( time - march , input , &result_1 );
		function( time , input , &result_2 );

		double beu = ( result_2[0] - result_1[0] ) / march;

		(*gradient).push_back( beu );

	}

	void centralEuler ( void (*function)( double , vector<double>* , vector<double>* ) , 
		double time , double march , vector<double>* input , vector<double>* gradient )
	{
		vector<double> result_1 , result_2;
		function( time - march , input , &result_1 );
		function( time + march , input , &result_2 );

		double ceu = ( result_2[0] - result_1[0] ) / (2.0*march);

		(*gradient).push_back( ceu );
	}

	void secondTaylor( void (*function)(double , vector<double>* , vector<double>* ) ,
		double time , double march , vector<double>* input , vector<double>* gradient )
	{
		vector<double> result_1 , result_2 , result_3;
		function( time , input , &result_1 );
		function( time + march , input , &result_2 );
		function( time + 2.0 * march , input , &result_3 );

		double sta = -(result_3[0]/(2.0*march)) - 3.0*(result_1[0]/(2.0*march)) + 2.0*result_2[0]/march;
		(*gradient).push_back(sta);	
	}
}


