#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>

/* Local Approximation Analysis */

using namespace std;

void generic_polynomial_function ( double time , vector<double>* input , vector<double>* result )
{
	double answer = 7.0*pow(time,3.0) - 8.0*pow(time,2.0) - 3.0*time + 3.0;
	(*result).push_back( answer );
}

namespace local_approximation {

	void forwardEuler( void (*function)(double time , vector<double>* , vector<double>* ) ,
		double time , double march , vector<double>* input , vector<double>* gradient )
	{
		vector<double> result_1 , result_2;

		function( time , input , &result_1 );
		function( time + march , input , &result_2 );

		double feu = ( result_2[0] - result_1[0] ) / march;

		(*gradient).push_back( feu );
	}

	void backwardEuler( void (*function)( double time , vector<double>* , vector<double>* ) ,
		double time , double march , vector<double>* input , vector<double>* gradient )
	{
		vector<double> result_1 , result_2;

		function( time - march , input , &result_1 );
		function( time , input , &result_2 );

		double beu = ( result_2[0] - result_1[0] ) / march;

		(*gradient).push_back( beu );

	}

	void centralEuler ( void (*function)( double time , vector<double>* , vector<double>* ) , 
		double time , double march , vector<double>* input , vector<double>* gradient )
	{
		vector<double> result_1 , result_2;
		function( time - march , input , &result_1 );
		function( time + march , input , &result_2 );

		double ceu = ( result_2[0] - result_1[0] ) / (2.0*march);

		(*gradient).push_back( ceu );
	}
}

int main( void )
{
	cout<<fixed;

	void (*function)( double , vector<double>* , vector<double>* ) = generic_polynomial_function;

	vector<double> input , gradient ;

	local_approximation::forwardEuler( function , 1.0 , 0.1 , &input , &gradient );
	local_approximation::backwardEuler( function , 1.0 , 0.1 , &input , &gradient );
	local_approximation::centralEuler( function , 1.0 , 0.1 , &input , &gradient );

	cout << " forward gradient at 1.0 = " << gradient[0] << endl;
	cout << " backward gradient at 1.0 = " << gradient[1] << endl;
	cout << " central gradient at 1.0 = " << gradient[2] << endl;

}