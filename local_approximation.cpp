#include "local_approximation.hpp"

/* Local Approximation Analysis */

void generic_polynomial_function ( double time , vector<double>* input , vector<double>* result )
{
    double answer = 7.0*pow(time,3.0) - 8.0*pow(time,2.0) - 3.0*time + 3.0;
    (*result).push_back( answer );
}

namespace differentiation_approximation 
{

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

	void error_calculation::richardsonEstimation( void (*function)(double , vector<double>* , vector<double>* ) ,
		double time , double march , vector<double>* input , vector<double>* error )
	{
		vector<double> gradient;
		differentiation_approximation::forwardEuler( function , time , march*4.0 , input , &gradient );
		differentiation_approximation::forwardEuler( function , time , march*2.0 , input , &gradient );
		differentiation_approximation::forwardEuler( function , time , march 	 , input , &gradient );

		double err = ( gradient[0] - gradient[1] ) / ( gradient[1] - gradient[2] );
		(*error).push_back( err );
	}
}

namespace integration_approximation {
    
    void rectangle( void(*function)( double time, vector<double>*, vector<double>*),
                   double time, double march, vector<double>* input, vector<double>* integration )
    {
        vector<double> result;
        function(time, input, &result);
        
        double rect_integation = march*result[0];
        
        (*integration).push_back( rect_integation );
    }
    
    void trapezoid(void(*function)( double time, vector<double>*, vector<double>*),
                   double time, double march, vector<double>* input, vector<double>* integration )
    {
        vector<double> result_1, result_2;
        function(time, input, &result_1);
        function(time+march, input, &result_2);
        
        double trap_integation = march*0.5*(result_1[0] + result_2[0]);
        
        (*integration).push_back( trap_integation );
        
    }
    
    void midpoint(void(*function)( double time, vector<double>*, vector<double>*),
                  double time, double march, vector<double>* input, vector<double>* integration )
    {
        vector<double> result;
        double midpt = time + (march/2.0);
        function(midpt, input, &result);
        
        double midpt_integation = march*result[0];
        
        (*integration).push_back( midpt_integation );
    }
    
    void simpson(void(*function)( double time, vector<double>*, vector<double>*),
                 double time, double march, vector<double>* input, vector<double>* integration )
    {
        double pt1 = time;
        double pt2 = time + (march/2.0);
        double pt3 = time + march;
        
        vector<double> result_1, result_2, result_3;
        
        function(pt1, input, &result_1);
        function(pt2, input, &result_2);
        function(pt3, input, &result_3);
        
        double simp_integation = march*(1.0/6.0)*(result_1[0] + 4.0*result_2[0] + result_3[0]);
        
        (*integration).push_back( simp_integation );
    }
    
    
    void gaus2pt(void(*function)( double time, vector<double>*, vector<double>*),
                 double time, double march, vector<double>* input, vector<double>* integration )
    {
        double midpt = time + (march/2.0);
        double pt1 = midpt - ((1.0/(2*sqrt(3)))*march);
        double pt2 = midpt + ((1.0/(2*sqrt(3)))*march);
        
        vector<double> result_1, result_2;
        
        function(pt1, input, &result_1);
        function(pt2, input, &result_2);
        
        double simp_integation = march*(0.5)*(result_1[0] + result_2[0]);
        
        (*integration).push_back( simp_integation );
    }
}

double retrieve_element( vector<double>* result ){
    double first_element = (*result)[0];
    (*result).erase(result->begin());
    return first_element;
}



