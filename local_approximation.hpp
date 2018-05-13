//
//  local_approx.hpp
//  
//
//  Created by Ariana Bruno on 5/12/18.
//

#ifndef local_approx_hpp
#define local_approx_hpp
#pragma once
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

void generic_polynomial_function ( double time , vector<double>* input , vector<double>* result );

namespace differentiation_approximation
{
    
    void forwardEuler( void (*function)(double , vector<double>* , vector<double>* ) ,
                      double time , double march , vector<double>* input , vector<double>* gradient );

    void backwardEuler( void (*function)( double , vector<double>* , vector<double>* ) ,
                       double time , double march , vector<double>* input , vector<double>* gradient );
    
    void centralEuler ( void (*function)( double , vector<double>* , vector<double>* ) ,
                       double time , double march , vector<double>* input , vector<double>* gradient );

    void secondTaylor( void (*function)( double , vector<double>* , vector<double>* ) , 
    				   double time , double march , vector<double>* input , vector<double>* gradient );

    namespace error_calculation
    {
    	void richardsonEstimation( void (*function)( double , vector<double>* , vector<double>* ) ,
    						  double time , double march , vector<double>* input , vector<double>* error );
    }

};

namespace integration_approximation {
    
    void rectangle( void(*function)( double time, vector<double>*, vector<double>*),
                   double time, double march, vector<double>* input, vector<double>* integration );
    
    void trapezoid(void(*function)( double time, vector<double>*, vector<double>*),
                   double time, double march, vector<double>* input, vector<double>* integration );
    
    void midpoint(void(*function)( double time, vector<double>*, vector<double>*),
                   double time, double march, vector<double>* input, vector<double>* integration );
    
    void simpson(void(*function)( double time, vector<double>*, vector<double>*),
                  double time, double march, vector<double>* input, vector<double>* integration );
    
    void gaus2pt(void(*function)( double time, vector<double>*, vector<double>*),
                 double time, double march, vector<double>* input, vector<double>* integration );
};

double retrieve_element( vector<double>* result );


#endif /* local_approx_hpp */
