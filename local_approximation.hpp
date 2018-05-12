//
//  local_approx.hpp
//  
//
//  Created by Ariana Bruno on 5/12/18.
//

#ifndef local_approx_hpp
#define local_approx_hpp

#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

void generic_polynomial_function ( double time , vector<double>* input , vector<double>* result );

namespace local_approximation
{
    
    void forwardEuler( void (*function)(double time , vector<double>* , vector<double>* ) ,
                      double time , double march , vector<double>* input , vector<double>* gradient );

    void backwardEuler( void (*function)( double time , vector<double>* , vector<double>* ) ,
                       double time , double march , vector<double>* input , vector<double>* gradient );
    
    void centralEuler ( void (*function)( double time , vector<double>* , vector<double>* ) ,
                       double time , double march , vector<double>* input , vector<double>* gradient );
}


#endif /* local_approx_hpp */
