//
//  main.cpp
//  FinalExamPrep
//
//  Created by Ariana Bruno on 5/12/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#include <iostream>
#include "local_approximation.hpp"
#include "utility_functions.hpp"

int main(int argc, const char * argv[]) {
    cout<<fixed;
    
    void (*function)( double , vector<double>* , vector<double>* ) = generic_polynomial_function;
    
    vector<double> input , gradient , error;
    
    differentiation_approximation::forwardEuler( function , 1.0 , 0.1 , &input , &gradient );
    differentiation_approximation::backwardEuler( function , 1.0 , 0.1 , &input , &gradient );
    differentiation_approximation::centralEuler( function , 1.0 , 0.1 , &input , &gradient );
    differentiation_approximation::secondTaylor( function , 1.0 , 0.1 , &input , &gradient );
    differentiation_approximation::error_calculation::richardsonEstimation( function , 1.0 , 0.1 , &input , &error );
    
    vector<double> integration;
    integration_approximation::rectangle( function ,-1.0 , 2.0 , &input , &integration );
    integration_approximation::trapezoid( function ,-1.0 , 2.0 , &input , &integration );
    integration_approximation::midpoint( function ,-1.0 , 2.0 , &input , &integration );
    integration_approximation::simpson( function ,-1.0 , 2.0 , &input , &integration );
    integration_approximation::gaus2pt( function ,-1.0 , 2.0 , &input , &integration );
    
    cout << " forward gradient at 1.0 = " << gradient[0] << endl;
    cout << " backward gradient at 1.0 = " << gradient[1] << endl;
    cout << " central gradient at 1.0 = " << gradient[2] << endl;
    cout << " 2nd order taylor grad at 1.0 = " << gradient[3] << endl;

    cout << " richardson error estimation = " << error[0] << endl;
    
    cout << " rectangle integration from -1.0 and 1.0 = " << retrieve_element(&integration) << endl;
    cout << " trapezoid integration at -1.0 and 1.0 = " << retrieve_element(&integration) << endl;
    cout << " midpoint integration at -1.0 and 1.0 = " << retrieve_element(&integration) << endl;
    cout << " simpson integration at -1.0 and 1.0 = " << retrieve_element(&integration) << endl;
    cout << " gaussian 2-pt integration at -1.0 and 1.0 = " << retrieve_element(&integration) << endl;
    
    return 0;
}
