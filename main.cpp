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
#include "fullSolver.hpp"
#include "non_linear_solvers.hpp"

int main(int argc, const char * argv[]) {
    cout<<fixed;
    
    void (*function)( double , vector<double>* , vector<double>* ) = generic_polynomial_function;
    void (*exp_func)( double , vector<double>* , vector<double>* ) = exponential_function;
    void (*nol_func)( double , vector<double>* , vector<double>* ) = generic_nonlinear_function;
    
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

    cout << endl << endl;

    test_jacobi_iter_solver();
    
    test_direct_solver();

    test_SOR_iter_solver();

    cout << " testing recurrence relation function : " ; 
    test_recurrence_relation();
    cout << " testing secant gradient function : " ; 
    test_secant_gradient();
    cout << " testing secant hessian function : " << endl; 
    test_secant_hessian();

    vector<double> initial_guess_1 , initial_guess_2 , parameter_solutions , perfomance_metrics;
    initial_guess_1.push_back(1.0);
    initial_guess_2.push_back(0.9);
    secantMinimization::secant_delta( exp_func , 0.0 , initial_guess_1 , initial_guess_2 , 
        &parameter_solutions , &perfomance_metrics );

    cout << " testing secant delta function : " ; 
    printVector( &parameter_solutions );
    cout << " perfomance_metrics : " ;
    printVector( &perfomance_metrics );


    vector<double> initial_guess_1b , initial_guess_2b , parameter_solutions_b , perfomance_metrics_b;
    initial_guess_1b.push_back(1.0);
    initial_guess_1b.push_back(1.0);
    initial_guess_2b.push_back(2.0);
    initial_guess_2b.push_back(3.0);

    secantMinimization::secant_delta( nol_func , 0.0 , initial_guess_1b , initial_guess_2b , 
        &parameter_solutions_b , &perfomance_metrics_b );

    cout << " testing secant delta function : " ; 
    printVector( &parameter_solutions_b );
    cout << " perfomance_metrics : " ;
    printVector( &perfomance_metrics_b );


    return 0;
}
