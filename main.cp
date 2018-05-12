//
//  main.cpp
//  FinalExamPrep
//
//  Created by Ariana Bruno on 5/12/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#include <iostream>
#include "local_approximation.hpp"

int main(int argc, const char * argv[]) {
    cout<<fixed;
    
    void (*function)( double , vector<double>* , vector<double>* ) = generic_polynomial_function;
    
    vector<double> input , gradient ;
    
    local_approximation::forwardEuler( function , 1.0 , 0.1 , &input , &gradient );
    local_approximation::backwardEuler( function , 1.0 , 0.1 , &input , &gradient );
    local_approximation::centralEuler( function , 1.0 , 0.1 , &input , &gradient );
    
    cout << " forward gradient at 1.0 = " << gradient[0] << endl;
    cout << " backward gradient at 1.0 = " << gradient[1] << endl;
    cout << " central gradient at 1.0 = " << gradient[2] << endl;
    
    return 0;
}
