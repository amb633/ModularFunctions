#include "utility_functions.hpp"

void printMatrix( vector<vector<double>>* matrix )
{
	// assuming square matrix
	int rank = (*matrix).size();
	for ( int i = 0 ; i < rank ; i++ ){
		for ( int j = 0 ; j < rank ; j++ ){
			cout << ((*matrix)[i])[j] << "	";
		}
		cout << endl;
	}

}

void printVector( vector<double>* vector_1 )
{
	int rank = (*vector_1).size();
	for ( int i = 0 ; i < rank ; i++ ){
		cout << (*vector_1)[i] << "   ";
	}
	cout << endl;
}


void identityMatrix( int rank , vector<vector<double>>* I )
{
	// function to create identity matrix of size rank

	for ( int c = 0 ; c < rank ; c++ ){
		vector<double> row;
		for ( int r = 0 ; r < rank ; r++ ){
			if ( r == c ) row.push_back(1.0);
			else row.push_back(0.0);
		}
		(*I).push_back(row);
	}

}

void zeroMatrix( int rank , vector<vector<double>>* z )
{
	// function to create zero matrix of size rank
	vector<double> row;
	for ( int r = 0 ; r < rank ; r++ ){
		row.push_back(0.0);
	}
	for ( int c = 0 ; c < rank ; c++ ){
		(*z).push_back(row);
	}
}

void zeroVector( int rank , vector<double>* v )
{
	// function to create zero vector of size rank
	for ( int i = 0 ; i < rank ; i++ ){
		(*v).push_back(0.0);
	}
}

void copyMatrix( vector<vector<double>>* original , vector<vector<double>>* copy )
{
	// creates a deep copy of the matrix
	int rank = (*original).size();
	zeroMatrix( rank , copy );
	for ( int i = 0 ; i < rank ; i++ ){
		for ( int j = 0 ; j < rank ; j++ ){
			((*copy)[i])[j] = ((*original)[i])[j];
		}
	}
}

void matrixProduct( vector<vector<double>>* matrix_1 , vector<vector<double>>* matrix_2 , vector<vector<double>>* result )
{	
	// calculates result = matrix_1*matrix_2 -> output is a matrix
	int rank = (*matrix_1).size();
	vector<vector<double>> temp;
	zeroMatrix( rank , &temp );
	for ( int i = 0; i < rank ; i++ ){
		for ( int j = 0 ; j < rank ; j++ ){
			for ( int k = 0 ; k < rank ; k++ ){
				temp[i][j] +=((*matrix_1)[i][k])*((*matrix_2)[k][j]);
			}
		}
	}
	(*result) = temp;
}

void vectorProduct( vector<vector<double>>* matrix_1 , vector<double>* vector_1 , vector<double>* result )
{
	// calculates result = matrix_1 * vector_1 -> output is a vector
	int rank = (*matrix_1).size();
	zeroVector( rank , result );
	for ( int i = 0 ; i < rank ; i++ ){
		for ( int j = 0 ; j < rank ; j++ ){
			(*result)[i] += ( ( (*matrix_1)[i][j] ) * (*vector_1)[j] );
		}
	}
}

void addVectors( vector<double>* vector_1 , vector<double>* vector_2 , vector<double>* sum )
{
	// returns element wise addition of sum[i] = vector_1[i] + vector_2[i]
	if ( (*vector_1).size() != (*vector_2).size() ){
        cout << "vectors have different dimensions... cannot add them togther! " << endl;
		return;
	}
	for ( int i = 0 ; i < (*vector_1).size() ; i++ ){
		(*sum).push_back( (*vector_1)[i] + (*vector_2)[i] );
	}
}

void multiplyVectors( vector<double>* vector_1 , vector<double>* vector_2 , vector<double>* product )
{
	// returns element wise multiplication of sum[i] = vector_1[i] * vector_2[i]
	if ( (*vector_1).size() != (*vector_2).size() ){
        cout << "vectors have different dimensions... cannot multiply them togther! " << endl;
		return;
	}
	for ( int i = 0 ; i < (*vector_1).size() ; i++ ){
		(*product).push_back( (*vector_1)[i] * (*vector_2)[i] );
	}
}

void scaleVector( double scalar, vector<double>* a,  vector<double>* result)
{
	// elementwise multiple: result = scalar * a;
    if( (*result) == (*a)){
    	return;
    }
    for( int i = 0; i<(*a).size(); i++){
        (*result).push_back((*a)[i]*scalar);
    }
}

void shiftVector( double scalar, vector<double>* a,  vector<double>* result)
{
	// adds a constant to all th elements in a vector
    for( int i = 0; i<(*a).size(); i++){
        (*result).push_back((*a)[i] + scalar);
    }
}

void expVector( vector<double>* input , vector<double>* result)
{
	// calculates the exponential of each element in a vector
	if( (*input).size() < 1 ){
        cout << "ERROR: cannot calculate the exponential of an empty vector " << endl;
        return;
    }
    for( int i = 0 ; i < (*input).size(); i++ ){
            (*result).push_back(exp((*input)[i]));
    }
}


void vectorNorm( vector<double>* vector_1 , double& result )
{
	// calculates norm of vector_1
	double squareSum = 0.0;
	for ( int i = 0 ; i < (*vector_1).size() ; i++ ){
		double temp = (*vector_1)[i];
		squareSum += temp*temp;
	}
	result = sqrt( squareSum );
}

void vectorNorm( vector<double>* vector_1 , vector<double>* vector_2 , double& result )
{
	// calculates norm of vector_1 - vector_2
	double squareSum = 0.0;
	for ( int i = 0 ; i < (*vector_1).size() ; i++ ){
		double temp = (*vector_1)[i] - (*vector_2)[i];
		squareSum += temp*temp;
	}
	result = sqrt( squareSum );
}