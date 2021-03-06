#include "fullSolver.hpp"
#include "utility_functions.hpp"

void fullSolver( vector<vector<double>>* matrix , vector<double>* b , vector<double>* solution ){
	
	// get some information about the matrix
	int rank = (*matrix).size();
	
	// condition the matrix through partial row pivoting
	conditionMatrix( matrix , b );
	//printMatrix( matrix );

	// create a dummy copy of matrix
	vector< vector<double>> dummy;
	copyMatrix( matrix , &dummy );

	// number of M matrices required = rank - 1;
	vector<vector<vector<double>>> M_matrices;
	for ( int i = 0 ; i < (rank - 1) ; i++ ){
		
		// first create an identity matrix
		vector<vector<double>> m;
		identityMatrix( rank , &m );

		// second figure out the off-diagonal elements 
		calcAtomicVector( &dummy , i , &m );
		M_matrices.push_back(m);

		// update dummy matrix for next iteration
		matrixProduct( &m , &dummy , &dummy );

	}

	// multiply all the m matrices in reverse order
	vector<vector<double>> M;
	identityMatrix( rank , &M );
	
	for ( int i = M_matrices.size() - 1 ; i >= 0 ; i-- ){
		matrixProduct( &M , &M_matrices[i] , &M );
	}
	
	// no need to explicity calculate lower triangular matrix
		// Ux = (L^-1)*b where L = (M1^-1)(M2^-1)(M3^-1)...
		// so (L^-1) = ...(M3)(M2)(M1) -> in reverse order (we're already calculated this for the U matrix)
		// so Ux = y = M*b
	vector<double> y;
	vectorProduct( &M , b , &y );
	
	// define upper triangular matrix
	vector<vector<double>> U;
	matrixProduct( &M , matrix , &U );

	// backward substitution -> x = U\y
	backwardSubstitution( &U , &y , solution );

	return;
}

void conditionMatrix( vector<vector<double>>* matrix , vector<double>* b ){
// function to condition matrix through partial row pivoting
	// assume square matrix
	int rank = (*matrix).size();

	// iterate over each column
	for ( int col = 0 ; col < rank ; col++ ){
		//printMatrix(matrix);

		double max = 0.0;
		int max_row = col;

		// search for the max value along current column
		for ( int i = col ; i < rank ; i++ ){
			double v = abs((*matrix)[i][col]);
			if ( v > max ) {
				max = v;
				max_row = i;
			}
		}
		// swap row with max_row
		vector<double> temp = ((*matrix)[col]);
		((*matrix)[col]) = ((*matrix)[max_row]);
		((*matrix)[max_row]) = temp;
		double t = (*b)[col];
		(*b)[col] = (*b)[max_row];
		(*b)[max_row] = t;
	}

}

void calcAtomicVector( vector<vector<double>>* matrix, int row , vector<vector<double>>* m ){
	// calculates the off diagonal elements for gaussian elimination where the current pivot is matrix[row][row]
	// the coefficients are stored in m
	int rank = (*matrix).size();
	double pivot = (*matrix)[row][row];
	for ( int i = (row + 1) ; i < rank ; i++ ){
		(*m)[i][row] = -(*matrix)[i][row]/pivot;
	}

}

void backwardSubstitution( vector<vector<double>>* U , vector<double>* y , vector<double>* x ){
	// backward sub to find x in Ux = y where U is upper triangular and y = (L^-1)b
	int rank = (*y).size();
	for ( int i = 0 ; i < rank ; i++ ){
		(*x).push_back(0.0);
	}
	int end = rank - 1;
	double temp = (*y)[end] /((*U)[end])[end];
	(*x)[end] = temp;
	double sum;

		for ( int i = end-1 ; i >= 0 ; i-- ){
		sum = 0.0;
		for ( int j = end ; j > i ; j-- ){
			sum += (*x)[j]*((*U)[i])[j];
		}
		temp = ( (*y)[i] - sum )/((*U)[i])[i];
		(*x)[i] = temp;
	}

}

void get_diagInv( vector<vector<double>>* input, vector<double>* diag ){
    for( int i=0; i<input->size(); i++){
        for( int j=0; j< (*input)[0].size(); j++){
            if( i == j){
                (*diag).push_back(1.0/(*input)[i][j]);
            }
        }
    }
}

/* decompose AF matrix into diagonal elements (stored in DF) and non-diagonal elements (stored in LUF) */
void decomposeMatrix( vector<double>* DF , vector< vector<double>>* LUF , vector< vector<double>>* AF ){
    
    copyMatrix( AF, LUF );
    
    for ( int i = 0 ; i < (*AF).size() ; i++ ){
        // extract the diagonal elements to DS
        (*DF).push_back( retrieveElement( AF , i , i ) );
        // change the diagonal elements in LUF to 0
        changeElement( LUF , i , i , 0.0 );
    }
    
    // negate all the elements in LUF
    scaleMatrix( LUF , -1.0 );
}

void jacobiIter( vector<double>* X , vector<double>* DF , vector< vector<double>>* LUF , vector<double>* BF ){
    
    vector<double> matPdt;
    vectorProduct( LUF, X, &matPdt);
    
    for ( int i = 0 ; i < (*DF).size() ; i++ ){
        double dInv = 1.0/((*DF)[i]);
        (*X)[i] = ((*BF)[i] + matPdt[i])*dInv;
    }
}

void jacobiSolver( vector<vector<double>>* A, vector<double>* b, vector<double>* x, double tol ){
    vector<double> D;
    vector<vector<double>> L;
    decomposeMatrix( &D , &L , A );
    
    double normPrev = 2;
    double normCurrent = 1;
    while ( abs(normCurrent - normPrev) > tol ) {
        normPrev = normCurrent;
        jacobiIter(x, &D, &L, b);
        vector<double> Ax;
        vectorProduct(A, x, &Ax);
        vectorNorm(b, &Ax, normCurrent);
    }
}

void SORIter( vector<vector<double>>* A, vector<double>* b, vector<double>* x, double omega){
    
    vector<double> r, product_Dr, product_Ax, neg_product_Ax, iter_adj;
    vector<double> vector_dummy = *x;
    
    //calculating A*x_k-1
    vectorProduct(A, x, &product_Ax);
    
    //solving for r_k-1
    scaleVector(-1.0, &product_Ax, &neg_product_Ax);
    addVectors(b, &neg_product_Ax, &r);
    
    vector<double> D_inv;
    get_diagInv(A, &D_inv);
    
    //solving for iterative adjustment for x_k-1 w * D^-1 * r_k-1
    vectorProduct(&D_inv, &r, &product_Dr);
    scaleVector( omega, &product_Dr, &iter_adj);
    
    //solving for the current x_k
    vectorSum(x, &iter_adj, x);
    vector_dummy = *x;
    
}

void SORSolver( vector<vector<double>>* A, vector<double>* b , vector<double>* x, double tol, double omega){

    vector<double> x_SOR = *x;
    double norm_SOR = 1.0;
    size_t size = x_SOR.size();
    
    
    while(norm_SOR > tol ){
        
        SORIter(A, b, &x_SOR, omega);
        
        //calculating the residual norm for (b - Ax_k)
        vector<double> new_product_Ax;
        vectorProduct(A, &x_SOR, &new_product_Ax);
        vectorNorm(&new_product_Ax, b, norm_SOR);
        
    }

}
