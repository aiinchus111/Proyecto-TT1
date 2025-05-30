#ifndef _MATRIX_
#define _MATRIX_

#include <iostream>
#include <vector>

using namespace std;

class Matrix {
public:
    int n_row, n_column;
    double **data;

    // Constructor de la clase Matrix
    Matrix(int n_row, int n_column);
    Matrix ();
     // Destructor
    // ~Matrix();

    // Constructor de la clase Matrix (n x n)
    Matrix(int n);

    double &operator()(int n);

    // Member operators
    double& operator () (int row, int column);
    

    Matrix& operator + (Matrix &m);
    Matrix& operator - (Matrix &m);
    Matrix& operator * (Matrix &m);
    Matrix& operator / (Matrix &m);
    Matrix& operator = (const Matrix& m);
    
     // Sobrecarga de operadores para operaciones con escalares
    Matrix& operator + (double scalar);
    Matrix& operator - (double scalar);
    Matrix& operator * (double scalar);
    Matrix& operator / (double scalar);

    // Métodos
   
    Matrix& transpose();
    Matrix& inv();
    double determinant();
    double cofactor(int row, int column) ;
    
    double norm();
    double dot(Matrix& m) ;
    Matrix&  cross( Matrix& m) ;
    Matrix&  extract_vector(int n, int m);
 
	Matrix& extract_row(int row) ;
    Matrix& extract_column(int column) ;
    Matrix& assign_row(int row, Matrix v);
    Matrix& assign_column(int column, Matrix v);

    // Non-member operators
    friend ostream& operator << (ostream &o, Matrix &m);
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

// Métodos
Matrix& zeros(const int n_row, const int n_column);
 Matrix& eye(int n);
Matrix& union_vector( Matrix v, Matrix v2);
Matrix& zeros(int n);
#endif
