#include "../include/matrix.h"
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

Matrix::Matrix () {
	this->n_column=0;
	this->n_row=0;
	this->data=nullptr;
}

// Constructor de la clase Matrix
Matrix::Matrix(int n_row, int n_column) : n_row(n_row), n_column(n_column) {
    if (n_row <= 0 || n_column <= 0) {
        cout << "Matrix create: error in n_row/n_column" << endl;
        exit(EXIT_FAILURE);
    }

    // Asignar memoria para la matriz
    data = new double*[n_row];
    if (data == nullptr) {
        cout << "Matrix create: error in data" << endl;
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n_row; ++i) {
        data[i] = new double[n_column];
        if (data[i] == nullptr) {
             cout << "Matrix create: error in data["<<i<<"]" << endl;
             exit(EXIT_FAILURE);
        }
    }
}


// // Destructor de la clase Matrix
// Matrix::~Matrix() {
//     if (data) {
//         for (int i = 0; i < n_row; ++i) {
//             delete[] data[i];
//         }
//         delete[] data;
//     }
// }

/**
 * @brief acceder a un elemento 
 * 
 * @param row la fila a la que accedo
 * @param column la columna a la que accedo
 * @return double& el elemento pedido
 */
double& Matrix::operator()(int row, int column) {
    if (row <= 0 || row > n_row || column <= 0 || column > n_column) {
        cout << "Matrix get: error in row/column" << endl;
        exit(EXIT_FAILURE);
    }
    return data[row - 1][column - 1];
}


// Sobrecarga del operador + para la suma de matrices
Matrix& Matrix::operator+(Matrix &m) {
    if (n_row != m.n_row || n_column != m.n_column) {
        cout << "Matrix sum: error in n_row/n_column" << endl;
        exit(EXIT_FAILURE);
    }

    Matrix *m_aux = new Matrix(n_row, n_column);
    for (int i = 1; i <= n_row; ++i) {
        for (int j = 1; j <= n_column; ++j) {
            (*m_aux)(i, j) = (*this)(i, j) + m(i, j);
        }
    }
    return *m_aux;
}

// Sobrecarga del operador - para la resta de matrices
Matrix& Matrix::operator-(Matrix &m) {
    if (n_row != m.n_row || n_column != m.n_column) {
        cout << "Matrix sub: error in n_row/n_column" << endl;
        exit(EXIT_FAILURE);
    }

    Matrix *m_aux = new Matrix(n_row, n_column);
    for (int i = 1; i <= n_row; ++i) {
        for (int j = 1; j <= n_column; ++j) {
            (*m_aux)(i, j) = (*this)(i, j) - m(i, j);
        }
    }
    return *m_aux;
}

// Sobrecarga del operador * para la multiplicación de matrices
Matrix& Matrix::operator*(Matrix &m) {
    if (n_column != m.n_row) {
        cout << "Matrix multiply: incompatible dimensions" << endl;
        exit(EXIT_FAILURE);
    }

    Matrix *m_aux = new Matrix(n_row, m.n_column);
    for (int i = 1; i <= n_row; ++i) {
        for (int j = 1; j <= m.n_column; ++j) {
            (*m_aux)(i, j) = 0;
            for (int k = 1; k <= n_column; ++k) {
                (*m_aux)(i, j) += (*this)(i, k) * m(k, j);
            }
        }
    }
    return *m_aux;
}

// Sobrecarga del operador / para la división de matrices (utiliza la inversa)
Matrix& Matrix::operator/(Matrix &m) {
    // La división de matrices A / B se define como A * inv(B)
    Matrix temp = m.inv();
    return *this * temp;
}

// Sobrecarga del operador = para la asignación de matrices
Matrix& Matrix::operator=(const Matrix& m) {
    if (this != &m) { // Evitar la autoasignación
        // Si las dimensiones son diferentes, hay que redimensionar la matriz
        if (n_row != m.n_row || n_column != m.n_column) {
            // Liberar la memoria actual
            for (int i = 0; i < n_row; ++i) {
                delete[] data[i];
            }
            delete[] data;

            // Asignar nueva memoria
            n_row = m.n_row;
            n_column = m.n_column;
            data = new double*[n_row];
            if (data == nullptr) {
                cout << "Matrix assign: error in data" << endl;
                exit(EXIT_FAILURE);
            }
            for (int i = 0; i < n_row; ++i) {
                data[i] = new double[n_column];
                 if (data[i] == nullptr) {
                    cout << "Matrix assign: error in data["<<i<<"]" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }

        // Copiar los datos
        for (int i = 0; i < n_row; ++i) {
            for (int j = 0; j < n_column; ++j) {
                data[i][j] = m.data[i][j];
            }
        }
    }
    return *this;
}

ostream& operator<<(ostream& o, Matrix& m) {
    for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
            o << fixed << setprecision(20) << m(i, j) << " ";
        o << "\n";
    }

    return o;
}

// Función para crear una matriz de ceros
Matrix& zeros(const int n_row, const int n_column) {
    Matrix *m_aux = new Matrix(n_row, n_column);
    for (int i = 1; i <= n_row; i++) {
        for (int j = 1; j <= n_column; j++) {
            (*m_aux)(i, j) = 0;
        }
    }
    return (*m_aux);
}

// Función para crear una matriz identidad
Matrix& eye(const int n) {
    Matrix *m_aux = new Matrix(n, n);
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (i == j)
                 (*m_aux)(i, j) = 1.0;
            else
                (*m_aux)(i, j) = 0.0;
        }
    }
    return *m_aux;
}

// Función para calcular la transpuesta de una matriz
Matrix& Matrix::transpose()  {
    Matrix *m_aux = new Matrix(n_column, n_row);
    for (int i = 1; i <= n_row; ++i) {
        for (int j = 1; j <= n_column; ++j) {
            (*m_aux)(j, i) = (*this)(i, j);
        }
    }
    return *m_aux;
}

double Matrix::determinant() {
    if (n_row != n_column) {
        printf("Matrix must be square to calculate determinant.\n");
		exit(EXIT_FAILURE);
    }
    if (n_row == 1) {
        return data[0][0];
    }
    if (n_row == 2) {
        return data[0][0] * data[1][1] - data[0][1] * data[1][0];
    }
    double det = 0.0;
    Matrix submatrix(n_row - 1, n_column - 1);
    for (int j = 0; j < n_column; ++j) {
        for (int i = 1; i < n_row; ++i) {
            for (int k = 0, l = 0; k < n_column; ++k) {
                if (k != j) {
                    submatrix.data[i - 1][l++] = data[i][k];
                }
            }
        }
        det += data[0][j] * cofactor(0, j);
    }
    return det;
}

double Matrix::cofactor(int row, int column) {
    if (n_row != n_column) {
        printf("Matrix must be square to calculate cofactor.\n");
		exit(EXIT_FAILURE);
    }
    Matrix submatrix(n_row - 1, n_column - 1);
    for (int i = 0, sub_i = 0; i < n_row; ++i) {
        if (i != row) {
            for (int j = 0, sub_j = 0; j < n_column; ++j) {
                if (j != column) {
                    submatrix.data[sub_i][sub_j++] = data[i][j];
                }
            }
            sub_i++;
        }
    }
    return pow(-1, row + column) * submatrix.determinant();
}

Matrix& Matrix::inv() {
    if (n_row != n_column) {
        printf("Matrix must be square to find its inverse.\n");
		exit(EXIT_FAILURE);
    }
    double det = determinant();
    if (std::abs(det) < 1e-9) {
        printf("Matrix is singular (determinant is zero), cannot find inverse.\n");
		exit(EXIT_FAILURE);
    }

    Matrix cofactorMatrix =  zeros(n_row, n_column);
    for (int i = 0; i < n_row; ++i) {
        for (int j = 0; j < n_column; ++j) {
            cofactorMatrix.data[i][j] = cofactor(i, j);
        }
    }

    Matrix adjugateMatrix = cofactorMatrix.transpose();
    double invDet = 1.0 / det;
	Matrix& result =  zeros(n_row, n_column);
    for (int i = 0; i < n_row; ++i) {
        for (int j = 0; j < n_column; ++j) {
            result.data[i][j] = adjugateMatrix.data[i][j] * invDet;
        }
    }
    return result;
}
// Sobrecarga del operador + para la suma de un escalar a la matriz
Matrix& Matrix::operator+(double scalar) {
    Matrix* result = new Matrix(n_row, n_column);
    for (int i = 1; i <= n_row; ++i) {
        for (int j = 1; j <= n_column; ++j) {
            (*result)(i, j) = (*this)(i, j) + scalar;
        }
    }
    return *result;
}

// Sobrecarga del operador - para la resta de un escalar a la matriz
Matrix& Matrix::operator-(double scalar) {
    Matrix* result = new Matrix(n_row, n_column);
    for (int i = 1; i <= n_row; ++i) {
        for (int j = 1; j <= n_column; ++j) {
            (*result)(i, j) = (*this)(i, j) - scalar;
        }
    }
    return *result;
}

// Sobrecarga del operador * para la multiplicación de la matriz por un escalar
Matrix& Matrix::operator*(double scalar) {
    Matrix* result = new Matrix(n_row, n_column);
    for (int i = 1; i <= n_row; ++i) {
        for (int j = 1; j <= n_column; ++j) {
            (*result)(i, j) = (*this)(i, j) * scalar;
        }
    }
    return *result;
}

// Sobrecarga del operador / para la división de la matriz por un escalar
Matrix& Matrix::operator/(double scalar) {
    if (scalar == 0.0) {
        cout << "Matrix division: division by zero" << endl;
        exit(EXIT_FAILURE);
    }
    Matrix* result = new Matrix(n_row, n_column);
    for (int i = 1; i <= n_row; ++i) {
        for (int j = 1; j <= n_column; ++j) {
            (*result)(i, j) = (*this)(i, j) / scalar;
        }
    }
    return *result;
}

// Constructor de matriz identidad de tamaño n
 Matrix::Matrix(int n) {
    n_row = n;
    n_column = n;
    data = new double*[n];
    if (!data) {
        cout << "Matrix: error al asignar memoria" << endl;
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n; ++i) {
        data[i] = new double[n];
        if (!data[i]) {
            cout << "Matrix: error al asignar memoria para la fila " << i << endl;
            exit(EXIT_FAILURE);
        }
        for (int j = 0; j < n; ++j) {
            data[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

}

// Sobrecarga del operador () para acceso a un elemento con un solo índice (como en un vector)
double& Matrix::operator()(int n) {
    if (n <= 0 || n > n_row * n_column) {
        cout << "Matrix operator(): índice fuera de rango" << endl;
        exit(EXIT_FAILURE);
    }
    int row = (n - 1) / n_column + 1;
    int col = (n - 1) % n_column + 1;
    return (*this)(row, col);
}



// Función para crear una matriz de ceros de tamaño n x n
Matrix& zeros(int n) {
    Matrix* result = new Matrix(n, n);
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            (*result)(i, j) = 0.0;
        }
    }
    return *result;
}

// Función para calcular la norma de la matriz (Frobenius)
double Matrix::norm() {
    double sum = 0.0;
    for (int i = 1; i <= n_row; ++i) {
        for (int j = 1; j <= n_column; ++j) {
            sum += (*this)(i, j) * (*this)(i, j);
        }
    }
    return sqrt(sum);
}

// Función para calcular el producto punto de dos matrices (consideradas como vectores)
double Matrix::dot( Matrix& m)  {
    if (n_row != m.n_row || n_column != m.n_column) {
        cout << "Matrix dot: las matrices deben tener las mismas dimensiones" << endl;
        exit(EXIT_FAILURE);
    }
    double result = 0.0;
    for (int i = 1; i <= n_row; ++i) {
        for (int j = 1; j <= n_column; ++j) {
            result += (*this)(i, j) * m(i, j);
        }
    }
    return result;
}

// Función para calcular el producto cruz de dos matrices (consideradas como vectores 3D)
Matrix& Matrix::cross( Matrix& m)  {
    if (n_row != 1 || n_column != 3 || m.n_row != 1 || m.n_column != 3) {
        cout << "Matrix cross: las matrices deben ser vectores fila 3D" << endl;
        exit(EXIT_FAILURE);
    }

    Matrix& result=zeros(1,3);
    result(1, 1) = (*this)(1, 2) * m(1, 3) - (*this)(1, 3) * m(1, 2);
    result(1, 2) = (*this)(1, 3) * m(1, 1) - (*this)(1, 1) * m(1, 3);
    result(1, 3) = (*this)(1, 1) * m(1, 2) - (*this)(1, 2) * m(1, 1);
    return result;
}

/// Función para extraer un vector de la matriz
Matrix& Matrix::extract_vector(int n, int m) {
	Matrix& result=zeros(1,m-n+1);
	for (int i=1;i<=m-n+1;i++){
		result(i)=(*this)(n-1+i);
	}
return result;
    }

// Función para unir un vector a la matriz
Matrix& union_vector( Matrix v, Matrix v2) {
	if (v.n_row!=v2.n_row || v.n_row!=1) {
		cout << "Error en union_vector: v o v2 no es un vector" << endl;
		exit(EXIT_FAILURE);
	}
	Matrix& result = zeros(1, v.n_column+v2.n_column);
	int i;
	for(i=1; i<=v.n_column; i++)  {
		result(i)=v(i);
	}
	
	for(int j = 1; j <= v2.n_column; j++) {
		result(i)=v2(j);
		i++;
	}
    
    return result;
}

// Función para extraer una fila de la matriz
Matrix& Matrix::extract_row(int row)  {
    if (row <= 0 || row > n_row) {
        cout << "Matrix extract_row: fila fuera de rango" << endl;
        exit(EXIT_FAILURE);
    }
    Matrix* result = new Matrix(1, n_column);
    for (int j = 1; j <= n_column; ++j) {
        (*result)(1, j) = (*this)(row, j);
    }
    *this = *result;
    return *this;
}

// Función para extraer una columna de la matriz
Matrix& Matrix::extract_column(int column)  {
    if (column <= 0 || column > n_column) {
        cout << "Matrix extract_column: columna fuera de rango" << endl;
        exit(EXIT_FAILURE);
    }
    Matrix* result = new Matrix(n_row, 1);
    for (int i = 1; i <= n_row; ++i) {
        (*result)(i, 1) = (*this)(i, column);
    }
    *this = *result;
    return *this;
}

// Función para asignar una fila a la matriz
Matrix& Matrix::assign_row(int row, Matrix v) {
    if (row <= 0 || row > n_row) {
        cout << "Matrix assign_row: fila fuera de rango" << endl;
        exit(EXIT_FAILURE);
    }
    if (v.n_row != 1 || v.n_column != n_column) {
        cout << "Matrix assign_row: el tamaño del vector no coincide con el número de columnas" << endl;
        exit(EXIT_FAILURE);
    }
    for (int j = 1; j <= n_column; ++j) {
        (*this)(row, j) = v(1, j);
    }
    return *this;
}

// Función para asignar una columna a la matriz
Matrix& Matrix::assign_column(int column, Matrix v) {
    if (column <= 0 || column > n_column) {
        cout << "Matrix assign_column: columna fuera de rango" << endl;
        exit(EXIT_FAILURE);
    }
    if (v.n_column != 1 || v.n_row != n_row) {
        cout << "Matrix assign_column: el tamaño del vector no coincide con el número de filas" << endl;
        exit(EXIT_FAILURE);
    }
    for (int i = 1; i <= n_row; ++i) {
        (*this)(i, column) = v(i, 1);
    }
    return *this;
}

