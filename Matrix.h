//
// Created by ZHI on 7/26/2018.
//

#ifndef MPI_SIO2_MATRIX_H
#define MPI_SIO2_MATRIX_H

#include <iostream>
#include <vector>
using namespace std;

class Matrix {
private:
    typedef vector< vector<long double> >matrix;
public:
    Matrix();
    matrix Sum(const matrix& a, const matrix& b);
    matrix Subtract(const matrix& a, const matrix& b);
    matrix Multiplication(const matrix& a, const matrix& b);
    long double SpecialMultiply(const matrix& a, const matrix& b);
    matrix Transpose(const matrix& m);
    matrix Inverse(const matrix& m);
    void getCofactor(const matrix& m, matrix& temp, int p, int q, size_t n);
    long double determinant(const matrix& m, size_t n);
    matrix Adjoint(const matrix& m);
    matrix NumberMultiply(const matrix& m, long double n);
    bool ZeroMatrix(const matrix& m);
};


#endif //MPI_SIO2_MATRIX_H
