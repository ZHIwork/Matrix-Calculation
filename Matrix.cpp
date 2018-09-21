//
// Created by ZHI on 7/26/2018.
//

#include "Matrix.h"
typedef vector< vector<long double> >matrix;

Matrix::Matrix(){};


matrix Matrix::Sum(const matrix& a, const matrix& b)
{
    size_t nrows = a.size();
    size_t ncols = a[0].size();
    matrix c(nrows, vector<long double>(ncols));
    for (int i = 0; i < nrows; ++i)
    {
        for (int j = 0; j < ncols; ++j)
        {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    return c;
};


matrix Matrix::Subtract(const matrix& a, const matrix& b)
{
    size_t nrows = a.size();
    size_t ncols = a[0].size();
    matrix c(nrows, vector<long double>(ncols));
    for (int i = 0; i < nrows; ++i)
    {
        for (int j = 0; j < ncols; ++j)
        {
            c[i][j] = a[i][j] - b[i][j];
        }
    }
    return c;
};


matrix Matrix::Multiplication(const matrix& a, const matrix& b)
{
    size_t n = a.size();
    size_t m = a[0].size();
    size_t p = b[0].size();
    matrix c(n, vector<long double>(p));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            long double sum = 0;
            for (int k = 0; k < m; ++k) {
                sum = sum + a[i][k]* b[k][j];
            }
            c[i][j] = sum;
        }
    }
    return c;
};

long double Matrix::SpecialMultiply(const matrix& a, const matrix& b)//a is 1*n; b is n*1
{
    size_t n = a.size();
    size_t m = a[0].size();
    size_t p = b[0].size();
    matrix c(n, vector<long double>(p));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            long double sum = 0;
            for (int k = 0; k < m; ++k) {
                sum = sum + a[i][k]* b[k][j];
            }
            c[i][j] = sum;
        }
    }
    long double value = c[0][0];
    return value;
};


matrix Matrix::Transpose(const matrix& m)
{
    size_t nrows = m.size();
    size_t ncols = m[0].size();

    //resize transpose matrix
    matrix n(ncols, vector<long double>(nrows));

    for (int i = 0; i < nrows; ++i)
    {
        for (int j = 0; j < ncols; ++j)
        {
            n[j][i] = m[i][j];
        }
    }

    return n;
};


matrix Matrix::Inverse(const matrix& m)
{
    size_t n = m.size();
    matrix inverse(n, vector<long double>(n));
    matrix adj(n, vector<long double>(n));
    // Find determinant
    long double det = determinant(m, n);

    // Find adjoint
    adj = Adjoint(m);

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            inverse[i][j] = adj[i][j]/det;

    return inverse;
}

void Matrix::getCofactor(const matrix& m, matrix& temp, int p, int q, size_t n)
{
    int i = 0, j = 0;
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q)
            {
                temp[i][j++] = m[row][col];

                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
};


long double Matrix::determinant(const matrix& m, size_t n)
{
    long double D = 0; // Initialize result

    //  Base case : if matrix contains single element
    if (n == 1)
        return m[0][0];

    matrix temp(n, vector<long double>(n)); // To store cofactors

    long double sign = 1;  // To store sign multiplier

    // Iterate for each element of first row
    for (int f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f]
        getCofactor(m, temp, 0, f, n);
        D += sign * m[0][f] * determinant(temp, n - 1);

        // terms are to be added with alternate sign
        sign = -sign;
    }

    return D;
};


matrix Matrix::Adjoint(const matrix& m)
{
    size_t n = m.size();

    // temp is used to store cofactors of A[][]
    long double sign;
    matrix temp(n, vector<long double>(n));
    matrix adj(n, vector<long double>(n));

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            // Get cofactor of A[i][j]
            getCofactor(m, temp, i, j, n);

            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = sign * determinant(temp, n-1);
        }
    }
    return adj;
}

matrix Matrix::NumberMultiply(const matrix& m, long double n)
{
    size_t nrow = m.size();
    size_t ncols = m[0].size();

    matrix temp(nrow, vector<long double>(ncols));
    for(int i=0; i<nrow; ++i)
    {
        for(int j=0; j<ncols; ++j) {
            temp[i][j] = n*m[i][j];
        }
    }
    return temp;
};


bool Matrix::ZeroMatrix(const matrix& m)
{
    size_t nrow = m.size();
    size_t ncols = m[0].size();
    for(int i=0; i<nrow; ++i)
    {
        for(int j=0; j<ncols; ++j) {
            if(m[i][j] !=0 && m[i][j] > 1e-6)
                return false;
            else
                continue;
        }
    }
    return true;
};
