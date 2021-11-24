#include "Decomposition.h"
#include <cassert>


// 1) Constructors:
Decomposition::Decomposition(Matrix& anyA, Matrix& anyB)
{
    // 0. Checking of sizes. If matrix isn't square then error out!
    assert((anyA.get_cSize() == anyA.get_rSize()) && "ERROR_MATRIX_IS_NOT_SQUARE");
    assert((anyA.get_cSize() == 0) && "ERROR_MATRIX_IS_EMPTY");

    // 1. The data is set there:
    this->size = anyA.get_cSize();
    this->A = anyA;
    this->B = anyB;

   // this->values = { 10.0, -0.1, 1.0, 10.1 };   //  <- Заменить на LU разложение
    // По умолчанию: (L\U) = { 10.0, -0.1, 1.0, 10.1 } - LU разложение для матрицы А = {10.0, -1.0, 1.0, 10.0}

}

// 2) Destructor:
Decomposition::~Decomposition()
{
    this->values.clear();
    this->values.shrink_to_fit();
}

// 3) Geters and seters:
const double Decomposition::get_elemL(unsigned int row, unsigned int col) const
{
    // 0. Checking of the indexes!
    assert(((row < this->size) && (col < this->size)) && "ERROR_MATRIX_INDEX_IS_OUT_SIZE");

    return 0.0;
}

const double Decomposition::get_elemU(unsigned int row, unsigned int col) const
{
    // 0. Checking of the indexes!
    assert(((row < this->size) && (col < this->size)) && "ERROR_MATRIX_INDEX_IS_OUT_SIZE");

    return 0.0;
}

const double Decomposition::get_size() const
{
    return this->size;
}

const Matrix Decomposition::get_L() const
{
    return Matrix();
}

const Matrix Decomposition::get_U() const
{
    return Matrix();
}

const Matrix Decomposition::LU_decomposition() {
    Matrix LU(A.get_rSize(), A.get_rSize());

    if (A.norm() != 0) {
        for (int k = 0; k < A.get_rSize(); k++) {
            // вычисление первой строки матрицы U
            LU.at(0, k) = A.at(k, 0);
        }
        for (int k = 1; k < A.get_rSize(); k++) {
            // вычисление первого столбца матрицы L
            LU.at(k, 0) = A.at(k, 0) / A.at(0, 0);
        }
        // вычисление остальной матрицы LU
        for (int i = 1; i < A.get_cSize(); i++) {

            // U
            for (int k = i; k < A.get_cSize(); k++) {
                double sum = 0.0;
                for (int p = 0; p < i; p++) {
                    sum += LU.at(i, p) * LU.at(p, k);
                }
                LU.at(i, k) = A.at(i, k) - sum;
            }
            // L
            for (int k = i + 1; k < A.get_cSize(); k++) {
                double sum = 0.0;
                for (int p = 0; p < i - 1; p++) {
                    sum += LU.at(i, p) * LU.at(p, k);
                }
                LU.at(i, k) = (A.at(i, k) - sum) / LU.at(i, i);
            }
        }
        this->LU = LU;
    }
    return LU;

}

const Matrix Decomposition::Gausse() const {

    Matrix y;
    y.at(0, 0) = B.at(0, 0);

    // Ly = b
    for (size_t k = 1; k < LU.get_rSize(); k++) {
        y.at(k, 0) = B.at(k, 0);
        for (size_t p = 0; p < k; p++) {
            y.at(k, 0) -= LU.at(k, p) * y.at(p, 0);
        }
    }

    // Ux=y
    Matrix x = Matrix(LU.get_rSize(), 0);
    int m = LU.get_rSize() - 1;
    
    x.at(m, 0) = y.at(m, 0) / LU.at(m, m);

    for (size_t k = LU.get_rSize() - 2; k >= 0; k--) {
        x.at(k, 0) = y.at(k, 0);
        for (size_t p = k + 1; p < LU.get_rSize() - 1; p++) {
            x.at(k, 0) -= LU.at(k, p) * x.at(p, 0);
        }
        x.at(k, 0) *= 1 / LU.at(k, k);
    }

    return x;
}