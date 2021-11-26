#include "Decomposition.h"
#include <cassert>


Decomposition::Decomposition() {
    size = 0;
    LU = Matrix();
}

// 1) Constructors:
Decomposition::Decomposition(Matrix& any)
{
    // 0. Checking of sizes. If matrix isn't square then error out!
    assert((any.get_cSize() == any.get_rSize()) && "ERROR_MATRIX_IS_NOT_SQUARE");
    assert((any.get_cSize() != 0) && "ERROR_MATRIX_IS_EMPTY");

    // 1. The data is set there:
    this->size = any.get_cSize();
    LU = Matrix();
    this->A = any;

    LU = Matrix(any.get_rSize(), any.get_rSize());
    LU_decomposition();

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

    Matrix L = get_L();
    return L.at(row, col);
}

const double Decomposition::get_elemU(unsigned int row, unsigned int col) const
{
    // 0. Checking of the indexes!
    assert(((row < this->size) && (col < this->size)) && "ERROR_MATRIX_INDEX_IS_OUT_SIZE");

    Matrix U = get_U();

    return U.at(row, col);
}

const double Decomposition::get_size() const
{
    return this->size;
}

const void Decomposition::LU_decomposition() {

    //if (A.norm() != 0) {
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
    //}
}
const Matrix Decomposition::get_L() const
{
    Matrix L(size, size);
    for (size_t i = 0; i < size; i++)
    {
        for (size_t k = 0; k < size; k++)
        {
            if (i > k) { L.at(i, k) = LU.at(i, k); continue; }
            if (i == k) { L.at(i, k) = 1; continue; }
            if (i < k) { L.at(i, k) = 0; }
        }
    }
    return L;
}
const Matrix Decomposition::get_U() const
{
    Matrix U(size, size);
    for (size_t i = 0; i < size; i++)
    {
        for (size_t k = 0; k < size; k++)
        {
            if (i <= k) { U.at(i, k) = LU.at(i, k); continue; }
            if (i > k) { U.at(i, k) = 0; }
        }
    }
    return U;
}