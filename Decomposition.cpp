#include "Decomposition.h"
#include <cassert>


// 1) Constructors:
Decomposition::Decomposition(const Matrix& any):LU   // LU можно убрать если что
{
    // 0. Checking of sizes. If matrix isn't square then error out!
    assert((any.get_cSize() == any.get_rSize()) && "ERROR_MATRIX_IS_NOT_SQUARE");
    assert((any.get_cSize() == 0) && "ERROR_MATRIX_IS_EMPTY");

    // 1. The data is set there:
    this->size = any.get_cSize();
    this->values = { 10.0, -0.1, 1.0, 10.1 };   //  <- Заменить на LU разложение
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