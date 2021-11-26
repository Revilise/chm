#pragma once
//выполнять разложение квадратной и невырожденной матрицы;
//хранить матрицы L и U(разложение) в компактном виде(L\U);
//уметь извлекать матрицы L и U из(L\U);
//уметь извлекать элементы матриц L и U.
#include <iostream>
#include <vector>
#include "Matrix.h"

class Decomposition
{
private:
	// 0) Values:
	std::vector <double> values;
	unsigned int size;

	Matrix A;
	Matrix LU;

public:
	// 1) Constructors:
	Decomposition(Matrix& any);
	Decomposition();

	// 2) Destructor:
	~Decomposition();

	// 3) Geters and seters:
	const double get_elemL(unsigned int row, unsigned int col) const;
	const double get_elemU(unsigned int row, unsigned int col) const;

	//const Matrix& U(unsigned int row, unsigned int col) const;
	//const Matrix& L(unsigned int row, unsigned int col) const;

	const Matrix get_L() const; // Можно ли подставить const Matrix& (т.е. ссылку на локальный объект)?
	const Matrix get_U() const; // Можно ли подставить const Matrix& (т.е. ссылку на локальный объект)?

	const double get_size() const;

	const void LU_decomposition();
};