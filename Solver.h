#pragma once

#include "Matrix.h"
#include "Decomposition.h"

class Solver
{
private:
	// 0) Values:
	// Матрица системы (матрица А [см. конспект])   || Эти два параметра можно определить в классе Task 
	// Правая часть (вектор b [см. конспект])		|| И Асоциировать класс Task c классом Solver.
	//
	// Такой подход позволит избавить нас от необходимости хранить декомпозицию или систему в солвере.
	// Можно также добавить флаг, хранящий способ задания матрицы системы или возможность инициализации объектом разложения.
	// А можно просто использовать объект разложения только в случае вызова солвера LU разложением.

	Matrix A; // оставить вот это и забить на остальное.
	Matrix B;
	Decomposition LU;

public:
	Solver() = delete;
	Solver(std::string way);

	Matrix Cramer();
	const Matrix FindXWithLU() const; // Найим х через Гаусса

	Matrix get_x_exist(const unsigned int rown);
	void print(const Matrix& Any);
	Matrix read(std::string fullway2data);
};