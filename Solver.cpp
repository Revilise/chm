#include "Solver.h"
#include <cassert>

Solver::Solver(const Matrix& any, const Matrix& b) // Конструктор.
{
	assert((any.get_cSize() == any.get_rSize()) && "ERROR_MATRIX_IS_NOT_SQUARE");
	assert((any.get_cSize() == 0) && "ERROR_MATRIX_IS_EMPTY");

	assert((b.get_cSize() == 0) && "ERROR_MATRIX_IS_EMPTY");

	this->A = any;
	this->B = b;
}


Matrix Solver::Cramer() // Метод крамера.
{
	Matrix x_Cramer(A.get_rSize(), 1);
	const double det = A.det();
	if(det != 0) // Если оперделитель не равен 0, то существует решение СЛАУ, и оно едиснтвенное
	{
		std::vector <double> values;
		for(int i = 0; i < A.get_cSize(); i++)
		{
			Matrix insult = A.set_column(i, B);
			values.push_back(insult.det());
		}
		for(int i = 0; i < x_Cramer.get_rSize(); i++)
		{
			x_Cramer.at(i, 0) = values[i] / det;
		}
		/*
		Нахожу определитель A (Det(A))
		Нахожу определитель для каждой матрицы с замененным столбцом (Det(A, b, i)) i - номер стобца матрицы A, который заменю на b
		Вычисляю корни xi (Определитель i / Определитель A) и все значения определяю в матрицу x_Cramer;
		*/
	}
	return x_Cramer;
} 