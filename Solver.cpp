#include "Solver.h"
#include <cassert>

Solver::Solver(const Matrix& any, const Matrix& b) // �����������.
{
	assert((any.get_cSize() == any.get_rSize()) && "ERROR_MATRIX_IS_NOT_SQUARE");
	assert((any.get_cSize() == 0) && "ERROR_MATRIX_IS_EMPTY");

	assert((b.get_cSize() == 0) && "ERROR_MATRIX_IS_EMPTY");

	this->A = any;
	this->B = b;
}


Matrix Solver::Cramer() // ����� �������.
{
	Matrix x_Cramer(A.get_rSize(), 1);
	const double det = A.det();
	if(det != 0) // ���� ������������ �� ����� 0, �� ���������� ������� ����, � ��� ������������
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
		������ ������������ A (Det(A))
		������ ������������ ��� ������ ������� � ���������� �������� (Det(A, b, i)) i - ����� ������ ������� A, ������� ������ �� b
		�������� ����� xi (������������ i / ������������ A) � ��� �������� ��������� � ������� x_Cramer;
		*/
	}
	return x_Cramer;
} 