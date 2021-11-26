#pragma once
//��������� ���������� ���������� � ������������� �������;
//������� ������� L � U(����������) � ���������� ����(L\U);
//����� ��������� ������� L � U ��(L\U);
//����� ��������� �������� ������ L � U.
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

	const Matrix get_L() const; // ����� �� ���������� const Matrix& (�.�. ������ �� ��������� ������)?
	const Matrix get_U() const; // ����� �� ���������� const Matrix& (�.�. ������ �� ��������� ������)?

	const double get_size() const;

	const void LU_decomposition();
};