#pragma once

#include "Matrix.h"
#include "Decomposition.h"

class Solver
{
private:
	// 0) Values:
	// ������� ������� (������� � [��. ��������])   || ��� ��� ��������� ����� ���������� � ������ Task 
	// ������ ����� (������ b [��. ��������])		|| � ������������ ����� Task c ������� Solver.
	//
	// ����� ������ �������� �������� ��� �� ������������� ������� ������������ ��� ������� � �������.
	// ����� ����� �������� ����, �������� ������ ������� ������� ������� ��� ����������� ������������� �������� ����������.
	// � ����� ������ ������������ ������ ���������� ������ � ������ ������ ������� LU �����������.

	Matrix A; // �������� ��� ��� � ������ �� ���������.
	Matrix B;
	Decomposition LU;

public:
	Solver() = delete;
	Solver(std::string way);

	Matrix Cramer();
	const Matrix FindXWithLU() const; // ����� � ����� ������

	Matrix get_x_exist(const unsigned int rown);
	void print(const Matrix& Any);
	Matrix read(std::string fullway2data);
};