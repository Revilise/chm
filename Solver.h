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
public:
	
};