#include "Solver.h"
#include "Decomposition.h"

#include <cassert>

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>


Solver::Solver(std::string way) // Конструктор.
{
	this->A = read(way);

	//assert((A.get_cSize() == A.get_rSize()) && "ERROR_MATRIX_IS_NOT_SQUARE");
	//assert((A.get_cSize() == 0) && "ERROR_MATRIX_IS_EMPTY");

	Matrix x_exist = get_x_exist(A.get_cSize());

	this->B = A * x_exist;
	this->LU = Decomposition(A);

	std::printf("\nreal root =\n");
	print(x_exist);

	std::printf("\nreal x =\n");
	print(x_exist);

	std::printf("\nDecomposition =\n");
	print(FindXWithLU());

	std::printf("\nCrumer =\n");
	print(Cramer());
}

Matrix Solver::get_x_exist(const unsigned int rown) {
	Matrix x_exist(rown, 1);

	for (size_t row = 0; row < rown; row++) {
		x_exist.at(row, 0) = std::sin(row + 1);
	}
	return x_exist;
}

Matrix Solver::read(std::string fullway2data)
{
	std::ifstream inputfile;
	inputfile.open(fullway2data);

	Matrix Res;

	if (inputfile.is_open())
	{
		std::string buff_s;
		double buff_d;
		std::vector <std::vector<double>> buff_data;
		std::vector <double> buff_data_row;

		while (getline(inputfile, buff_s))
		{
			std::istringstream buff_ss(buff_s);

			while (buff_ss >> buff_d)
			{
				buff_data_row.push_back(buff_d);
			}

			buff_data.push_back(buff_data_row);
			buff_data_row.clear();
		}

		Res = Matrix(buff_data.size(), buff_data.at(0).size());

		for (size_t row = 0; row < Res.get_rSize(); row++)
		{
			assert((buff_data.at(row).size() == Res.get_cSize()) && "ERROR_COPIED_MATRIX_COLUMNS_SIZES_SHOULD_BE_EQUAL");

			if (buff_data.at(row).size() != Res.get_cSize())
			{
				std::cout << "ERROR: copying matrix is failed! Process was stopped!" << std::endl;

				return Res;
			}

			for (size_t col = 0; col < Res.get_cSize(); col++)
			{
				Res.at(row, col, buff_data.at(row).at(col));
			}
		}
	}
	else
	{
		std::cout << "ERROR: copying matrix is failed! File isn't opened!" << std::endl;
	}

	return Res;
}

void Solver::print(const Matrix& Any)
{
	int precicion = 16; // TODO: изменить на %;
	if ((Any.get_rSize() == 0) || (Any.get_cSize() == 0))
	{
		std::cout << "WARNING: printed matrix is empty!" << std::endl;
	}

	for (size_t i = 0; i < Any.get_rSize(); i++)
	{
		for (size_t j = 0; j < Any.get_cSize(); j++)
		{
			std::cout << std::setprecision(precicion) << std::scientific << Any.at(i, j) << "		";
		}
		std::cout << std::endl;
	}
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

const Matrix Solver::FindXWithLU() const {

	int size = LU.get_size();

	Matrix y(size, 1);
	y.at(0, 0) = B.at(0, 0);

	for (size_t k = 1; k < size; k++)
	{
		y.at(k, 0) = B.at(k, 0);
		for (size_t p = 0; p < k; p++)
		{
			y.at(k, 0) -= LU.get_elemL(k, p) * y.at(p, 0);
		}
	}

	Matrix x(size, 1);
	int m = size - 1;
	x.at(m, 0) = y.at(m, 0) / LU.get_elemU(m, m);

	for (size_t k = size - 2; k >= 0 && k < size; k--)
	{
		x.at(k, 0) = y.at(k, 0);
		for (size_t p = k + 1; p < m; p++)
		{
			x.at(k, 0) -= LU.get_elemU(k, p) * x.at(p, 0);
		}
		x.at(k, 0) *= 1 / LU.get_elemU(k, k);
	}
	return x;
}
