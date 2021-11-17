//#define NDEBUG    // инструкция дл препроцессора отключающая макрос assert

#include "Matrix.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include <sstream>
#include <iomanip>

Matrix read(std::string fullway2data)
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

Matrix get_x_exist(const unsigned int rown) {
	Matrix x_exist(rown, 1);

	for (size_t row = 0; row < rown; row++) {
		x_exist.at(row, 0) = std::sin(row + 1);
	}
	return x_exist;
}

void print(const Matrix& Any)
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

int main()
{
    // В релизной версии программы необходимо отключить все макросы проверки (assert).
    // Все провеки, связанные с загрузкой матрицы из файла и корректности размеров,
    // должны быть явно реализованы в main. Релизная версия программы всегда должна завершаться
    // с кодом 0. Другой код, свидетельствующий об ошибке,в отлаженной программе возникать не должен.
    // 
    // Отключить все макросы необходимо добавив дерективу дял препроцессора в этом файле.
    // Добавление директивы пропусти добавление всех строк содержащих assert.

    std::string way = "I:/чм/9.txt";
    
    Matrix A = read(way); // ошибка тут
    Matrix x_exist = get_x_exist(A.get_cSize());
    Matrix b = A * x_exist;

	std::printf("x_exist =\n");
	print(x_exist);

	std::printf("A =\n");
	print(A);

	//std::printf("b =\n");
	//print(b);

    return 0;
}
