#include "Matrix.h"
#include <cassert>
#include <fstream>

// Stuff functions:
unsigned int get_row2swap(const unsigned int index_diag, const Matrix& Any)
// Ôóíêöèÿ îñóùåñòâëÿåò ïîèñê íîìåðà ñòðîêè íà êîòîðóþ íåîáõîäèìî çàìåíèòü ðàññìàòðèâàåìóþ,
// ñ íóëåâûì äèàãîíàëüíûì ýëåìåíòîì (ïîèñê íåíóëåâîãî ýëåìåíòà â ñòîëáöå ñ íîìåðîì index_diag).
// 
// Åñëè âåðíóâøååñÿ çíà÷åíèå èç ôóíêöèè index2swap ðàâíÿåòñÿ êîëè÷åñòâó ñòðê â ìàòðèöå =>
// => îïðåäåëèòåëü ìàòðèöû íîëü!
{
	unsigned int index2swap = index_diag;
	while ((index2swap < Any.get_rSize()) && (Any.at(index2swap, index_diag) == 0.0))
	{
		index2swap = index2swap + 1;
	}

	return index2swap;
};

bool swap_rows(const unsigned int index_diag, Matrix& Any)
// Ôóíêöèÿ âîçâðàùàåò ôëàã bool:
//		true, åñëè áûëà âûïîëíåíà ïåðåñòàíîâêà ñòðîê
//		false, åñëè ïåðåñòàíîâêè ñòðîê íå áûëî
{
	bool swap_flag = false;

	const unsigned int index2swap = get_row2swap(index_diag, Any);

	if (index2swap != Any.get_rSize())
	{
		// Ïåðåêëþ÷åíèå ôëàãà:
		swap_flag = true;

		// Ïåðåñòàíîâêà ñòðîê:
		double buffer_value;
		for (size_t col = index_diag; col < Any.get_cSize(); col++)
		{
			buffer_value = Any.at(index_diag, col);
			Any.at(index_diag, col) = Any.at(index2swap, col);
			Any.at(index2swap, col) = buffer_value;
		}
	}

	return swap_flag;
}

void column_reset(const unsigned int index_diag, Matrix& Any)
{
	double swaped_value;
	for (size_t row = index_diag + 1; row < Any.get_rSize(); row++)
	{
		swaped_value = Any.at(row, index_diag);
		Any.at(row, index_diag) = 0.0;

		for (size_t col = index_diag + 1; col < Any.get_cSize(); col++) // Áûëà îïå÷àòêà â èíäåêñàõ. Ñòàðòîâûé èíäåêñ áûë col = row, ò.å. êàæäûé ðàç ïðîèñõîèëî ñìåùåíèå íà + 1
		{
			//add(Any, row, col, -swaped_value * Any.get_elem(index_diag, col));
			Any.at(row, col) -= swaped_value * Any.at(index_diag, col);
		}
	}
}

// Äîáàâëåííà ôóíêöèÿ, êîòîðàÿ äåëèò ñòðîêó íà äèàãîíàëüíûé ýëåìåíò è êîððåêòèðóåò çíà÷åíèå îïðåäåëåèòåëÿ
void row_sub(const unsigned int index_diag, double& det_value, Matrix& Copy)
{
	const double value_diag = Copy.at(index_diag, index_diag);
	if (value_diag != 1.0)
	{
		det_value = det_value * value_diag;

		// Äåëåíèå ñòðîêè íà äèàãîíàëüíûé ýëåìåíò:
		for (size_t col = index_diag; col < Copy.get_cSize(); col++)
		{
			Copy.at(index_diag, col) /= value_diag;
		}
	}
}


// -1) The private geter gets a linear index:
unsigned int Matrix::get_index(unsigned int row, unsigned int col) const
{
	// n = i - 1 + (j - 1) * rown    â ñëó÷àå, åñëè i in [1, rown], à j in [1, coln]  =>
	// => return row - 1 + (col - 1) * this->rown;
	// n = i + j * rown		â ñëó÷àå, åñëè i in [0, rown-1], à j in [0, coln-1] =>
	// => return row + col * this->rown;
	// Ñìîòðè ïîäðîáíîå îïèñàíèå â exel ôàéëå çàäàíèÿ.

	//assert((col < this->coln) && "ERROR_MATRIX_INDEX_IS_OUT_SIZE"); // assert(bool = true)
	//assert((row < this->rown) && "ERROR_MATRIX_INDEX_IS_OUT_SIZE");

	return row + col * this->rown;
}

// 1) Ñonstructors:
Matrix::Matrix() : values({ 10.0, -1.0, 1.0, 10.0 }), rown(2), coln(2) {}

Matrix::Matrix(unsigned int rown, unsigned int coln) : coln(coln), rown(rown), values(coln* rown) {}

// 2) Destructior:
Matrix::~Matrix()
{
	values.clear();
	values.shrink_to_fit();
}

// 3) Geters and seters:
const unsigned int Matrix::get_rSize() const
{
	return this->rown;
}

const unsigned int Matrix::get_cSize() const
{
	return this->coln;
}

void Matrix::at(unsigned int row, unsigned int col, const double value)
{
	assert((row < this->get_rSize()) && (col < this->get_cSize()) && "KEK");
	this->at(row, col) = value;
}

double& Matrix::at(unsigned int row, unsigned int col)
{
	assert((row < this->get_rSize()) && (col < this->get_cSize()) && "KEK");
	return this->values.at(get_index(row, col));
}

const double& Matrix::at(unsigned int row, unsigned int col) const
{
	assert((row < this->get_rSize()) && (col < this->get_cSize()) && "KEK");
	return this->values.at(get_index(row, col));
}

const Matrix Matrix::set_column(const unsigned int col, const Matrix& column) const
{
	// 0. Checking of the indexes:
	assert((col < this->coln) && "ERROR_MATRIX_INDEX_IS_OUT_SIZE");
	assert((column.rown == this->rown) && "ERROR_MATRIXES_SIZES_SHOULD_BE_EQUAL");
	assert((column.coln == 1) && "ERROR_MATRIX_SHOULD_BE_A_COLUMN");
	Matrix Res = *this;
	// 1. The column is inserted there:
	Res.values.erase(Res.values.begin() + col * Res.rown, Res.values.begin() + (col + 1) * Res.rown);
	Res.values.insert(Res.values.begin() + col * Res.rown, column.values.begin(), column.values.end());

	return Res;
}

// Ôóíêöèÿ ðàññ÷¸òà îïðåäåëèòåëÿ ìåòîäîì èñêëþ÷åíèÿ Ãàóññà.
// Åñëè ïðèèñàòü - ãîòîâàÿ ôóíêöèÿ äëÿ ìåòîäà ïîäñòàíîâîê (èñêëþ÷åíèÿ)
// äëÿ ìåòîäà LU ðàçëîæåíèÿ â ñëîâåð.
// Åñëè áðàòü â äàëüíå¸øóþ ðàîòó, òî íåîáõîäèìî ïðè÷åñàòü, ðàçáèâ íà ôóíêöèè.
const double Matrix::det() const
{
	// 0. Checking of the sizes:
	assert((this->coln == this->rown) && "ERROR_MATRIX_IS_NOT_SQUARE");
	assert((this->coln != 0) && "ERROR_MATRIXES_SIZES_SHOULD_BE_NO_ZERO");

	// 1. If matrix is number:
	if ((this->coln == 1) && (this->rown == 1))
	{
		return this->at(0, 0);
	}

	// 2. If matrix is suqare:
	Matrix Copy = *this;

	double det_value = 1.0;

	unsigned int index_diag = 0;
	double value_diag;

	while (index_diag < Copy.get_rSize())
	{
		// Ïðîâåðêà íà íóëåâîé äèàãîíàëüíûé ýëåìåíò:
		if (Copy.at(index_diag, index_diag) == 0.0)
		{
			// Ïîèñê íåíóëåâîãî ýëåíòà â ñòîëáöå, ëåæàùåãî íèæå:
//**********// call swap_rows(...)
			unsigned int index2swap = index_diag;
			while ((index2swap < Copy.get_rSize()) && (Copy.at(index2swap, index_diag) == 0.0))
			{
				index2swap = index2swap + 1;
			}

			// Ïðîâåðêà èññêëþ÷åíèÿ, åñëè âñå ýëåìåíòû - íóëåâûå => det = 0
			// swap_rows(...) -> false => return det_value = 0.0;
			// swap_rows(...) -> true => det_value = det_value * (-1.0);
			if (index2swap == Copy.get_rSize())
			{
				return det_value = 0.0;
			}
			// Ïåðåñòàíîâêà ñòðîê ìåñòàìè:
			else
			{
				double buffer_value;
				for (size_t col = index_diag; col < Copy.get_cSize(); col++)
				{
					buffer_value = Copy.at(index_diag, col);
					Copy.at(index_diag, col) = Copy.at(index2swap, col);
					Copy.at(index2swap, col) = buffer_value;
				}

				det_value = det_value * (-1.0); // ò.ê. ïðè ïåðåñòàíîâå ñòðîê íåîáõîäèìî ïîìåíÿòü îïðåäåëèòåëü ìåñòàìè
			}
		}

		// Ïðîöåññ èñêëþ÷åíèÿ (áóäåò çàïóùåí, òîëüêî åñëè) det != 0
//******// call row_sub(...)
		value_diag = Copy.at(index_diag, index_diag);

		det_value = det_value * value_diag;

		// Äåëåíèå ñòðîêè íà äèàãîíàëüíûé ýëåìåíò:
		for (size_t col = index_diag; col < Copy.get_cSize(); col++)
		{
			Copy.at(index_diag, col) /= value_diag;
		}
		//******// end call row_sub(...)

				// Èñêëþ÷åíèå ýëåìåíòîâ ëåæàùèõ íèæå äèàãîíàëüíûõ:
		//******// call column_reset(...)
		double swaped_value;
		for (size_t row = index_diag + 1; row < Copy.get_rSize(); row++)
		{
			swaped_value = Copy.at(row, index_diag);
			Copy.at(row, index_diag) = 0.0;

			for (size_t col = index_diag + 1; col < Copy.get_cSize(); col++)
			{
				Copy.at(row, col) -= -swaped_value * Copy.at(index_diag, col);
			}
		}
		//******// end call column_reset(...)

		index_diag = index_diag + 1;
	}

	return det_value;
}

const double Matrix::norm() const
{
	return 0.0;
}

Matrix& Matrix::operator=(const Matrix& Any)
{
	// Èñïîëüçîâàíèå ïåðåãðóçêè îïðåàòîðà Matrix& Matrix::operator=(const Matrix& Any),
	// âîçâðàùàþùåãî ññûëêó íà îáúåêò, à íå îáúåêò, ïîçâîëÿåò âûïîëíÿòü öåïî÷êó ïðèñâîåíèé!
	//
	// Äëÿ òîãî, ÷òîáû âåðíóòü ññûëêó íà îáúåêò, îïðåäåëÿåìûé â òåëå èñïîëüçóéòå return *this;

	// 0. Ïðîâåðêà íà ñàìîïðèñâîåíèå.
	// ×òîáû íå âûïîëíÿòü ëèøíåå êîïèðîâàíèå.
	// Âîçâðàùàåò ññûëêó íà òåêùèé îáúåêò.
	if (this == &Any)
	{
		return *this;
	}

	// 1. The copying of the object values:
	this->coln = Any.coln;
	this->rown = Any.rown;
	this->values = Any.values;

	return *this;
}

Matrix operator+(const Matrix& left, const Matrix& right)
{
	// 0. Checking of the sizes:
	assert((left.get_cSize() == right.get_cSize()) && "ERROR_MATRIXES_SIZES_SHOULD_BE_EQUAL");
	assert((left.get_rSize() == right.get_rSize()) && "ERROR_MATRIXES_SIZES_SHOULD_BE_EQUAL");
	assert((left.get_cSize() != 0) && "ERROR_MATRIXES_SIZES_SHOULD_BE_NO_ZERO");
	assert((left.get_rSize() != 0) && "ERROR_MATRIXES_SIZES_SHOULD_BE_NO_ZERO");

	// 1. The matrix result is created there:
	Matrix result(left.get_rSize(), left.get_cSize());

	for (size_t j = 0; j < right.get_cSize(); j++)
	{
		for (size_t i = 0; i < right.get_rSize(); i++)
		{
			result.at(i, j) = left.at(i, j) + right.at(i, j);
		}
	}

	return result;
}

Matrix operator-(const Matrix& left, const Matrix& right)
{
	// 0. Checking of the sizes:
	assert((left.get_cSize() == right.get_cSize()) && "ERROR_MATRIXES_SIZES_SHOULD_BE_EQUAL");
	assert((left.get_rSize() == right.get_rSize()) && "ERROR_MATRIXES_SIZES_SHOULD_BE_EQUAL");
	assert((left.get_cSize() != 0) && "ERROR_MATRIXES_SIZES_SHOULD_BE_NO_ZERO");
	assert((left.get_rSize() != 0) && "ERROR_MATRIXES_SIZES_SHOULD_BE_NO_ZERO");
	Matrix result(left.get_rSize(), left.get_cSize());

	for (size_t j = 0; j < right.get_cSize(); j++)
	{
		for (size_t i = 0; i < right.get_rSize(); i++)
		{
			result.at(i, j) = left.at(i, j) - right.at(i, j);
		}
	}

	return result;
	//return Matrix();
}

Matrix operator*(const Matrix& left, const Matrix& right)
{
	
	
	// 0. Checking of the sizes:
	assert((left.get_rSize() != right.get_cSize()) && "ERROR_MATRIXES_SIZES_AT_ROW_AND_AT_COLUMN_SHOULD_BE_EQUAL");
	assert((left.get_cSize() != 0) && "ERROR_MATRIXES_SIZES_SHOULD_BE_NO_ZERO");
	assert((left.get_rSize() != 0) && "ERROR_MATRIXES_SIZES_SHOULD_BE_NO_ZERO");
	
	Matrix result(left.get_rSize(), right.get_cSize());
	for (size_t j = 0; j < left.get_rSize(); j++)
	{
		for (size_t i = 0; i < right.get_cSize(); i++)
		{

			for (size_t k = 0; k < left.get_cSize(); k++) 
			{
				result.at(j, i) += left.at(j, k) * right.at(k, i);
			}
			
		}
	}
 
	return result;

	//return Matrix();
}
