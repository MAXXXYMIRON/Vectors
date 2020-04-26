#pragma once
#include<iostream>
#include <vector>
#include <math.h>
#include <ctime>

typedef std::vector<std::vector<double>> matrix;

class Matrix
{
private:
	matrix Matr;
public:
	//Конструкторы класса
	Matrix();
	~Matrix();
	Matrix(unsigned Row, unsigned Col);

	enum ERRORS
	{
		IndexOutsideMatrix = 1, //Индекс вне границ матрицы
		MatricesNotEqual = 2, //Порядки матриц не равны
		LineNotEqualColoumn = 3, //Строка матрицы не равна стобцу другой матрицы
		DivByZero = 4, //Деление на ноль
		OffLineIndex = 5, //Индекс для доступа к строке, вне границ строк матрицы
		NotRectMatrix = 6, //Матрица не прямоугольная
		NotSquareMatric = 7, //Матрица не квадратная
		DeterminateIsZero = 8, //Определитель равен нулю
		OddElementsMatrix = 9 //Матрица с нечетным количеством эл.
	};

	//Задать и получить матрицу
	matrix GetMatrix();
	void SetMatrix(matrix Matr1);
	void Set(std::vector<std::vector<float>> &matr);

	//Получить или установить эл. матрицы
	double GetElement(unsigned n, unsigned m);
	void SetElement(unsigned n, unsigned m, double Value);

	//Задание значений эл. матрицы
	//Все эл. одним значением
	void FixValue(double Value);
	//Верхняя и нижняя грань рандомных значений для эл. матрицы
	//Все эл. разные
	void RandValue(int FirstLim, int LastLim);

	//Присвоит полученое значение новой матрице
	Matrix operator + (const Matrix& Matr2) const;
	Matrix operator - (const Matrix& Matr2) const;
	Matrix operator * (const Matrix& Matr2) const;

	//Присвоит полученое значение новой матрице
	Matrix Add(double Value) const;
	Matrix Sub(double Value) const;
	Matrix Mul(double Value) const;
	Matrix Div(double Value) const;

	//Альтернатива вшеопределённым
	Matrix operator - () const; //Домножит на -1
	Matrix operator - (double Value) const;
	Matrix operator + (double Value) const;
	Matrix operator * (double Value) const;
	Matrix operator / (double Value) const;
	Matrix operator ++ () const; //+1
	Matrix operator -- () const; //-1

	//Текущая матрица преобразуется в транспонированную или обратную
	void Transpose();
	void Inverse();

	//Приведение к обратной методом разбиения на блоки (Нужно для Вычмата)
	void InverseBloking();

	//Вернет определитель матрицы, если та - квадратная
	double Detr() const;

	//Получить или установить строку матрицы
	std::vector<double> GetString(unsigned Row);
	void SetString(unsigned Row, std::vector<double> Str);

	//Полученное значение запишет в текущую матрицу
	void operator += (const Matrix& Matr2);
	void operator -= (const Matrix& Matr2);
	void operator *= (const Matrix& Matr2);

	//Полученное значение запишет в текущую матрицу
	void operator += (double Value);
	void operator -= (double Value);
	void operator *= (double Value);
	void operator /= (double Value);

private:
	//Проверка на не выход индекса за границы матрицы
	bool WithinBorders(unsigned IndexRow, unsigned IndexCol) const;
	//Проверка на равенстнво матриц
	bool EqualMatrix(const matrix& Matr1, const matrix& Matr2) const;
	//Проверка на равенство строки столбцу
	bool EqualRowCol(const matrix& Matr1, const matrix& Matr2) const;
	//Проверка на не ступенчатость матрицы
	bool RectMatrix(const matrix& Matr1) const;
	//Проверка матрицы на квадратность
	bool SquareMatrix(const matrix& Matr1) const;
	//Проверка на четное кол-во эл. в матрице
	bool EvenNumbOfEl(const matrix& Matr1) const;

	//Определитель матрицы 2 на 2
	double Detr(const matrix& M) const;
	//Нахождение определителей миноров матрицы
	double Detr(const matrix& M, unsigned Row, unsigned Col) const;

	//Проверка на повторное использование строки, для нахождения обратной матрицы
	bool RepeatUsedLine(const std::vector<double>& ArrayIndexLine, unsigned IndexLine);
	//Деление строки матрицы на число
	std::vector<double> LineDivNumber(const std::vector<double>& Line, double Number);
	//Вычитание строк, одна из которых домножена на Index эл. первой
	std::vector<double> SubLine(const std::vector<double>& Line1, const std::vector<double>& Line2, unsigned Index);
	//Выравнивание еденичной матрицы, т.к. еденицы могут находится не на главной диагонали
	void Unit(matrix& Matr1);

	//Вернет определённую четверть матрицы
	matrix Block(unsigned IndexRow, unsigned IndexCol) const;
	//Запишет указанную матрицу в матрицу экземпляра
	void Block(unsigned IndexRow, unsigned IndexCol, matrix TempMatr);
};
