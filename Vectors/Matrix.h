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
	//������������ ������
	Matrix();
	~Matrix();
	Matrix(unsigned Row, unsigned Col);

	enum ERRORS
	{
		IndexOutsideMatrix = 1, //������ ��� ������ �������
		MatricesNotEqual = 2, //������� ������ �� �����
		LineNotEqualColoumn = 3, //������ ������� �� ����� ������ ������ �������
		DivByZero = 4, //������� �� ����
		OffLineIndex = 5, //������ ��� ������� � ������, ��� ������ ����� �������
		NotRectMatrix = 6, //������� �� �������������
		NotSquareMatric = 7, //������� �� ����������
		DeterminateIsZero = 8, //������������ ����� ����
		OddElementsMatrix = 9 //������� � �������� ����������� ��.
	};

	//������ � �������� �������
	matrix GetMatrix();
	void SetMatrix(matrix Matr1);
	void Set(std::vector<std::vector<float>> &matr);

	//�������� ��� ���������� ��. �������
	double GetElement(unsigned n, unsigned m);
	void SetElement(unsigned n, unsigned m, double Value);

	//������� �������� ��. �������
	//��� ��. ����� ���������
	void FixValue(double Value);
	//������� � ������ ����� ��������� �������� ��� ��. �������
	//��� ��. ������
	void RandValue(int FirstLim, int LastLim);

	//�������� ��������� �������� ����� �������
	Matrix operator + (const Matrix& Matr2) const;
	Matrix operator - (const Matrix& Matr2) const;
	Matrix operator * (const Matrix& Matr2) const;

	//�������� ��������� �������� ����� �������
	Matrix Add(double Value) const;
	Matrix Sub(double Value) const;
	Matrix Mul(double Value) const;
	Matrix Div(double Value) const;

	//������������ ��������������
	Matrix operator - () const; //�������� �� -1
	Matrix operator - (double Value) const;
	Matrix operator + (double Value) const;
	Matrix operator * (double Value) const;
	Matrix operator / (double Value) const;
	Matrix operator ++ () const; //+1
	Matrix operator -- () const; //-1

	//������� ������� ������������� � ����������������� ��� ��������
	void Transpose();
	void Inverse();

	//���������� � �������� ������� ��������� �� ����� (����� ��� �������)
	void InverseBloking();

	//������ ������������ �������, ���� �� - ����������
	double Detr() const;

	//�������� ��� ���������� ������ �������
	std::vector<double> GetString(unsigned Row);
	void SetString(unsigned Row, std::vector<double> Str);

	//���������� �������� ������� � ������� �������
	void operator += (const Matrix& Matr2);
	void operator -= (const Matrix& Matr2);
	void operator *= (const Matrix& Matr2);

	//���������� �������� ������� � ������� �������
	void operator += (double Value);
	void operator -= (double Value);
	void operator *= (double Value);
	void operator /= (double Value);

private:
	//�������� �� �� ����� ������� �� ������� �������
	bool WithinBorders(unsigned IndexRow, unsigned IndexCol) const;
	//�������� �� ���������� ������
	bool EqualMatrix(const matrix& Matr1, const matrix& Matr2) const;
	//�������� �� ��������� ������ �������
	bool EqualRowCol(const matrix& Matr1, const matrix& Matr2) const;
	//�������� �� �� ������������� �������
	bool RectMatrix(const matrix& Matr1) const;
	//�������� ������� �� ������������
	bool SquareMatrix(const matrix& Matr1) const;
	//�������� �� ������ ���-�� ��. � �������
	bool EvenNumbOfEl(const matrix& Matr1) const;

	//������������ ������� 2 �� 2
	double Detr(const matrix& M) const;
	//���������� ������������� ������� �������
	double Detr(const matrix& M, unsigned Row, unsigned Col) const;

	//�������� �� ��������� ������������� ������, ��� ���������� �������� �������
	bool RepeatUsedLine(const std::vector<double>& ArrayIndexLine, unsigned IndexLine);
	//������� ������ ������� �� �����
	std::vector<double> LineDivNumber(const std::vector<double>& Line, double Number);
	//��������� �����, ���� �� ������� ��������� �� Index ��. ������
	std::vector<double> SubLine(const std::vector<double>& Line1, const std::vector<double>& Line2, unsigned Index);
	//������������ ��������� �������, �.�. ������� ����� ��������� �� �� ������� ���������
	void Unit(matrix& Matr1);

	//������ ����������� �������� �������
	matrix Block(unsigned IndexRow, unsigned IndexCol) const;
	//������� ��������� ������� � ������� ����������
	void Block(unsigned IndexRow, unsigned IndexCol, matrix TempMatr);
};
