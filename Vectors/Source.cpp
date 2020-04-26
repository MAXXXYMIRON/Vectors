#include<iostream>
#include<cmath>
#include<vector>
#include "Matrix.h"

using namespace std;

vector<vector<float>> E(int size)
{
	vector<vector<float>> E;
	E.resize(size);
	for (int i = 0; i < size; i++)
	{
		E[i].resize(size);
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			E[i][j] = (i != j) ? 0 : 1;
		}
	}
	return E;
}

void MUL(vector<vector<float>> &A, vector<vector<float>> &B, Matrix &Inv, int size)
{
	vector<vector<float>> temp;
	temp.resize(size);
	for (int i = 0; i < size; i++)
	{
		temp[i].resize(size);
		for (int j = 0; j < size; j++)
			temp[i][j] = 0;
	}


	for (unsigned i = 0; i < size; i++)
		for (unsigned j = 0; j < size; j++)
			for (unsigned k = 0; k < size; k++)
				temp[i][j] += Inv.GetElement(i, k) * A[k][j];

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			A[i][j] = 0;
	}

	for (unsigned i = 0; i < size; i++)
		for (unsigned j = 0; j < size; j++)
			for (unsigned k = 0; k < size; k++)
				A[i][j] += temp[i][k] * B[k][j];
}

vector<vector<float>> MUL(vector<vector<float>> A, vector<vector<float>> B)
{
	int row = A.size();
	int col = B[0].size();
	vector<vector<float>> temp;

	temp.resize(row);
	for (int i = 0; i < row; i++)
	{
		temp[i].resize(col);
		for (int j = 0; j < col; j++)
			temp[i][j] = 0;
	}

	for (unsigned i = 0; i < row; i++)
		for (unsigned j = 0; j < col; j++)
			for (unsigned k = 0; k < row; k++)
				temp[i][j] += A[i][k] * B[k][j];
	return temp;
}

void Frobenius(vector<vector<float>> &A, vector<vector<vector<float>>> &B)
{
	int n = A.size();
	B.resize(n - 1);

	Matrix Invers = Matrix(n, n);

	for (int i = 0; i < n - 1; i++)	
	{
		B[i] = E(n);
		for (int j = 0; j < n; j++)
		{
			B[i][n - (i + 2)][j] = (j != n - (i + 2))

				? -A[n - (i + 1)][j] / A[n - (i + 1)][n - (i + 2)]

				: 1 / A[n - (i + 1)][n - (i + 2)];
		}
		Invers.Set(B[i]);
		Invers.Inverse();
		MUL(A, B[i], Invers, n);
	}
}

void Out(vector<vector<float>> array)
{
	for (int i = 0; i < array.size(); i++)
	{
		for (int j = 0; j < array[0].size(); j++)
		{
			cout << array[i][j] << "  ";
		}
		cout << endl;
	}
}

int main()
{
	setlocale(0, "");

	vector<vector<float>> A = 
	{
		{4, -2, 2, 1, 3},
		{-2, 1, 0, -4, -3},
		{2, -3, -2, 0, -2},
		{4, 5, 0, 5, 5},
		{-5, -2, -2, 5, 3}
	};
	vector<vector<vector<float>>> B;
	Frobenius(A, B);

	for (int i = 0; i < B.size(); i++)
	{
		cout << "B" << i + 1 << endl;
		Out(B[i]);
	}
	cout << "A" << endl;
	Out(A);


	//Cмотрите фото в этом же каталоге
	vector<vector<float>> Y =
	{
		{29.9611},
		{-12.8061},
		{5.47368},
		{-2.339559},
		{1}
	};
	cout << "X" << endl;
	Out(MUL(MUL(MUL(MUL(B[0], B[1]), B[2]), B[3]), Y));


	system("pause");
	return 0;
}