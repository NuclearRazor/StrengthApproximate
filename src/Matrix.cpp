// Matrix.cpp

#ifndef MATRIX_CPP
#define MATRIX_CPP

#include "../include/Matrix.h"


Matrix::Matrix(int rows, int cols)
{
  this->rows = rows;
  this->cols = cols;

  for (int i = 0; i<rows; i++)
  {
    std::vector<long double> curRow;
    for (int j = 0; j<cols; j++)
      curRow.push_back(0.0);
    matrix.push_back(curRow);
  }
}


Matrix::Matrix(Matrix& A_im2, Matrix& A_im1, Matrix& A_i, Matrix& A_ip1, Matrix& A_ip2, int& dimension)
{
  this->rows = dimension;
  this->cols = dimension;

  //how many matrices will be filled in current tensor, The implementation is correct, but it is necessary to take into the formulation of the problem
  int fill = dimension / DIM_A;

  //counter for filled matrices
  int push = 0;

  Matrix& cur = A_i;

  for (int i = 0; i < this->rows; i++)
  {
    if (push == 0 || push % 5 == 0)
      cur = A_im2;
    else if (push == 1 || push % 6 == 0)
      cur = A_im1;
    else if (push == 2 || push % 7 == 0)
      cur = A_i;
    else if (push == 3 || push % 8 == 0)
      cur = A_ip1;
    else if (push == 4 || push % 9 == 0)
      cur = A_ip2;

    //vector <-> string that we fill in the matrix
    std::vector <long double> curRow;

    for (int j = 0; j < this->cols; j++)
    {
      //fill with zeros "left" (these are zeros under the main diagonal of the tensor)
      if (j < push*DIM_A)
        curRow.push_back(0.0);

      //fill tensor by values
      curRow.push_back(cur.getElement(i - (push*DIM_A), j - (push*DIM_A)));
    }

    matrix.push_back(curRow);

    if (i != 0 && i % 7 == 0)
      push++;
  }
}


int Matrix::getRowsCount()
{
  return this->rows;
}


int Matrix::getColumnsCount()
{
  return this->cols;
}


std::vector<std::vector<long double>> Matrix::getMatrix()
{
  return matrix;
}


std::vector < std::vector <long double> > Matrix::getUnitMatrix()
{
  Matrix local_unit_matrix(this->rows, this->cols);

  for (int i = 0; i < this->rows; i++)
  {
    for (int j = 0; j < this->cols; j++)
    {
      //for the main diagonal // || j == DIM - i - 1) // for the secondary diagonal
      if (i == j)
        local_unit_matrix.setElement(i, j, 1.0);
    }
  }

  matrix = local_unit_matrix.getMatrix();

  return matrix;
}


std::vector < std::vector <long double> > Matrix::div(double x)
{
  Matrix local_unit_matrix(this->rows, this->cols);

  //you can use assert

  //use simple logic to catch division by zero
  if (x == 0)
  {
    std::cout << "\nDivision by zero!\n" << std::endl;
  }
  //else divide each matrix element by scalar value
  else
  {
    for (int i = 0; i < this->rows; i++)
      for (int j = 0; j < this->cols; j++)
        local_unit_matrix.setElement(i, j, matrix[i][j] / x);
  }

  matrix = local_unit_matrix.getMatrix();

  return matrix;
}


Matrix Matrix::mSum(Matrix obj)
{
  Matrix local_matrix(this->rows, this->cols);

  std::vector<std::vector<long double>> cur = obj.getMatrix();

  double lastVal = 0.0;

  for (int i = 0; i < this->rows; i++)
  {
    for (int j = 0; j < this->cols; j++)
    {
      lastVal = 0.0;
      lastVal = cur[i][j] + matrix[i][j];
      local_matrix.setElement(i, j, lastVal);
    }
  }

  //update matrix
  matrix = local_matrix.getMatrix();

  return *this;
}


void Matrix::power_matrix(int order)
{
  Matrix local_matrix(this->rows, this->cols);

  std::vector<std::vector<long double>> cur = matrix;

  double lastVal = 0.0;

  for (int k = 0; k < order - 1; k++)
  {
    if (k > 0)
      cur = local_matrix.getMatrix();

    for (int i = 0; i < this->rows; i++)
    {
      for (int j = 0; j < this->cols; j++)
      {
        lastVal = 0.0;

        // min - check for the smallest value, the matrix must be square
        for (int l = 0; l < std::min(this->rows, this->cols); l++)
        {
          lastVal += cur[i][l] * matrix[l][j];
        }
        local_matrix.setElement(i, j, lastVal);
      }
    }
  }

  if (order > 1)
  {
    matrix = local_matrix.getMatrix();
  }
}


double Matrix::getElement(int rowIndex, int colIndex)
{
  if ((rowIndex < 0 || rowIndex >= this->rows) || (colIndex < 0 || colIndex >= this->cols))
  {
    return 0.0;
  }
  else
  {
    return this->matrix.at(rowIndex).at(colIndex);
  }
}


void Matrix::setElement(int rowIndex, int colIndex, double val)
{
  if ((rowIndex >= 0 && rowIndex < this->rows) && (colIndex >= 0 && colIndex < this->cols))
  {
    this->matrix.at(rowIndex).at(colIndex) = val;
  }
}


void Matrix::print_matrix()
{
  for (std::vector<std::vector<int>>::size_type i = 0; i < matrix.size(); i++)
  {
    for (std::vector<std::vector<int>>::size_type j = 0; j < matrix[i].size(); j++) { std::cout << "\t" << matrix[i][j]; }
    std::cout << "\n";
  }
}


#endif