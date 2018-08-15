// Matrix.h

#ifndef MATRIX_H
#define MATRIX_H

const int DIM_A = 8;

class Matrix
{

private:

  int rows;
  int cols;
  std::vector<std::vector<long double>> matrix;

public:


  /*Matrix Constructor*/
  Matrix(int rows, int cols);


  /*Diagonal tensor constructor*/
  Matrix(Matrix& A_im2, Matrix& A_im1, Matrix& A_i, Matrix& A_ip1, Matrix& A_ip2, int& dimension);


  ~Matrix() {};


  int getRowsCount();


  int getColumnsCount();


  /*-------------------GET CURRENT MATRIX---------------------*/
  std::vector < std::vector <long double> > getMatrix();


  /*-------------------CREATING A UNIT MATRIX---------------------*/
  std::vector < std::vector <long double> > getUnitMatrix();


  /*-------------------DIVISION OF THE MATRIX TO NUMBER-----------------------*/
  std::vector < std::vector <long double> > div(double x);


  /*------------------------SUMMING OF MATRIXES-----------------------*/
  Matrix mSum(Matrix obj);


  /*-------------------RAISE MATRIX TO DEGREE-------------------*/
  Matrix operator ^ (int order);


  /*-------------------GET MATRIX ELEMENT BY COLUMN AND ROW INDEX-------------------*/
  double getElement(int rowIndex, int colIndex);


  /*-------------------ADD ELEMENT INTO MATRIX BY COLUMN AND ROW INDEX-------------------*/
  void setElement(int rowIndex, int colIndex, double val);


  /*-------------------PRINT OUT MATRIX-------------------*/
  void printMatrix();

};

#endif