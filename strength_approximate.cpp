//MIT License
//
//Copyright(c) 2017 Ivan Blagopoluchnyy
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files(the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions :
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

#include "strength_approximate.h"

/*-----WINAPI CONFIG START-----*/
void SetColor(int text, int background)
{
  HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
  SetConsoleTextAttribute(hStdOut, (WORD)((background << 4) | text));
}

void SetColorBgTex(int Bg = 0, int Tex = 15)
{
  HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
  Bg = Bg << 4;
  Bg = Bg | Tex;
  SetConsoleTextAttribute(hConsole, Bg);
}
/*-----WINAPI CONFIG END-----*/

/*-------------------MATRIX CALCULATIONS--------------------------*/
void evaluate_matrix(double& radius, double& x_val, double& phi_val,
	double& h_val, double& E_mod_val, double& nu_val, double& delta_phi_val, double& delta_x_val)
{
	//initialization non - zero values of matrix A1, A2, A3, A4, A5
  //Non - zero matrix elements for the canonical system of differential equations 
  //that characterize the state of the cross section along the normal of the cylindrical shell
	double A1_a12 = -nu_val / radius, 
		   A1_a21 = -1.0 / radius, 
		   A1_a42 = -nu_val / (radius*radius), 
		   A1_a56 = -1.0 / radius,
		   A1_a63 = -(E_mod_val*h_val) / (radius*radius), 
		   A1_a65 = -nu_val / radius, 
		   A1_a68 = -nu_val / (radius*radius),
		   A1_a72 = (E_mod_val*h_val) / (radius*radius); //A1_ij != 0

	double A2_a43 = nu_val / (radius*radius), 
		   A2_a51 = -(E_mod_val*h_val*h_val*h_val) / (6.0*(1.0 + nu_val)*radius*radius*radius*radius),
		   A2_a54 = (E_mod_val*h_val*h_val*h_val) / (6.0*(1.0 + nu_val)*radius*radius*radius), 
		   A2_a62 = -(E_mod_val*h_val) / (radius*radius),
		   A2_a78 = -nu_val / (radius*radius), 
		   A2_a81 = (E_mod_val*h_val*h_val*h_val) / (6.0*(1.0 + nu_val)*radius*radius*radius),
		   A2_a84 = (E_mod_val*h_val*h_val*h_val) / (6.0*(1.0 + nu_val)*radius*radius); //A2_ij != 0

	double A3_a63 = (E_mod_val*h_val*h_val*h_val) / (12.0*radius*radius*radius*radius), 
		   A3_a72 = -(E_mod_val*h_val*h_val*h_val) / (12.0*radius*radius*radius*radius); //A3_ij != 0

	double A4_a73 = (E_mod_val*h_val*h_val*h_val) / (12.0*radius*radius*radius*radius); //A4_ij != 0 

	double A5_a13 = -nu_val / radius, 
		   A5_a15 = (1.0 + nu_val*nu_val) / (E_mod_val*h_val), 
		   A5_a26 = (2.0 * (1.0 + nu_val)) / E_mod_val*h_val,
		   A5_a34 = -1.0, 
		   A5_a48 = (12.0 * (1.0 - nu_val*nu_val)) / E_mod_val*h_val*h_val*h_val, 
		   A5_a73 = (E_mod_val*h_val) / (radius*radius),
		   A5_a75 = nu_val / radius, 
		   A5_a87 = 1.0; //A5_ij != 0

	/*------------CALCULATION TENSOR COMPONENTS A1, A2, A3, A4, A5-----------*/
	MatrixXd A1(DIM_A, DIM_A);
	A1.setZero(); //set all components to zero
	A1(0, 1) = A1_a12, A1(1, 0) = A1_a21, A1(3, 1) = A1_a42, A1(4, 5) = A1_a56;
	A1(5, 2) = A1_a63, A1(5, 4) = A1_a65, A1(5, 7) = A1_a68, A1(6, 1) = A1_a72;
	std::cout << "\nMatrix А1:\n" << A1 << std::endl;

	MatrixXd A2(DIM_A, DIM_A);
	A2.setZero(); //set all components to zero
	A2(3, 2) = A2_a43, A2(4, 0) = A2_a51, A2(4, 3) = A2_a54, A2(5, 1) = A2_a62;
	A2(6, 7) = A2_a78, A2(7, 0) = A1_a65, A2(7, 3) = A2_a84;
	std::cout << "\nMatrix А2:\n" << A2 << std::endl;

	MatrixXd A3(DIM_A, DIM_A);
	A3.setZero(); //set all components to zero
	A3(5, 2) = A3_a63;
	A3(6, 1) = A3_a72;
	std::cout << "\nMatrix А3:\n" << A3 << std::endl;

	MatrixXd A4(DIM_A, DIM_A);
	A4.setZero(); //set all components to zero
	A4(6, 2) = A4_a73;
	std::cout << "\nMatrix А4:\n" << A4 << std::endl;

	MatrixXd A5(DIM_A, DIM_A);
	A5.setZero(); //set all components to zero
	A5(0, 2) = A5_a13, A5(0, 4) = A5_a15, A5(1, 5) = A5_a26, A5(2, 3) = A5_a34;
	A5(3, 7) = A5_a48, A5(6, 2) = A5_a73, A5(6, 4) = A5_a75, A5(7, 6) = A5_a87;
	std::cout << "\nMatrix А5:\n" << A5 << std::endl;

	/*-------CALCULATION OF FINITE DIFFERENCES A_i+2, A_i+1, A_i, A_i-1, A_i-2---------*/

	MatrixXd A_ip2(DIM_A, DIM_A), A_ip1(DIM_A, DIM_A), A_i(DIM_A, DIM_A), A_im1(DIM_A, DIM_A), A_im2(DIM_A, DIM_A);
	A_ip2.setZero(), A_ip1.setZero(), A_i.setZero(), A_im1.setZero(), A_im2.setZero();

	double a1_block = (1.0 / (1.0*pow(delta_phi_val, 3.0))), a2_block = (1.0 / (pow(delta_phi_val, 4.0)));
	double b1_block = (1.0 / (2.0*delta_phi_val)), b2_block = (1.0 / (pow(delta_phi_val, 2.0))), b3_block = (1.0 / (pow(delta_phi_val, 3.0)));
	double b4_block = (4.0 / (pow(delta_phi_val, 4.0))), c1_block = (2.0 / (pow(delta_phi_val, 2.0))), c2_block = (6.0 / (pow(delta_phi_val, 4.0)));
	double d1_block = (1.0 / (pow(delta_phi_val, 3.0))), e1_block = (1.0 / (2.0*pow(delta_phi_val, 3.0)));

	for (int i = 0; i < DIM_A; i++)
	{
		for (int j = 0; j < DIM_A; j++)
		{
			A_ip2(i, j) = a1_block*A3(i, j) + a2_block*A4(i, j);
			A_ip1(i, j) = b1_block*A1(i, j) + b2_block*A2(i, j) - b3_block*A3(i, j) - b4_block*A4(i, j);
			A_i(i, j) = (-c1_block)*A2(i, j) + c2_block*A4(i, j) + A5(i, j);
			A_im1(i, j) = (-b1_block)*A1(i, j) + b2_block*A2(i, j) - d1_block*A3(i, j) - b4_block*A4(i, j);
			A_im2(i, j) = (-e1_block)*A3(i, j) + a2_block*A4(i, j);
		}
	}

	std::cout << "\nMatrix А i+2:\n" << A_ip2 << std::endl;
	std::cout << "\nMatrix А i+1:\n" << A_ip1 << std::endl;
	std::cout << "\nMatrix А i:\n" << A_i << std::endl;
	std::cout << "\nMatrix А i-1:\n" << A_im1 << std::endl;
	std::cout << "\nMatrix А i-2:\n" << A_ip2 << std::endl;

	int num_n = phi_val / delta_phi_val; //amody > a > y, n is order for SDE matrix

	std::cout << "\nOrder for SDE:\n\nn = " << num_n << std::endl;

	///*---------------------------------CREATION OF THE SDE MATRIX----------------------------------------*/

	MatrixXd SDE(num_n + 1, num_n + 1);
	SDE.setZero();

  //TODO[1]
	/*------TEST BEGIN-----*/

	//Matrix SDE = Matrix(A_im2, A_im1, A_i, A_ip1, A_ip2, num_n);

	std::cout << "\nSDE matrix:\n" << SDE << std::endl;


	//Eigen::Tensor<double, 3> epsilon(3, 3, 3);
	//epsilon.setZero();
	//epsilon(0, 1, 2) = 1;
	//epsilon(1, 2, 0) = 1;
	//epsilon(2, 0, 1) = 1;
	//epsilon(1, 0, 2) = -1;
	//epsilon(2, 1, 0) = -1;
	//epsilon(0, 2, 1) = -1;
	//Eigen::Tensor<double, 4> grassmannIdentity(3, 3, 3, 3);
	//grassmannIdentity.setZero();

	//// this is not the most efficient way to write such a product,
	//// but is the only way possible with the current feature set
	//for (int i = 0; i < 3; i++) {
	//	for (int j = 0; j < 3; j++) {
	//		for (int k = 0; k < 3; k++) {
	//			for (int l = 0; l < 3; l++) {
	//				for (int m = 0; m < 3; m++) {
	//					grassmannIdentity(i, j, l, m) += epsilon(i, j, k) * epsilon(k, l, m);
	//				}
	//			}
	//		}
	//	}
	//}

	//// verify
	//for (int i = 0; i < 3; i++) {
	//	for (int j = 0; j < 3; j++) {
	//		for (int l = 0; l < 3; l++) {
	//			for (int m = 0; m < 3; m++) {
	//				assert(grassmannIdentity(i, j, l, m) == (int(i == l) * int(j == m) - int(i == m) * int(j == l)));
	//			}
	//		}
	//	}
	//}

	//std::cout << grassmannIdentity << std::endl;

	//// dimensionalities
	//assert(epsilon.dimension(0) == 3);
	//assert(epsilon.dimension(1) == 3);
	//assert(epsilon.dimension(2) == 3);
	//auto dims = epsilon.dimensions();
	//assert(dims[0] == 3);
	//assert(dims[1] == 3);
	//assert(dims[2] == 3);

	//int num_m = x_val / delta_x_val; //amody > a > y, m is order for Cauchy-Krylov matrix

	//std::cout << "\nOrder for Cauchy-Krylov matrix:\n\nm = " << num_m << std::endl;

	///*---------------------------------CALCULATION OF THE MATRIX EXPONENTIAL---------------------------------*/

	////dimension of the output matrix must be equal SDE matrix dimension
	//Matrix unit_matrix = Matrix(SDE.getRowsCount(), SDE.getColumnsCount());

	////initialization of the unit matrix
	//unit_matrix.getUnitMatrix();

	//std::cout << "\nUnit matrix:\n" << unit_matrix << std::endl;

	//double value_x = 0; //delta_x^m
	//
	////calculate series sum where m <= n
	//for (int t = 1; t <= num_m; t++)
	//{
	//	SDE.power_matrix(t); //pow SDE matrix due to t: SDE^t
	//	value_x = value_x + pow(delta_x_val, t); //calculate delta_x^m 
	//	value_x = value_x / t; //division current value by t: delta_x^m/t
	//	SDE.div(t*value_x); //division 
	//	unit_matrix.mSum(SDE); //result of summing with unit matrix
	//}

	//std::cout << "\nCauchy-Krylov matrix:\n" << unit_matrix << std::endl;

	///*---------------------CALCULATION P & Q MATRICES BY LU - ALGORITHM---------------------------*/

	////matrices U and L must have dimension that equals to Cauchy-Krylov matrix dimension
	//Matrix U_matrix = Matrix(unit_matrix.getRowsCount(), unit_matrix.getColumnsCount());
	//Matrix L_matrix = Matrix(unit_matrix.getRowsCount(), unit_matrix.getColumnsCount());

	////make U unit matrix
	//U_matrix.getUnitMatrix();

	//int i = 0, j = 0, k = 0, p = 0;
	//double sum = 0.0;

	///*---------------------LU - ALGORITHM---------------------------*/
	//for (k = 1; k <= unit_matrix.getRowsCount(); k++)
	//{

	//	for (i = k; i <= unit_matrix.getColumnsCount(); i++)
	//	{
	//		sum = 0.0;

	//		for (p = 1; p <= k - 1; p++)
	//			sum += L_matrix.getElement(i, p) * U_matrix.getElement(p, k);
	//		L_matrix.setElement(i, k, unit_matrix.getElement(i, k) - sum);

	//	}

	//	for (j = k + 1; j <= unit_matrix.getColumnsCount(); j++)
	//	{
	//		sum = 0.0;

	//		for (p = 1; p <= k - 1; p++)
	//			sum += L_matrix.getElement(k, p)*U_matrix.getElement(p, j);
	//		U_matrix.setElement(k, j, (unit_matrix.getElement(k, j) - sum) / L_matrix.getElement(k, k));
	//	}

	//}

	//std::cout << "\nP - matrix is (U - matrix):\n" << U_matrix << std::endl;
	//std::cout << "\nQ - matrix is (L - matrix):\n" << L_matrix << std::endl;

	////create matrix for verification
	//Matrix check_matrix = Matrix(unit_matrix.getRowsCount(), unit_matrix.getColumnsCount());

	//double check_value = 0.0;

	////check it by the matrix product
	//for (int i = 0; i < unit_matrix.getRowsCount(); i++)
	//{
	//	for (int j = 0; j < unit_matrix.getRowsCount(); j++)
	//	{
	//		check_value = 0.0;

	//		for (int l = 0; l < unit_matrix.getRowsCount(); l++)
	//			check_value += (L_matrix.getElement(i, l) * U_matrix.getElement(l, j));
	//		check_matrix.setElement(i, j, check_value);

	//	}
	//}

	//std::cout << "\nA = P*Q:\n" << check_matrix << std::endl;

	///*--------------SUMMING WITH UNIT MATRIX-----------------*/
	////for the ideal case matrix of the external forces it is equal unit matrix

	//Matrix ExteriorForcesMatrix = Matrix(check_matrix.getRowsCount(), check_matrix.getRowsCount());

	//ExteriorForcesMatrix.getUnitMatrix();

	//check_matrix.mSum(ExteriorForcesMatrix);

	//std::cout << "\nResult of summing with matrix of external forces:\n" << check_matrix << std::endl;

	///*---------------------CALCULATE CAUCHY-KRYLOV MATRIX FOR RESULT MATRIX--------------------------*/
	//Matrix unit_matrix_kk = Matrix(check_matrix.getRowsCount(), check_matrix.getColumnsCount());

	//unit_matrix_kk.getUnitMatrix();

	//double value_variable = 0;

	////it is similar where is calculation of the matrix exponential located
	//for (int t = 1; t <= num_m; t++) 
	//{
	//	check_matrix.power_matrix(t);
	//	value_variable = value_variable + pow(delta_x_val, t);
	//	value_variable = value_variable / t;
	//	check_matrix.div(t*value_variable);
	//	unit_matrix_kk.mSum(check_matrix);
	//}

	//std::cout << "\nCauchy-Krylov matrix:\n" << unit_matrix_kk << std::endl;

	////save result matrix elements in data.txt
	//save_data(unit_matrix_kk);
  /*------TEST END-----*/
}

char GetInput()
{
  char choice;
  std::cin >> choice;
  return choice;
}

void DisplayMainMenu()
{
  std::cout << "\n1 - enter values" << std::endl;
  std::cout << "2 - use built in balues" << std::endl;
  std::cout << "0 - exit\n";
}

/*--------------------------------DATA INPUT-------------------------------------*/
void enter_data()//use and initialize input values
{
  std::cout << "\nInput of initial data\n" << std::endl;

  double R = 0.0, x = 0.0, phi = 0.0, h = 0.0, E_mod = 0.0, nu = 0.0,
    phi_buf = 0.0, delta_phi = 0.0, delta_pbuf = 0.0, delta_x = 0.0;
  bool flag = true;
  //double step = 0.0; //phi step - you can use or ignore it

  std::cin.ignore(); //ignore 'buffer'

  do
  {
    DisplayMainMenu();
    switch (GetInput())
    {
    case '1':
    {

      std::cout << "\nEnter values:" << std::endl;
      std::cout << "\nR = "; std::cin >> R;
      std::cout << "\nx = "; std::cin >> x;
      std::cout << "\nphi (°C) = "; std::cin >> phi;
      std::cout << "\nh = "; std::cin >> h;
      std::cout << "\nE = "; std::cin >> E_mod;
      std::cout << "\nnu = "; std::cin >> nu;
      //std::cout << "\nStep = "; std::cin >> step;
      std::cout << "\ndelta(phi) (°C) = "; std::cin >> delta_phi;
      std::cout << "\ndelta(x) = "; std::cin >> delta_x; //in meters

      phi_buf = phi * (M_PI / 180.0); //phi in radians
      delta_pbuf = delta_phi * (M_PI / 180.0);

      //calculate each strength tensor component
      evaluate_matrix(R, x, phi_buf, h, E_mod, nu, delta_pbuf, delta_x);

      break;
    }
    case '2':
    {

      std::cout << "\nDefault values:" << std::endl;
      //in SI (m, Pa), degree (phi) in radians

      //is given
      R = 2.0; x = 5.0; phi = 60.0; h = 0.1; nu = 0.34; E_mod = 6180;

      std::cout << "\nR = " << R << std::endl;
      std::cout << "\nx = " << x << std::endl;
      std::cout << "\nphi = " << phi << " °C" << std::endl; //in °C
      std::cout << "\nh = " << h << std::endl;
      std::cout << "\nE = " << E_mod << std::endl;
      std::cout << "\nnu = " << nu << std::endl;
      //std::cout << "\nStep = " << step << std::endl;
      std::cout << "\nenter delta(phi) (°C) = "; std::cin >> delta_phi; //enter in °C
      std::cout << "\ndelta(x) = "; std::cin >> delta_x; //enter in m

      phi_buf = phi * (M_PI / 180.0); //phi in radians
      delta_pbuf = delta_phi * (M_PI / 180.0);

      //calculate each strength tensor component
      evaluate_matrix(R, x, phi_buf, h, E_mod, nu, delta_pbuf, delta_x);

      break;
    }
    case '0':
      flag = false;
      break;

    default:
    {
      SetColor(12, 0);
      std::cout << "\nIncorrect choise!" << std::endl;
      SetColor(15, 0);
      break;
    }
    }
  } while (flag != false);
}

//TODO[2]
//void save_data(MatrixXd unit_matrix_kk)
//{
//	std::ofstream save_data("data.txt");
//	//second parameter: ios_base::app - to append data to file
//
//	if (!save_data.is_open())
//		std::cout << "\nCannot create data.txt" << std::endl;
//	else
//	{
//		save_data << "Result matrix M is:" << std::endl;
//
//		for (int i = 0; i < unit_matrix_kk.rows(); i++)
//		{
//			for (int j = 0; j < unit_matrix_kk.cols(); j++)
//			{
//				save_data << "\nM[" << i + 1 << "][" << j + 1 << "] = " << unit_matrix_kk(i, j);
//			}
//		}
//
//		std::cout << "\nValues of M matrix are stored in data.txt" << std::endl;
//	}
//}

/*-------------------------------MAIN FUNCTION--------------------------------------*/
int main(void)
{

	setlocale(LC_ALL, "");
	//if you have an a 1251/866 encoding in yours console output
	//you can use it winapi cheats (check it with chcp cmd):
	//SetConsoleCP(1251);
	//SetConsoleOutputCP(1251);

	std::streamsize ss = std::cout.precision();

	std::cout << "\nEnter the number of decimal places:\n" << std::endl;

	int k_enter = 0;

	std::cin >> k_enter;

	std::cout.precision(k_enter + 1);

	std::cin.clear();

	enter_data();

	std::cin.get();
}