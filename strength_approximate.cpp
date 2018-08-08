//MIT License
//
//Copyright(c) 2017 - 2018 Ivan Blagopoluchnyy
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

#include "include/strength_approximate.h"


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


/*---------------UI LOGIC START---------------*/
char GetInput()
{
  char choice;
  std::cin >> choice;
  return choice;
}


void DisplayMainMenu()
{
  std::cout << "\n1 - enter values" << std::endl;
  std::cout << "2 - use built in values" << std::endl;
  std::cout << "0 - exit\n";
}

/*-----PYTHON SCRIPT CALLABLE START-----*/
void CallPyPlotter()
{
  char filename[] = "momentum.py";
  FILE* fp;

  // initialuze python interpreter
  Py_Initialize();

  // pass file name in read mode
  fp = _Py_fopen(filename, "r");

  // run python script
  PyRun_SimpleFile(fp, filename);

  // finalize callables
  Py_Finalize();
}
/*-----PYTHON SCRIPT CALLABLE END-----*/
/*---------------UI LOGIC START---------------*/



/*-------------------MATRIX CALCULATIONS START--------------------------*/
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
    A5_a15 = (1.0 + nu_val * nu_val) / (E_mod_val*h_val),
    A5_a26 = (2.0 * (1.0 + nu_val)) / E_mod_val * h_val,
    A5_a34 = -1.0,
    A5_a48 = (12.0 * (1.0 - nu_val * nu_val)) / E_mod_val * h_val*h_val*h_val,
    A5_a73 = (E_mod_val*h_val) / (radius*radius),
    A5_a75 = nu_val / radius,
    A5_a87 = 1.0; //A5_ij != 0

                  /*------------CALCULATION TENSOR COMPONENTS A1, A2, A3, A4, A5-----------*/

  Matrix A1 = Matrix(DIM_A, DIM_A);
  A1.setElement(0, 1, A1_a12);
  A1.setElement(1, 0, A1_a21);
  A1.setElement(3, 1, A1_a42);
  A1.setElement(4, 5, A1_a56);
  A1.setElement(5, 2, A1_a63);
  A1.setElement(5, 4, A1_a65);
  A1.setElement(5, 7, A1_a68);
  A1.setElement(6, 1, A1_a72);
  std::cout << "\nMatrix A1:\n" << std::endl;
  A1.print_matrix();

  Matrix A2 = Matrix(DIM_A, DIM_A);
  A2.setElement(3, 2, A2_a43);
  A2.setElement(4, 0, A2_a51);
  A2.setElement(4, 3, A2_a54);
  A2.setElement(5, 1, A2_a62);
  A2.setElement(6, 7, A2_a78);
  A2.setElement(7, 0, A1_a65);
  A2.setElement(7, 3, A2_a84);
  std::cout << "\nMatrix A2:\n" << std::endl;
  A2.print_matrix();

  Matrix A3 = Matrix(DIM_A, DIM_A);
  A3.setElement(5, 2, A3_a63);
  A3.setElement(6, 1, A3_a72);
  std::cout << "\nMatrix A3:\n" << std::endl;
  A3.print_matrix();

  Matrix A4 = Matrix(DIM_A, DIM_A);
  A4.setElement(6, 2, A4_a73);
  std::cout << "\nMatrix A4:\n" << std::endl;
  A4.print_matrix();

  Matrix A5 = Matrix(DIM_A, DIM_A);
  A5.setElement(0, 2, A5_a13);
  A5.setElement(0, 4, A5_a15);
  A5.setElement(1, 5, A5_a26);
  A5.setElement(2, 3, A5_a34);
  A5.setElement(3, 7, A5_a48);
  A5.setElement(6, 2, A5_a73);
  A5.setElement(6, 4, A5_a75);
  A5.setElement(7, 6, A5_a87);
  std::cout << "\nMatrix A5:\n" << std::endl;
  A5.print_matrix();

  /*-------CALCULATION OF FINITE DIFFERENCES A_i+2, A_i+1, A_i, A_i-1, A_i-2---------*/

  //fill matrices with zeros
  Matrix A_ip2 = Matrix(DIM_A, DIM_A);
  Matrix A_ip1 = Matrix(DIM_A, DIM_A);
  Matrix A_i = Matrix(DIM_A, DIM_A);
  Matrix A_im1 = Matrix(DIM_A, DIM_A);
  Matrix A_im2 = Matrix(DIM_A, DIM_A);

  for (int i = 0; i < DIM_A; i++)
  {
    for (int j = 0; j < DIM_A; j++)
    {
      A_ip2.setElement(i, j, (1.0 / (1.0*pow(delta_phi_val, 3.0)))*A3.getElement(i, j) + (1.0 / (pow(delta_phi_val, 4.0)))*A4.getElement(i, j));
      A_ip1.setElement(i, j, (1.0 / (2.0*delta_phi_val))*A1.getElement(i, j) + (1.0 / (pow(delta_phi_val, 2.0)))*A2.getElement(i, j)
        - (1.0 / (pow(delta_phi_val, 3.0)))*A3.getElement(i, j) - (4 / (pow(delta_phi_val, 4.0)))*A4.getElement(i, j));
      A_i.setElement(i, j, -(2.0 / (pow(delta_phi_val, 2.0)))*A2.getElement(i, j)
        + (6.0 / (pow(delta_phi_val, 4.0)))*A4.getElement(i, j) + A5.getElement(i, j));
      A_im1.setElement(i, j, -(1.0 / (2.0*delta_phi_val))*A1.getElement(i, j) + (1.0 / pow(delta_phi_val, 2.0))*A2.getElement(i, j)
        - (1.0 / (pow(delta_phi_val, 3.0)))*A3.getElement(i, j) - (4.0 / (pow(delta_phi_val, 4.0)))*A4.getElement(i, j));
      A_im2.setElement(i, j, -(1.0 / (2.0*pow(delta_phi_val, 3.0)))*A3.getElement(i, j) + (1.0 / pow(delta_phi_val, 4.0))*A4.getElement(i, j));
    }
  }

  std::cout << "\n А i+2:\n" << std::endl;
  A_ip2.print_matrix();

  std::cout << "\nА i+1:\n" << std::endl;
  A_ip1.print_matrix();

  std::cout << "\nА i:\n" << std::endl;
  A_i.print_matrix();

  std::cout << "\nА i-1:\n" << std::endl;
  A_im1.print_matrix();

  std::cout << "\nА i-2:\n" << std::endl;
  A_im2.print_matrix();

  int num_n;
  num_n = static_cast<int> (phi_val / delta_phi_val); //amody > a > y, SDE order

  std::cout << "\nSDE order:\n\nn = " << num_n << std::endl;

  /*---------------------------------SDE EVALUATION----------------------------------------*/

  //SDE fills by finite differences
  Matrix SDE = Matrix(A_im2, A_im1, A_i, A_ip1, A_ip2, num_n);

  std::cout << "\nSDE tensor:\n" << std::endl;

  SDE.print_matrix();

  int num_m;
  num_m = static_cast<int> (x_val / delta_x_val); //amody> a> y, the order m for the Cauchy-Krylov matrix

  std::cout << "\nThe exponent series for calculating the Cauchy-Krylov matrix:\n\nm = " << num_m << std::endl;

  /*---------------------------------MATRIX EXPONENT EVALUATION---------------------------------*/

  //Create identity matrix, where stores final calculations
  Matrix unit_matrix = Matrix(SDE.getRowsCount(), SDE.getColumnsCount());

  unit_matrix.getUnitMatrix();

  std::cout << "\nThe identity matrix:\n" << std::endl;

  unit_matrix.print_matrix();

  double value_x = 0;

  for (int t = 1; t <= num_m; t++) //m <= n
  {
    SDE.power_matrix(t);

    value_x = value_x + pow(delta_x_val, t);

    value_x = value_x / t; //devide current value by delta_x^m/t

    SDE.div(t*value_x);

    unit_matrix.mSum(SDE);
  }

  std::cout << "\nCauchy-Krylov matrix for unit forces matrix is:\n" << std::endl;

  unit_matrix.print_matrix();

  /*---------------------COMPUTATION OF MATRICES P AND Q BY MEANS OF LU - DECOMPOSITIONS---------------------------*/

  Matrix U_matrix = Matrix(unit_matrix.getRowsCount(), unit_matrix.getColumnsCount());

  Matrix L_matrix = Matrix(unit_matrix.getRowsCount(), unit_matrix.getColumnsCount());

  U_matrix.getUnitMatrix();

  int i = 0, j = 0, k = 0, p = 0;
  double sum = 0.0;

  /*---------------------LU START---------------------------*/
  for (k = 1; k <= unit_matrix.getRowsCount(); k++)
  {

    for (i = k; i <= unit_matrix.getColumnsCount(); i++)
    {
      sum = 0.0;

      for (p = 1; p <= k - 1; p++)
        sum += L_matrix.getElement(i, p) * U_matrix.getElement(p, k);

      L_matrix.setElement(i, k, unit_matrix.getElement(i, k) - sum);

    }

    for (j = k + 1; j <= unit_matrix.getColumnsCount(); j++)
    {
      sum = 0.0;

      for (p = 1; p <= k - 1; p++)
        sum += L_matrix.getElement(k, p)*U_matrix.getElement(p, j);

      U_matrix.setElement(k, j, (unit_matrix.getElement(k, j) - sum) / L_matrix.getElement(k, k));
    }

  }
  /*---------------------LU END---------------------------*/

  std::cout << "\nP - matrix (U):\n" << std::endl;

  U_matrix.print_matrix();

  std::cout << "\nQ - matrix (L):\n" << std::endl;

  L_matrix.print_matrix();

  //check calculations
  Matrix check_matrix = Matrix(unit_matrix.getRowsCount(), unit_matrix.getColumnsCount());

  double check_value = 0.0;

  for (int i = 0; i < unit_matrix.getRowsCount(); i++)
  {
    for (int j = 0; j < unit_matrix.getRowsCount(); j++)
    {
      check_value = 0.0;

      for (int l = 0; l < unit_matrix.getRowsCount(); l++)
        check_value += (L_matrix.getElement(i, l) * U_matrix.getElement(l, j));
      check_matrix.setElement(i, j, check_value);

    }
  }

  std::cout << "\nA = P*Q:\n" << std::endl;

  check_matrix.print_matrix();

  /*--------------SUM CAUCHY-KRYLOV MATRIX WITH MATRIX OF EXTERIOR FORCES-----------------*/
  //it can also be initialized, for an ideal case we take it as a single

  Matrix ExteriorForcesMatrix = Matrix(check_matrix.getRowsCount(), check_matrix.getRowsCount());

  ExteriorForcesMatrix.getUnitMatrix();

  check_matrix.mSum(ExteriorForcesMatrix);

  std::cout << "\n:\n" << std::endl;

  check_matrix.print_matrix();

  /*---------------------FINAL EVALUATION (DECOMPOSE BY EXP SERIES) FOR (CAUCHY-KRYLOV + EXTERIOR FORCES) MATRIX--------------------------*/
  Matrix unit_matrix_kk = Matrix(check_matrix.getRowsCount(), check_matrix.getColumnsCount());

  unit_matrix_kk.getUnitMatrix();

  double value_variable = 0;

  for (int t = 1; t <= num_m; t++) //m <= n
  {
    check_matrix.power_matrix(t);

    value_variable = value_variable + pow(delta_x_val, t);

    value_variable = value_variable / t;

    check_matrix.div(t*value_variable);

    unit_matrix_kk.mSum(check_matrix);
  }

  std::cout << "\nFinal matrix is:\n" << std::endl;

  unit_matrix_kk.print_matrix();

  std::ofstream save_data("data.txt");

  if (!save_data.is_open())
    std::cout << "\nError while creating data.txt" << std::endl;
  else
  {

    for (int i = 0; i < unit_matrix_kk.getRowsCount(); i++)
    {
      for (int j = 0; j < unit_matrix_kk.getColumnsCount(); j++)
      {
        //std::cout << "\nM[" << i + 1 << "][" << j + 1 << "] = " << unit_matrix_kk.getElement(i, j);
        save_data << "\nM[" << i + 1 << "][" << j + 1 << "] = " << unit_matrix_kk.getElement(i, j);
      }
    }

    std::cout << "\nData stored into data.txt" << std::endl;

    // plot data by python api - matplotlib
    CallPyPlotter();

  }
}
/*-------------------MATRIX CALCULATIONS END--------------------------*/



/*--------------------------------DATA INPUT START-------------------------------------*/
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
  } while (flag);
}
/*--------------------------------DATA INPUT END-------------------------------------*/


/*-------------------------------MAIN ENTRY START--------------------------------------*/
int main(int argc, char **argv)
{

  setlocale(LC_ALL, "");
  //if you have an a 1251/866 encoding in yours console output
  //you can use it winapi cheats (check it with chcp cmd):
  //SetConsoleCP(1251);
  //SetConsoleOutputCP(1251);

  //or use own config
  std::streamsize ss = std::cout.precision();

  std::cout << "\nEnter the number of decimal places:\n" << std::endl;

  int k_enter = 0;

  std::cin >> k_enter;

  std::cout.precision(k_enter + 1);

  std::cin.clear();

  enter_data();

  std::cin.get();
}
/*-------------------------------MAIN ENTRY END--------------------------------------*/