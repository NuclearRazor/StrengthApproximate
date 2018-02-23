#pragma once
#define EIGEN_RUNTIME_NO_MALLOC // Define this symbol to enable runtime tests for allocations
#define _USE_MATH_DEFINES //preprocessor directive to use mathematical constants
#include <cmath> 
#include <iostream> 
#include <fstream>
#include <Windows.h>
#include <C:\LIBS_ETC\vcpkg\packages\eigen3_x86-windows\include\eigen3\Eigen\Dense>
#include <C:\LIBS_ETC\vcpkg\packages\eigen3_x86-windows\include\eigen3\unsupported\Eigen\CXX11\Tensor>

/*--------------------FONT COLOR--------------------------*/
BOOL PrintText(LPCTSTR szText);
BOOL SetConsoleAttrib(WORD wAttrib);
BOOL SetCurrentPos(SHORT x, SHORT y);
BOOL SetConsoleSize(SHORT x, SHORT y);

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
//using Eigen::MatrixXd;
const int DIM_A = 8;

//functions declaration
void enter_data();
void evaluate_matrix(double& radius, double& x_val, double& phi_val, double& h_val, double& E_mod_val, double& nu_val, double& delta_phi_val, double& delta_x_val);
char GetInput();
void DisplayMainMenu();
void save_data(MatrixXd unit_matrix_kk);