#pragma once
#define _USE_MATH_DEFINES //preprocessor directive to use mathematical constants
#include <cmath> 
#include <iostream> 
#include <fstream>
#include <vector>
#include <numeric>
#include <Windows.h>

//Call for use implementation of Matrix class
#include "Matrix.h"

/*--------------------FONT COLOR WIN API CONFIG START--------------------------*/
BOOL PrintText(LPCTSTR szText);
BOOL SetConsoleAttrib(WORD wAttrib);
BOOL SetCurrentPos(SHORT x, SHORT y);
BOOL SetConsoleSize(SHORT x, SHORT y);
/*--------------------FONT COLOR WIN API CONFIG END--------------------------*/


//functions declaration
void enter_data();
void evaluate_matrix(double& radius, double& x_val, double& phi_val, double& h_val, double& E_mod_val, double& nu_val, double& delta_phi_val, double& delta_x_val);
char GetInput();
void DisplayMainMenu();