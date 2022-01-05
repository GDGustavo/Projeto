// param.h

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <iterator>
#include <regex>

double read_param(std::string var);
int write_param(std::string var, float var_value);
int update_parameters();
void read_all_params();
double lamb();
double Gamma();
double update_Gamma(double new_value);
double U();
double z_m();
double update_z_m(double new_value);
int fN_max();
int update_N_max(int new_value);
double E_d();
double V_0();
double E_uv();
double W_1();
double W_2(); 
