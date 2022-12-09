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
double E_hec();
double U();
double z_m();
double update_z_m(double new_value);
int fN_max();
int fN_c();
int update_N_max(int new_value);
double E_d();
double V_0();
double update_V(double new_value);
double E_uv();
double W_1();
double W_2(); 
double E_ref();
double update_Eref(double new_value);
