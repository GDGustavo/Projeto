#include "param.h"

static double Lamb_ = 0.0; // Discretization parameter Lambda, to be read from file "param.txt"
static double E_hec_ = 0.0; // 
static double z_m_ = 0.0; // 
static double U_ = 0.0; // 
static int N_max_ = 0; //
static int N_c_ = 0; 
static double E_d_ = 0.0; // 
static double E_uv_ = 0.0; //  
static double W_1_ = 0.0; // 
static double W_2_ = 0.0; // 
static double V_0_ = 0.0; //
static double E_ref_ = 0.0; //  

double  read_param(std::string var){

  double var_read = 0;
  std::string line;

  std::ifstream file;
  file.open("parameters.txt");
  
  //Reading on file and Processing the data. 
  
  while (std::getline(file, line)) {
  std::string cstr = line; 
  std::string s = var + "==";
  std::regex e ("^([a-zA-Z_]\\w*)\\s*==\\s*([^#]+?)$");
  std::smatch cm;
  std::regex_match (cstr, cm, e, std::regex_constants::match_default );
  
  if(var == cm[1]){
  var_read = (std::stod(cm[2]));}
  }
 
  file.close(); 
  return var_read;
}

//Modifying the value of a variable var in parameters.txt
//var_value = var_value + 1;
//write_param(var, var_value); 
//update_parameters(); 
//printf("%.3lf\n",read_param(var));

int write_param(std::string var, float var_value){

  std::ifstream file;
  file.open("parameters.txt");
  std::ofstream file_aux;
  file_aux.open("aux.txt");
  
  //Reading on file and Processing the data. 
  std::string line; 
  std::string str1;
  
  while (std::getline(file, line)) {
  std::string cstr = line; 
  std::string s = var + "==";
  std::regex e ("^([a-zA-Z_]\\w*)\\s*==\\s*([^#]+?)$");
  std::smatch cm;
  std::regex_match (cstr, cm, e, std::regex_constants::match_default);
 
  str1 = cm[2];     
  if(var == cm[1]){
  str1  = std::to_string(var_value);
  }
 
  file_aux << cm[1];
  file_aux << "==";
  file_aux << str1;
  file_aux << std::endl;
  }
 
  //Closing the file
  file.close();
  file_aux.close();
  printf("Writing in file part 1: Done!\n");
  printf("Don't forget to update the parameters.txt file\n");
  return 0;
}


int update_parameters(){

  std::ifstream file;
  file.open("aux.txt");
  std::ofstream file_aux;
  file_aux.open("parameters.txt");
  
  //Reading on file and Processing the data. 
  std::string line;
  std::string var ("var");
  std::string str1; 
  
  while (std::getline(file, line)) {
  std::string cstr = line; 
  std::string s = var + "==";
  std::regex e ("^([a-zA-Z_]\\w*)\\s*==\\s*([^#]+?)$");
  std::smatch cm;
  std::regex_match (cstr, cm, e, std::regex_constants::match_default );
 
  str1 = cm[2];
   
  file_aux << cm[1];
  file_aux << "==";
  file_aux << str1;
  file_aux << std::endl;
  }
 
  file.close();
  file_aux.close();
  printf("Update data in parameters.txt: Done!\n");
  return 0;
}


//Defining functions for read all importants Variables 

void read_all_params(){

std::cout << "Reading all parameters..." << std::endl;

Lamb_ = read_param("Lambda");
std::cout << "Lambda = "<< Lamb_ << std::endl;

//Gamma_ = read_param("Gamma");
//std::cout << "Gamma = "<< Gamma_ << std::endl;

N_max_ = (int) static_cast<int>(read_param("N_max"));
std::cout << "N_max="<< '\t' << N_max_ << std::endl;

N_c_ = (int) static_cast<int>(read_param("N_c"));
std::cout << "N_c="<< '\t' << N_c_ << std::endl;

E_d_ = read_param("E_d");
std::cout << "E_d=" << '\t'<< E_d_ << std::endl;

U_ = read_param("U");
std::cout << "U="<< '\t' << U_ << std::endl;

V_0_ = read_param("V_0");
std::cout << "V_0="<< '\t' << V_0_ << std::endl;

W_1_ = read_param("W_1");
std::cout << "W_1="<<  '\t' << W_1_ << std::endl;

W_2_ = read_param("W_2");
std::cout << "W_2="<< '\t' << W_2_ << std::endl;

z_m_ = read_param("z_m");
std::cout << "z_m=" << '\t' << z_m_ << std::endl;

E_uv_ = read_param("E_uv");
std::cout << "E_uv=" << '\t' << E_uv_ << std::endl;

E_hec_ = read_param("E_hec");
std::cout << "E_hec=" << '\t' << E_hec_ << std::endl;

E_ref_ = read_param("E_ref");
std::cout << "E_ref=" << '\t' << E_ref_ << std::endl;

std::cout << "Done! All parameters have already been read!" << std::endl;

}


// Defining the functions var and update var for each variable. 

//Lambda
double lamb()
{
      return Lamb_;
}

//E_ref
double E_ref()
{
      return E_ref_;
}

double update_Eref(double new_value)
{
     E_ref_ = new_value;
     return E_ref_;
}

//U
double U()
{
      	return U_;
}

//z_m
double z_m()
{
      return z_m_;
}

double update_z_m(double new_value)
{
     z_m_ = new_value;
     return z_m_;
}

//N
int fN_max()
{
	N_max_ = (int) static_cast<int>(read_param("N_max"));
	return N_max_;
}

int fN_c()
{
	N_c_ = (int) static_cast<int>(read_param("N_c"));
	return N_c_;
}

//E_d
double E_d()
{
      return E_d_;
}


//V_0
double V_0()
{
      return V_0_;
}

double update_V(double new_value)
{
     V_0_ = new_value;
     return V_0_;
}

//E_uv
double E_uv()
{
      return E_uv_;
}

//E_hec
double E_hec()
{
	return E_hec_;
}

// W_2
double W_1()
{
      return W_1_;
}

// W_1
double W_2()
{
      return W_2_;
}
