#ifndef DIELECTRIC_SCATTERING_DIEL_SCAT_FUNC_H
#define DIELECTRIC_SCATTERING_DIEL_SCAT_FUNC_H

#include <iostream>
#include <cmath>
#include <chrono>
#include <complex>
#include <fstream>
using namespace std;

///////   Константы   //////
const int C = 299792458 ;    // скорость света
const double omega = 0.15 * (1e+9) * 2 * M_PI;
const double k = omega/C;
//const double EPS_0 = 8.85 * (0.000000000001);    // электрическая постоянная ??????????
//const double MU_0 = 4 * M_PI * (0.0000001);    // магнитная постоянная ??????????
//const double k = omega * sqrt(EPS_0 * MU_0);

///////   Комплексная экспонента   ///////
complex<double> Comp_Exp(double fi);

///////   Модуль комплексного вектора   ///////
double Comp_abs(const complex<double>* h, int n);

///////   Скалярное произведение 1   ///////
complex<double> Scal_Prod(const complex<double>* x, complex<double>* y, int n);

///////   Скалярное произведение 2   ///////
complex<double> Scal_Prod(const complex<double>* x, double* y, int n);

///////   Скалярное произведение 3   ///////
double Scal_Prod(const double* x, const double* y, int n);

complex<double>** Transpose(complex<double>** A, int m, int n);

///////   Дискретизация области   ///////
double** Area_Discretization(const int* N, const double *L, const double *O);

///////   Заданная подынтегральная функция   ///////
double F(double **colloc_points, double p);

///////   Вычисление интеграла   ///////
double Count_Integral(int n, double **colloc_points, const double *dX);

///////   Функция Эпсилон   ///////
double Eps(double *x);

///////   Функция для Мю   ///////
double Mu(double *x);

///////   Функция Хи   ///////
int Chi(const double *x, const double *x_0, double R);

///////   Средний эпсилон на ячейке   ///////
double Eps_on_Cell(const double *points, const double *L, const double *O_sphere, double R);

///////   Функция К   ///////
complex<double>* K_func (complex<double>* g, const double* x, const double* y);

///////   Подсчет одного блока матрицы   ///////
void Count_Matrix_Block(complex<double>** C_ij, int I, int J, double** colloc_points, const double* E_cell, double V);

///////   Подсчет матрицы   ///////
void Get_Matrix(complex<double>** Matrix, double** colloc_points, const double* E_cell, double V, int n);

///////   Печать матрицы mxn  ///////
void Print_Matrix(complex<double>** Matrix, int m, int n, const string& str);

///////   Подсчет правой части СЛАУ  ///////
complex<double>* Znach_Func(double** colloc_points, int n);

///////   Перевод вектора в матрицу  ///////
complex<double>** Vec2Matr(complex<double> *x, int n);

///////   Решение СЛАУ   ///////
extern "C" void zgesv(int *n, int *nrhs, complex<double>  *A, int *lda, int* ipiv, complex<double> *b, int *ldb, int *info);

///////   Подсчет ЭПР  ///////
void Effective_Scattering_Area(double* sigma, double** colloc_points, complex<double>** g, const double* E_cell, double V, int n);

void Check_K_func_1();

void Check_K_func_2();

void Check_points(double** colloc_points, int n);

void Check_matrix_3nx3n(int I, int J, complex<double> **matrix);

void Check_g(complex<double>** g, int n);

void Check_lapack();

#endif //DIELECTRIC_SCATTERING_DIEL_SCAT_FUNC_H
