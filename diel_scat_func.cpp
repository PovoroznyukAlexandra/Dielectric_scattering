#include "diel_scat_func.h"


///////   Комплексная экспонента   ///////
complex<double> Comp_Exp(double fi){
    complex<double> z(cos(fi), sin(fi));
    //complex<double> z(cos(fi*M_PI/180), sin(fi*M_PI/180));
    return z;
}


///////   Модуль комплексного вектора   ///////
double Comp_abs(const complex<double>* h, int n = 3){
    double comp_abs = 0;
    for (int i = 0; i < n; i++){
        comp_abs += abs(h[i]) * abs(h[i]);
    }
    return sqrt(comp_abs);
}


///////   Скалярное произведение 1   ///////
complex<double> Scal_Prod(const complex<double>* x, complex<double>* y, int n = 3){
    complex<double> res = 0;
    for (int i = 0; i < n; i++){
        res += x[i]*y[i];
    }
    return res;
}


///////   Скалярное произведение 2   ///////
complex<double> Scal_Prod(const complex<double>* x, double* y, int n = 3){
    complex<double> res = 0;
    for (int i = 0; i < n; i++){
        res += x[i]*y[i];
    }
    return res;
}


///////   Скалярное произведение 3   ///////
double Scal_Prod(const double* x, const double* y, int n = 3){
    double res = 0;
    for (int i = 0; i < n; i++){
        res += x[i]*y[i];
    }
    return res;
}


complex<double>** Transpose(complex<double>** A, int m, int n){
    complex<double> **A_T = new complex<double>* [n];
    for (int i = 0; i < n; i++){
        A_T[i] = new complex<double> [m];
    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            A_T[i][j] = A[j][i];
        }
    }
    return A_T;
}


///////   Дискретизация области   ///////
double** Area_Discretization(const int* N, const double *L, const double *O){
    int cnt = 0;
    int n = N[0] * N[1] * N[2];
    double **colloc_points = new double* [n];    // Создаем массив указателей colloc_points
    for (int i = 0; i < n; i++){
        colloc_points[i] = new double [3];    // Создаем элементы colloc_points
        for (int j = 0; j < 3; j++){
            colloc_points[i][j] = 0;
        }
    }
    double dX[3];
    dX[0] = L[0] / N[0];
    dX[1] = L[1] / N[1];
    dX[2] = L[2] / N[2];
    // точка А - левый нижний угол параллелепипеда
    double A[3];
    A[0] = O[0] - L[0]/2;
    A[1] = O[1] - L[1]/2;
    A[2] = O[2] - L[2]/2;

    // заполняем матрицу colloc_points (i - номер точки, j - координат на оси(x, y, z))
    for (int i = 0; i < N[0]; i++){
        for (int j = 0; j < N[1]; j++){
            for (int h = 0; h < N[2]; h++){
                colloc_points[cnt][0] = A[0] + (i + 0.5) * dX[0];
                colloc_points[cnt][1] = A[1] + (j + 0.5) * dX[1];
                colloc_points[cnt][2] = A[2] + (h + 0.5) * dX[2];
                cnt++;
            }
        }
    }
    return colloc_points;
}


///////   Заданная подынтегральная функция   ///////
double F(double **colloc_points, double p){
    return 1;
}


///////   Вычисление интеграла   ///////
double Count_Integral(int n, double **colloc_points, const double *dX){
    double I = 0;
    double p = 1;

    for (int i = 0; i < n; i++){
        I += F(colloc_points, p);
    }
    return (I * dX[0] * dX[1] * dX[2]);
}


///////   Функция Эпсилон   ///////
double Eps(double *x){
    return 4;
}


///////   Функция для Мю   ///////
double Mu(double *x){
    return 1;
}


///////   Функция Хи   ///////
int Chi(const double *x, const double *x_0, double R){
    //if ((x[0] - x_0[0])*(x[0] - x_0[0]) + (x[1] - x_0[1])*(x[1] - x_0[1]) + (x[2] - x_0[2])*(x[2] - x_0[2]) <= R*R)
    double y[3];
    for (int i = 0; i < 3; i++){
        y[i] = x[i] - x_0[i];
    }
    if (Scal_Prod(y, y) <= R*R)
        return 1;
    else
        return 0;
}


///////   Средний эпсилон на ячейке   ///////
double Eps_on_Cell(const double *points, const double *L, const double *O_sphere, double R){
    int M[3];
    M[0] = 30;
    M[1] = 30;
    M[2] = 30;
    int m = M[0] * M[1] * M[2];

    double **one_cell = new double* [m];    // Создаем массив указателей
    for (int i = 0; i < m; i++) {
        one_cell[i] = new double[3];    // Создаем элементы
    }
    one_cell = Area_Discretization(M, L, points);    // заполняем

    double eps = 0;
    for (int i = 0; i < m; i++){
        eps += (Eps(one_cell[i])-1) * Chi(one_cell[i], O_sphere, R);
    }
    return eps/m + 1;
}


///////   Функция К   ///////
complex<double>* K_func (complex<double>* g, const double* x, const double* y){
    double v[3];
    v[0] = x[0] - y[0];
    v[1] = x[1] - y[1];
    v[2] = x[2] - y[2];
    double r = sqrt(Scal_Prod(v, v));
    complex<double>* res = new complex<double>[3];

    complex<double> coeff_g(-1/(r*r*r) + (k*k/r), k/(r*r));
    complex<double> coeff_v(3/(r*r*r) - (k*k/r), -3*k/(r*r));
    for (int i = 0; i < 3; i++){
        res[i] = g[i]*coeff_g + v[i]*Scal_Prod(g, v)*coeff_v/(r*r);
        res[i] *= Comp_Exp(k*r)/(4*M_PI);
    }
    return res;
}


///////   Подсчет одного блока матрицы   ///////
void Count_Matrix_Block(complex<double>** C_ij, int I, int J, double** colloc_points, const double* E_cell, double V){
    complex<double> e_q[3] = {0, 0, 0};
    complex<double> e_p[3] = {0, 0, 0};
    complex<double>* k_res = new complex<double>[3];
    if (I != J){
        for (int p = 0; p < 3; p++){
            e_p[0] = 0; e_p[1] = 0; e_p[2] = 0;
            e_p[p] = 1;
            for (int q = 0; q < 3; q++){
                e_q[0] = 0; e_q[1] = 0; e_q[2] = 0;
                e_q[q] = 1;
                k_res = K_func(e_q, colloc_points[I], colloc_points[J]);
                C_ij[p][q] = Scal_Prod(k_res, e_p) * V * (E_cell[J] - 1);
            }
        }
    }
    else{
        for (int q = 0; q < 3; q++){
            for (int p = 0; p < 3; p++){
                if (p == q){
                    C_ij[q][p] = (-1./3) * (E_cell[I] - 1) - 1;
                }
                else{
                    C_ij[q][p] = 0;
                }
            }
        }
    }
}


///////   Подсчет матрицы   ///////
void Get_Matrix(complex<double>** Matrix, double** colloc_points, const double* E_cell, double V, int n){
    complex<double> **C_ij = new complex<double>* [3];    // Создаем массив указателей Matrix
    for (int i = 0; i < 3; i++){
        C_ij[i] = new complex<double> [3];    // Создаем элементы Matrix
        for (int p = 0; p < 3; p++){
            C_ij[i][p] = 0;
        }
    }
    // Набор матрицы
    for (int I = 0; I < n; I++){
        for (int J = 0; J < n; J++){
            Count_Matrix_Block(C_ij, I, J, colloc_points, E_cell, V);
            for (int q = 0; q < 3; q++){
                for (int p = 0; p < 3; p++){
                    Matrix[I*3 + q][J*3 + p] = C_ij[q][p];
                }
            }
        }
    }
}


///////   Печать матрицы mxn  ///////
void Print_Matrix(complex<double>** Matrix, int m, int n, const string& str = "\n"){
    cout << "Матрица: " << str << endl;
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            cout << Matrix[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
}


///////   Подсчет правой части СЛАУ  ///////
complex<double>* Znach_Func(double** colloc_points, int n){
    complex<double> *E_inc = new complex<double> [3*n];
    for (int i = 0; i < n; i++){
        E_inc[3*i + 0] = 0;
        E_inc[3*i + 1] = -Comp_Exp(-k*colloc_points[i][0]); // ?????? E_inc
        E_inc[3*i + 2] = 0;
    }
    return E_inc;
}

///////   Перевод вектора в матрицу  ///////
complex<double>** Vec2Matr(complex<double> *x, int n){
    complex<double> **A = new complex<double> *[n];
    for (int i = 0; i < n; i++){
        A[i] = new complex<double> [3];
        for (int p = 0; p < 3; p++){
            A[i][p] = x[3*i + p];
        }
    }
    return A;
}


///////   Подсчет ЭПР   ///////
void Effective_Scattering_Area(double* sigma, double** colloc_points, complex<double>** g, const double* E_cell, double V, int n){
    double tau[181][3];
    for (int alpha = 0; alpha < 181; alpha++){
        tau [alpha][0] = cos(alpha*M_PI/180);
        tau [alpha][1] = sin(alpha*M_PI/180);
        tau [alpha][2] = 0;
    }

    //complex<double> sigma[181];
    complex<double> h[181][3];
    double E_inc = 1;
    for (int i = 0; i < 181; i++){
        for (int j = 0; j < 3; j++){
            h[i][j] = 0;
        }
    }
    /*
    for (int alpha = 0 ; alpha < 181; alpha++){
        for (int p = 0; p < 3; p++){
            for (int i = 0; i < n; i++){
                complex<double> exp = Comp_Exp(-k * Scal_Prod(tau[alpha], colloc_points[i]));
                complex<double> coeff = exp * (E_cell[i] - 1) * k * k;
                h[alpha][p] += coeff * (g[i][p] - Scal_Prod(g[i], tau[alpha]) * tau[alpha][p]) * 0.008;  //0.008 - объем одной ячейки
            }
        }
        sigma[alpha] = Comp_abs(h[alpha]) * Comp_abs(h[alpha]) / (E_inc * 4 * M_PI);
    }
     */
    complex<double> exp, coeff;
    for (int alpha = 0 ; alpha < 181; alpha++){
        for (int i = 0; i < n; i++){
            exp = Comp_Exp(-k * Scal_Prod(tau[alpha], colloc_points[i]));
            coeff = exp * (E_cell[i] - 1) * k * k * V;
            for (int p = 0; p < 3; p++){
                h[alpha][p] += coeff * (g[i][p] - Scal_Prod(g[i], tau[alpha]) * tau[alpha][p]);
            }
        }
        sigma[alpha] = Comp_abs(h[alpha]) * Comp_abs(h[alpha]) / (E_inc * 4 * M_PI);
    }
}

void Check_K_func_1(){
    complex<double> vec0(1, 1);
    complex<double> vec1(2, -1);
    complex<double> vec2(-1, 1);
    complex<double> vec[3] = {vec0, vec1, vec2};
    complex<double> *res_vec = new complex<double>[3];
    double x_vec[3] = {-1, 2, 3};
    double y_vec[3] = {1, -1, 2};

    res_vec = K_func(vec, x_vec, y_vec);
    cout << "For g = (1+i, 2-i, -1+i),  x = (-1, 2, 3), y = (1, -1, 2) :" << endl;
    cout << res_vec[0] << ' ' << res_vec[1] << ' ' << res_vec[2] << ' ' << endl;

}

void Check_K_func_2(){
    complex<double> vec[3] = {1, 0, 0};
    complex<double> *res_vec = new complex<double>[3];
    double x_vec[3] = {0, 0, 0};
    double y_vec[3] = {1, 1, 1};

    res_vec = K_func(vec, x_vec, y_vec);
    cout << "For g = (1, 0, 0),  x = (0, 0, 0), y = (1, 1, 1) :" << endl;
    cout << res_vec[0] << ' ' << res_vec[1] << ' ' << res_vec[2] << ' ' << endl;

}

void Check_points(double** colloc_points, int n){
    std::ofstream out;          // поток для записи
    out.open("Grid.gr");      // открываем файл для записи
    if (out.is_open()){
        out << "2\t" << n/2 << endl;
        for (int i = 0; i < n; i++){
            out << colloc_points[i][0] << "\t" << colloc_points[i][1] << "\t" << colloc_points[i][2] << std::endl;
        }
    }
    out.close();
    std::cout << "File Grid.gr has been written" << std::endl;
}

void Check_matrix_3nx3n(int I, int J, complex<double> **matrix){
    for (int p = 0; p < 3; p++){
        for (int q = 0; q < 3; q++){
            cout << matrix[3*I+p][3*J+q] << ' ';
        }
        cout << endl;
    }
}

void Check_g(complex<double>** g, int n){
    ofstream out;          // поток для записи
    out.open("Re_g.gv");      // открываем файл для записи
    if (out.is_open()){
        out << "2\t" << n/2 << endl;
        for (int i = 0; i < n; i++){
            out << g[i][0].real() << " \t " << g[i][1].real() << " \t " << g[i][2].real() << endl;
        }
    }
    out.close();
    cout << "File Re_g.gv has been written" << endl;

    out.open("Im_g.gv");      // открываем файл для записи
    if (out.is_open()){
        out << "2\t" << n/2 << endl;
        for (int i = 0; i < n; i++){
            out << g[i][0].imag() << " \t " << g[i][1].imag() << " \t " << g[i][2].imag() << endl;
        }
    }
    out.close();
    cout << "File Im_g.gv has been written" << endl;
}

void Check_lapack(){
    int dim = 3;
    int nrhs = 1;
    int LDA = dim;
    int *ipiv = new int[dim];
    int LDB = dim;
    int info;

    complex<double> *g = new complex<double>[3];
    g[0] = 0; g[1] = -7; g[2] = 13;

    complex<double> *Matrix_vect = new complex<double>[9];
    Matrix_vect[0] = 1; Matrix_vect[1] = 1; Matrix_vect[2] = 2;
    Matrix_vect[3] = 1; Matrix_vect[4] = -1; Matrix_vect[5] = 5;
    Matrix_vect[6] = -1; Matrix_vect[7] = -2; Matrix_vect[8] = 1;
    zgesv(&dim, &nrhs, Matrix_vect, &LDA, ipiv, g, &LDB, &info);

    cout << "x =" << endl;
    for (int i = 0; i < dim; i++){
        cout << g[i] << ' ';
    }
    cout << endl;
    cout << "Check: info = " << info << endl;

}
