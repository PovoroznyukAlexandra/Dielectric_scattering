#include <iostream>
#include <cmath>
#include <chrono>
#include <complex>
#include "diel_scat_func.h"
using namespace std;

int main() {
    // Точка O_sphere - центр шара, R - радиус шара
    double O_sphere[3] = {0, 0 ,0};
    double R = 1;
    cout << "Введите координаты центра шара: (0, 0, 0)\n";
    //cin >> O_sphere[0] >> O_sphere[1] >> O_sphere[2];
    cout << "Введите радиус шара: 1\n";
    //cin >> R;

    // Точка О - центр параллелепипеда
    double O[3] = {0, 0 ,0};
    cout << "Введите координаты центра параллелепипеда: (0, 0, 0)\n";
    //cin >> O[0] >> O[1] >> O[2];


    // L - Длины сторон параллелепипеда
    double L[3] = {2, 2 ,2};
    cout << "Введите длины сторон параллелепипеда (L1, L2, L3 соответственно): (2, 2, 2)" << endl;
    //cin >> L[0] >> L[1] >> L[2];
    if ((L[0] <= 0) || (L[1] <= 0) || (L[2] <= 0)){
        cout << "Стороны параллелепипеда должны быть неотрицательными!" << endl;
        return 1;
    }


    // N - Количество точек на отрезке
    int N[3] = {10, 10 ,10};
    cout << "Введите количество точек на сторонах параллелепипеда (N1, N2, N3 соответственно): (10, 10, 10) \n";
    //cin >> N[0] >> N[1] >> N[2];
    if ((N[0] <= 1) && (N[1] <= 1) && (N[2] <= 1)){
        cout << "Произведение разбиений сторон должно быть больше 1!" << endl;
        return 1;
    }

///--------------------------------------------------------------------------------------------------------------------///

    ///////   Начинаем замер времени   ///////
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

///--------------------------------------------------------------------------------------------------------------------///

    ///////   Нахождение точек коллокации   ///////
    // dx - Длина каждого сегмента на оси
    double dX[3];
    dX[0] = L[0] / N[0];
    dX[1] = L[1] / N[1];
    dX[2] = L[2] / N[2];

    // n - Количество строк в матрице
    int n = N[0] * N[1] * N[2];

    // Создаем и заполняем матрицу colloc_points

    double **colloc_points = new double* [n];    // Создаем массив указателей colloc_points
    for (int i = 0; i < n; i++){
        colloc_points[i] = new double [3];    // Создаем элементы colloc_points
    }

    colloc_points = Area_Discretization(N, L, O);    // заполняем colloc_points

///--------------------------------------------------------------------------------------------------------------------///

    ///////   Считаем значение интеграла   ///////
    double Integral = Count_Integral(n, colloc_points, dX);
    cout << "Значение интеграла I = " << Integral << endl;

///--------------------------------------------------------------------------------------------------------------------///

    ///////   Средний эпсилон на каждой ячейке   ///////
    double sum = 0;
    double *E_cell = new double [n];
    for (int i = 0; i < n; i++){
        E_cell[i] = Eps_on_Cell(colloc_points[i], dX, O_sphere, R);
        sum += E_cell[i];
    }
    for (int i = 0; i < n; i++){
        cout << E_cell[i] << ' ';
    }
    cout << endl;

///--------------------------------------------------------------------------------------------------------------------///

    ///////   Проверка   ///////
    double V = dX[0] * dX[1] * dX[2];
    double E = sum * V;
    cout << "E = " << E << endl;
    cout << "Погрешность Eps: " << abs(E - (8 + 4*M_PI)) << endl; //


///--------------------------------------------------------------------------------------------------------------------///

    ///////   Набор матрицы   ///////
    complex<double> **Matrix = new complex<double>* [3*n];    // Создаем массив указателей Matrix
    for (int i = 0; i < 3*n; i++){
        Matrix[i] = new complex<double> [3*n];    // Создаем элементы Matrix
        for (int j = 0; j < 3*n; j++){
            Matrix[i][j] = 0;
        }
    }
    Get_Matrix(Matrix, colloc_points, E_cell, V, n);

    complex<double> *Matrix_vect = new complex<double> [3*n * 3*n];
    for (int i = 0; i < 3*n; i++){
        for (int j = 0; j < 3*n; j++){
            Matrix_vect[i*3*n + j] = Matrix[j][i];
        }
    }

///--------------------------------------------------------------------------------------------------------------------///

    ///////   Решение СЛАУ   ///////
    complex<double> *g = new complex<double> [3*n]; // Создаем массив указателей
    for (int i = 0; i < 3*n; i++){
        g[i] = 0;
    }
    g = Znach_Func(colloc_points, n);

    cout << "g =" << endl;
    for (int i = 0; i < 3*n; i++){
        cout << g[i] << ' ';
    }
    cout << endl;

    int dim = 3*n;
    int nrhs = 1;
    int LDA = dim;
    int *ipiv = new int[dim];
    int LDB = dim;
    int info;
    zgesv(&dim, &nrhs, Matrix_vect, &LDA, ipiv, g, &LDB, &info);

    cout << "x =" << endl;
    for (int i = 0; i < dim; i++){
        cout << g[i] << ' ';
    }
    cout << endl;
    cout << "Check: info = " << info << endl;


///--------------------------------------------------------------------------------------------------------------------///

    ///////   Проверка   ///////

    cout << endl;
    cout << "K = " << k << endl;
    cout << endl;

    cout << "Check K_func #1" << endl;
    Check_K_func_1();
    cout << endl;

    cout << "Check K_func #2" << endl;
    Check_K_func_2();
    cout << endl;

    cout << "Check collocation points" << endl;
    Check_points(colloc_points, n);
    cout << endl;

    cout << "Check matrix 3nx3n" << endl;
    cout << "For block with I = 35, J = 37:" << endl;
    Check_matrix_3nx3n(35, 37, Matrix);
    cout << "Coordinates of collocation points:" << endl;
    cout << "35: " << colloc_points[35][0] << ' ' << colloc_points[35][1] << ' ' << colloc_points[35][2] << endl;
    cout << "37: " << colloc_points[37][0] << ' ' << colloc_points[37][1] << ' ' << colloc_points[37][2] << endl;
    cout << endl;

    cout << "Check g" << endl;
    Check_g(Vec2Matr(g, n), n);
    cout << endl;

    cout << "Check LAPACK" << endl;
    Check_lapack();
    cout << endl;

///--------------------------------------------------------------------------------------------------------------------///

    ///////   Подсчет ЭПР  ///////
    double sigma[181];
    Effective_Scattering_Area(sigma, colloc_points, Vec2Matr(g, n), E_cell, V, n);
    cout << "ЭПР:" << endl;
    for (int i = 0; i < 181; i++){
        cout << sigma[i] << ' ';
    }
    cout << endl;

    ///////   Перевод ЭПР в Децибелы (dB) ///////
    cout << "ЭПР в dB:" << endl;
    double sigma_dB[181];
    for (int i = 0; i < 181; i++){
        sigma_dB[i] = log10(sigma[i] / (M_PI)) * 10.; // ??????/ как будто надо разделить на 128 и ответ сойдется
        cout << sigma_dB[i] << ", ";
    }
    cout << endl << endl;



///--------------------------------------------------------------------------------------------------------------------///

    ///////   Освобождаем память   ///////
    for (int i = 0; i < 3; i++){
        delete[] colloc_points[i]; // Удаляем каждый элемент
    }
    for (int i = 0; i < 3*n; i++){
        delete[] Matrix[i]; // Удаляем каждый элемент
    }
    delete [] colloc_points; // Удаляем массив
    delete [] E_cell;
    delete [] Matrix;
    delete [] g;
    delete [] ipiv;


    ///////   Вывод времени работы программы   ///////
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Время выполнения: ";
    if (chrono::duration_cast<chrono::milliseconds>(end - begin).count() >= 60000) {
        cout << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << " millisec = ";
        cout << chrono::duration_cast<chrono::minutes>(end - begin).count() << " minutes" << endl;
    }
    else if (chrono::duration_cast<chrono::milliseconds>(end - begin).count() >= 1000){
        cout << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << " millisec = ";
        cout << chrono::duration_cast<chrono::seconds>(end - begin).count() << " sec" << endl;
    }
    else
        cout << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << " millisec" << endl;
    return 0;
}

