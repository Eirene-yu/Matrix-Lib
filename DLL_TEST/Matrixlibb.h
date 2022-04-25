#pragma once

#define DLL_API __declspec (dllexport)

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <initializer_list>

using namespace std;

class Matrix {
private:
    int N_line;
    int N_column;
    vector<vector<double>> MATRIX;
public:
    // конструкторы
    DLL_API Matrix() {
        N_line = 0;
        N_column = 0;
        MATRIX.assign(0, vector<double>(0));
    };
    // c параметрами
    DLL_API Matrix(int n, int m) {
        N_line = n;
        N_column = m;
        MATRIX.assign(n, vector<double>(m));
    };
    //c initializer_list<double> - матрица 1 * m
    DLL_API Matrix(initializer_list<int> list) {
        N_line = 1;
        N_column = list.size();
        MATRIX.assign(N_line, vector<double>(N_column));
        int count = 0;
        for (auto element : list)
        {
            MATRIX[0][count] = element;
            ++count;
        }
    };
    //c initializer_list<initializer_list<double>> - матрица n * m
    DLL_API Matrix(initializer_list<initializer_list<double>> list) {
        N_line = list.size();
        if (N_line == 0)throw runtime_error("Нельзя создать пустую матрицу");
        N_column = 0; int count = 0;
        for (auto element : list) {
            if (N_column == 0) N_column = element.size();
            for (auto el : element) ++count;
            if (N_column != count)throw runtime_error("Нельзя создать матрицу, в которой разное количество элементов в строках!");
            count = 0;
        }
        MATRIX.assign(N_line, vector<double>(N_column));
        int col = 0, row = 0;
        for (auto element : list) {
            for (auto el : element) {
                MATRIX[row][col] = el;
                ++col;
            }
            if (col == N_column) {
                col = 0;
                ++row;
            }
        }
    };
    //количество строк
    DLL_API int Num_line() const {
        return N_line;
    }
    //количество столбцов
    DLL_API int Num_column() const {
        return N_column;
    }
    //элемент на месте i, j
    DLL_API const double& Element_ij(int i, int j) const {
        return MATRIX.at(i).at(j);
    }

    DLL_API double& Element_ij(int i, int j) {
        return MATRIX[i][j];
    }

    DLL_API const vector<vector<double>>& GetMatrix() const {
        return MATRIX;
    }
    DLL_API vector<vector<double>>& GetMatrix() {
        return MATRIX;
    }
    DLL_API int& Num_line() {
        return N_line;
    }
    DLL_API int& Num_column() {
        return N_column;
    }
    //1 часть
    DLL_API friend Matrix operator+(const Matrix& a, const Matrix& b);//сложение матриц
    DLL_API friend Matrix operator-(const Matrix& a, const Matrix& b);//разность матриц
    DLL_API friend Matrix operator*(const Matrix& a, const Matrix& b);// произведение матриц
    DLL_API friend Matrix operator*(const Matrix& a, const int b); //умножение на число
    DLL_API friend Matrix Hadamart(const Matrix& a, const Matrix& b);//произведение Адамара
    DLL_API friend ostream& operator<< (ostream&, const Matrix&);
    DLL_API friend istream& operator>>(istream&, Matrix&);

    //транспонированная матрица
    DLL_API friend Matrix Transposition(const Matrix& A);
};

//Подкласс единичная
class OneMatrix : public Matrix {
public:
    DLL_API OneMatrix(int N) : Matrix(N, N) {
        for (int i = 0; i < N; ++i) {
            this->Element_ij(i, i) = 1;
        }
    }
};

//Подкласс диагональная
class DiagonalMatrix : public OneMatrix {
public:
    DLL_API DiagonalMatrix(int N) : OneMatrix(N) {
        for (int i = 0; i < N; ++i) {
            cout << "Ввод элемента: " << "[" << i << "]" << "[" << i << "]" << " ";
            double el;
            cin >> el;
            this->Element_ij(i, i) = el;
        }
    }
};

//Подкласс верхняя треугольная
class UpperTriangularMatrix : public Matrix {
public:
    DLL_API UpperTriangularMatrix(int N) : Matrix(N, N) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (j >= i) {
                    cout << "Ввод элемента: " << "[" << i << "]" << "[" << j << "]" << " ";
                    double el;
                    cin >> el;
                    this->Element_ij(i, j) = el;
                }
            }
        }
    }
};

//Подкласс нижняя треугольная
class LowerTriangularMatrix : public Matrix {
public:
    DLL_API LowerTriangularMatrix(int N) : Matrix(N, N) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (j <= i) {
                    cout << "Ввод элемента: " << "[" << i << "]" << "[" << j << "]" << " ";
                    double el;
                    cin >> el;
                    this->Element_ij(i, j) = el;
                }
            }
        }
    }
};

class SymmetricalMatrix : public Matrix {
public:
    DLL_API SymmetricalMatrix(int N) : Matrix(N, N) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double el;
                if (j == i) {
                    cout << "Ввод элемента: " << "[" << i << "]" << "[" << j << "]" << " ";
                    double el;
                    cin >> el;
                    this->Element_ij(i, j) = el;
                }
                else if (i > j) {
                    cout << "Ввод элементов: " << "[" << i << "]" << "[" << j << "]" << " и " << "[" << j << "]" << "[" << i << "]" << " ";
                    cin >> el;
                    this->Element_ij(i, j) = el;
                    this->Element_ij(j, i) = el;
                }
            }
        }
    }
};

//равенство матриц
DLL_API bool operator==(const Matrix& a, const Matrix& b);
//неравенство матриц
DLL_API bool operator!=(const Matrix& a, const Matrix& b);

//2 часть
//след матрицы
DLL_API double Trace(const Matrix& a);
//скалярное произведение векторов
DLL_API double Scalar_product(const Matrix& a, const Matrix& b);
//евклидова норма вектора
DLL_API double Euclidean_norm(const Matrix& a);
//максимальная норма вектора
DLL_API double Maximum_norm(const Matrix& a);
//Норма матрицы (норма Фробениуса)
DLL_API double Frobenius_norm(const Matrix& a);
//определитель матрицы используя метод Гаусса
DLL_API double Determinant(const Matrix& a);
//обработка по методу Гаусса верхнетреугольная матрица
DLL_API void Gauss_method(Matrix& a, int cur, int& num_permutations);

//3 часть
//угол между векторами
DLL_API double Angle(const Matrix& a, const Matrix& b);//транспонированная матрица
DLL_API Matrix Transposition(const Matrix& A);
//ранг матрицы используя алгоритм Гаусса
DLL_API int Rank(const Matrix& A);
//обратная матрица
DLL_API Matrix Reverse(const Matrix& A);
//для обнуления треугольника над диагональю
DLL_API void Lower_Gauss(Matrix& a, int cur);

//4 часть
//запись в текстовый файл
DLL_API void Write_text_file(const string& file, const Matrix& M);
//чтение из текстового файла
DLL_API Matrix Read_text_file(const string& file);
//запись в бинарный файл
DLL_API void Write_bin_file(const string& file, const Matrix& M);
//чтение из бинарного файла
DLL_API Matrix Read_bin_file(const string& file);


class PCA {
private:
    Matrix MATRIX, P, T, E;
    bool nip;
public:
    DLL_API PCA(const Matrix& a) {
        MATRIX.Num_line() = a.Num_line();
        MATRIX.Num_column() = a.Num_column();
        MATRIX.GetMatrix() = a.GetMatrix();
    }
    DLL_API const Matrix& GetPCA() const {
        return MATRIX;
    }
    DLL_API const Matrix& GetP() const {
        return P;
    }
    DLL_API const Matrix& GetT() const {
        return T;
    }
    DLL_API const Matrix& GetE() const {
        return E;
    }
    DLL_API void Centering();
    DLL_API void Scaling();
    DLL_API tuple<Matrix, Matrix, Matrix> NIPALS(size_t A = -10, double epsilon = 0.00000001);
    DLL_API void PRINT();
    DLL_API vector<double> Scope() const;
    DLL_API vector<double> Deviations() const;
    DLL_API double Average_distance();
    DLL_API vector<double> Dispersion();
    DLL_API double Full_dispersion();
    DLL_API double Explained_dispersion();
};
