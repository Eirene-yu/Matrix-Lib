#include "pch.h" // use stdafx.h in Visual Studio 2017 and earlier
#include "Matrixlibb.h"

#include <cmath>
#include <algorithm> 
#include <fstream>
#include <tuple>

#define PI 3.14159265

//перегрузка оператора вывода
ostream& operator<< (ostream& output, const Matrix& O)
{
    for (int i = 0; i < O.Num_line(); ++i)
    {
        for (int j = 0; j < O.Num_column(); ++j)
        {
            output << O.Element_ij(i, j);
            if (j != O.Num_column() - 1) output << ' ';
        }
        output << '\n';
    }
    return output;
}

//перегрузка оператора считывания матрицы
istream& operator>>(istream& is, Matrix& O)
{
    for (int i = 0; i < O.Num_line(); ++i)
    {
        for (int j = 0; j < O.Num_column(); ++j)
        {
            double current_num;
            is >> current_num;
            O.MATRIX[i][j] = current_num;
            is.ignore(1);
        }
    }
    return is;
}

//оператор проверки равенства матриц
bool operator==(const Matrix& a, const Matrix& b) {
    if (a.Num_line() != b.Num_line() || a.Num_column() != b.Num_column()) {
        throw false;
    }
    for (int i = 0; i < a.Num_line(); ++i) {
        for (int j = 0; j < b.Num_line(); ++j) {
            if (a.Element_ij(i, j) != b.Element_ij(i, j))
                return false;
        }
    }
    return true;
}

//оператор проверки неравенства матриц
bool operator!=(const Matrix& a, const Matrix& b) {
    return !(a == b);
}

//перегрузка оператора сложения
Matrix operator+(const Matrix& a, const Matrix& b) {
    if (a.Num_line() != b.Num_line() || a.Num_column() != b.Num_column()) throw runtime_error("Матрицы нельзя сложить, так как они разной размерности.");
    int N = a.Num_line(), M = a.Num_column();
    Matrix result(N, M);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            result.MATRIX[i][j] = a.Element_ij(i, j) + b.Element_ij(i, j);
        }
    }
    return result;
}

//перегрузка оператора вычитания
Matrix operator-(const Matrix& a, const Matrix& b) {
    if (a.Num_line() != b.Num_line() || a.Num_column() != b.Num_column()) throw runtime_error("Нельзя найти разность матриц, так как у них разные размерности.");
    int N = a.Num_line(), M = a.Num_column();
    Matrix result(N, M);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            result.MATRIX[i][j] = a.Element_ij(i, j) - b.Element_ij(i, j);
        }
    }
    return result;
}

//перегрузка оператора умножения
Matrix operator*(const Matrix& a, const Matrix& b) {
    if (a.Num_column() != b.Num_line()) throw runtime_error("Нельзя найти произведение матриц так как количество столбцов первой матрице не равно количеству строк во второй.");
    else {
        int N = a.Num_line(), M = b.Num_column();
        Matrix result(N, M);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                double res = 0;
                for (int k = 0; k < a.Num_column(); k++) {
                    res += a.Element_ij(i, k) * b.Element_ij(k, j);
                }
                if (abs(res) < 0.000001) res = 0;
                result.MATRIX[i][j] = res;
            }
        }
        return result;
    }
}

//перегрузка оператора умножения на число
Matrix operator*(const Matrix& a, const int b) {
    int N = a.Num_line(), M = a.Num_column();
    Matrix result(N, M);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            result.MATRIX[i][j] = b * a.Element_ij(i, j);
        }
    }
    return result;
}

//произведение Адамара
Matrix Hadamart(const Matrix& a, const Matrix& b) {
    if (a.Num_column() != b.Num_column() || a.Num_line() != b.Num_line()) throw runtime_error("Нельзя найти так как матрицы имеют разную размерность.");
    int N = a.Num_line(), M = a.Num_column();
    Matrix result(N, M);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            result.MATRIX[i][j] = a.Element_ij(i, j) * b.Element_ij(i, j);
        }
    }
    return result;
}

//след матрицы
double Trace(const Matrix& a) {
    if (a.Num_line() != a.Num_column()) throw runtime_error("Нельзя найти след, так как матрица не квадратная, количество столбцов не равно количеству строк.");
    double trace = 0;
    for (int i = 0; i < a.Num_line(); i++)
    {
        for (int j = 0; j < a.Num_column(); j++)
        {
            if (i == j) trace += a.Element_ij(i, j);
        }
    }
    return trace;
}

//скалярное произведение векторов
double Scalar_product(const Matrix& a, const Matrix& b) {
    if ((a.Num_line() != 1 && a.Num_column() != 1) || (b.Num_line() != 1 && b.Num_column() != 1)) throw runtime_error("Нельзя найти скалярное произведение, так как один из элементов не вектор.");
    double scalar = 0;
    if ((a.Num_line() == b.Num_line() && a.Num_column() == b.Num_column()) || (a.Num_line() == b.Num_column() && b.Num_line() == a.Num_column())) {
        vector<double> v1, v2;
        for (int i = 0; i < a.Num_line(); i++)
            for (int j = 0; j < a.Num_column(); j++)
                v1.push_back(a.Element_ij(i, j));
        for (int i = 0; i < b.Num_line(); i++)
            for (int j = 0; j < b.Num_column(); j++)
                v2.push_back(b.Element_ij(i, j));
        int N1 = a.Num_line();
        if (N1 < a.Num_column()) N1 = a.Num_column();
        for (int i = 0; i < N1; i++) {
            scalar += v1[i] * v2[i];
        }
    }
    else throw runtime_error("Нельзя найти скалярное произведение, так как вектора имеют разную длину.");
    return scalar;
}

//евклидова норма вектора
double Euclidean_norm(const Matrix& a) {
    if (a.Num_line() != 1 && a.Num_column() != 1) throw runtime_error("Нельзя найти евклидову норму вектора, так как это матрица, а не вектор.");
    double norm = 0;
    for (int i = 0; i < a.Num_line(); i++)
    {
        for (int j = 0; j < a.Num_column(); j++)
        {
            norm += a.Element_ij(i, j) * a.Element_ij(i, j);
        }
    }
    return sqrt(norm);
}

//максимальная норма вектора
double Maximum_norm(const Matrix& a) {
    if (a.Num_line() != 1 && a.Num_column() != 1) throw runtime_error("Нельзя найти максимальную норму вектора, так как это матрица, а не вектор.");
    double norm = 0;
    for (int i = 0; i < a.Num_line(); i++)
    {
        for (int j = 0; j < a.Num_column(); j++)
        {
            if (norm < abs(a.Element_ij(i, j))) norm = abs(a.Element_ij(i, j));
        }
    }
    return norm;
}

//Норма матрицы (норма Фробениуса)
double Frobenius_norm(const Matrix& a) {
    double F = 0;
    for (int i = 0; i < a.Num_line(); i++)
        for (int j = 0; j < a.Num_column(); j++)
            F += a.Element_ij(i, j) * a.Element_ij(i, j);
    return sqrt(F);
}

//обработка по методу Гаусса верхнетреугольная матрица
void Gauss_method(Matrix& a, int cur, int& num_permutations) {
    int ROW = a.Num_line(), COL = a.Num_column(), i_max = -1;
    auto& M = a.GetMatrix(); double max_v;
    bool f = true;
    for (int i = cur; i < ROW; i++) {
        if (f) {
            i_max = i;
            max_v = a.Element_ij(i, cur);
            f = false;
        }
        else {
            if (a.Element_ij(i, cur) > max_v) {
                max_v = a.Element_ij(i, cur);
                i_max = i;
            }
        }
    }
    if (cur != i_max) { //перемещаем строки
        swap(M[cur], M[i_max]);
        num_permutations++;
    }
    for (int i = cur + 1; i < ROW; i++) {
        if (M[cur][cur] == 0)
            break;
        double coeff = M[i][cur] / M[cur][cur];
        for (int j = 0; j < COL; j++) {
            M[i][j] -= M[cur][j] * coeff;
            if (abs(M[i][j] - 0) < 0.00001) M[i][j] = 0;
        }
    }
}

//определитель матрицы используя метод Гаусса
double Determinant(const Matrix& a) {
    if (a.Num_column() != a.Num_line()) throw runtime_error("Нельзя посчитать определитель методом Гаусса, так как матрица не квадратная.");
    if (a.Num_column() == 0 || a.Num_line() == 0) throw runtime_error("Нельзя посчитать определитель для пустой матрицы.");
    int N = a.Num_line(), num_permutations = 0;
    Matrix A = a;
    //вызов обработки строк
    for (int i = 0; i < N; i++) {
        Gauss_method(A, i, num_permutations);
    }
    double determ = 1;
    for (int i = 0; i < A.Num_line(); i++) {
        determ *= A.Element_ij(i, i);
    }
    if (determ == 0) return determ;
    else {
        return (pow(-1, num_permutations) * determ);
    }
}

//угол между векторами
double Angle(const Matrix& a, const Matrix& b) {
    double rad = 0;
    if ((a.Num_line() != 1 && a.Num_column() != 1) || (b.Num_line() != 1 && b.Num_column() != 1)) throw runtime_error("Нельзя найти угол, так как один из элементов не вектор.");
    if ((a.Num_line() == b.Num_line() && a.Num_column() == b.Num_column()) || (a.Num_line() == b.Num_column() && b.Num_line() == a.Num_column())) {
        vector<double> v1, v2;
        for (int i = 0; i < a.Num_line(); i++)
            for (int j = 0; j < a.Num_column(); j++)
                v1.push_back(a.Element_ij(i, j));
        for (int i = 0; i < b.Num_line(); i++)
            for (int j = 0; j < b.Num_column(); j++)
                v2.push_back(b.Element_ij(i, j));
        int N1 = a.Num_line();
        if (N1 < a.Num_column()) N1 = a.Num_column();
        double num = 0, den1 = 0, den2 = 0;
        for (int i = 0; i < N1; i++) {
            num += v1[i] * v2[i];
            den1 += v1[i] * v1[i];
            den2 += v2[i] * v2[i];
        }
        rad = num / (sqrt(den1) * sqrt(den2));
    }
    else throw runtime_error("Нельзя найти угол, так как вектора имеют разную длину.");
    return (acos(rad) * 180.0 / PI);
}

//транспонирование матрицы
Matrix Transposition(const Matrix& A) {
    int N = A.Num_line();
    int M = A.Num_column();
    Matrix result(M, N);
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            result.MATRIX[i][j] = A.Element_ij(j, i);
        }
    }
    return result;
}

//ранг матрицы используя алгоритм Гаусса
int Rank(const Matrix& A) {
    int N = A.Num_line(), M = A.Num_column();
    int rank = N;
    Matrix a = A;
    int R2 = 0;
    for (int i = 0; i < N; i++) {
        int t = 0;
        for (int j = 0; j < M; j++) {
            if (a.Element_ij(i, j) == 0) t += 1;
            if (t == M) {
                ++R2;
            }
        }
    }
    if (R2 == rank) return 0;
    else {
        //вызов обработки строк
        int num_permutations = 0;
        for (int i = 0; i < min(N, M); i++) {
            Gauss_method(a, i, num_permutations);
        }
        int h = N;
        for (int i = 0; i < N; i++) {
            int t = 0;
            for (int j = 0; j < M; j++) {
                if (a.Element_ij(i, j) == 0) t += 1;
                if (t == M) {
                    --rank;
                }
            }
        }
        return rank;
    }
}

//для обнуления треугольника над диагональю
void Lower_Gauss(Matrix& a, int cur) {
    auto& M = a.GetMatrix();
    int ROW = a.Num_line(), COL = a.Num_column();
    for (int i = cur - 1; i >= 0; i--) {
        if (M[cur][cur] == 0)
            break;
        double coeff = M[i][cur] / M[cur][cur];
        for (int j = 0; j < COL; j++) {
            M[i][j] -= M[cur][j] * coeff;
            if (abs(M[i][j] - 0) < 0.00001)
                M[i][j] = 0;
        }
    }
}

//Обратная матрица
Matrix Reverse(const Matrix& A) {
    if (A.Num_line() != A.Num_column()) throw runtime_error("Нельзя найти обратную, так как матрица не квадратная.");
    if (Determinant(A) == 0) throw runtime_error("Нельзя найти обратную матрицу, так как вырожденная, ее определитель = 0");
    int N = A.Num_line();
    Matrix result(N, 2 * N);
    int t = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 2 * N; j++) {
            if (j < N) result.GetMatrix()[i][j] = A.Element_ij(i, j);
            else {
                if (t == i && j == t + N) {
                    result.GetMatrix()[i][j] = 1;
                    ++t;
                }
                else result.GetMatrix()[i][j] = 0;
            }
        }
    }
    int num_permutations = 0;
    for (int i = 0; i < N; i++) {
        Gauss_method(result, i, num_permutations);
    }
    for (int i = N - 1; i >= 0; i--)
        Lower_Gauss(result, i);
    for (int i = 0; i < N; ++i) {
        double del = 1;
        for (int j = 0; j < 2 * N; ++j) {
            if (i == j) del = result.Element_ij(i, j);
            result.GetMatrix()[i][j] /= del;
        }
    }
    Matrix res(N, N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            res.GetMatrix()[i][j] = result.GetMatrix()[i][j + N];
        }
    }
    return res;
}

//функция записи в текстовый файл
void Write_text_file(const string& file, const Matrix& M) {
    ofstream out(file);
    if (!out.is_open()) throw runtime_error("Файл не был открыт для записи.");
    else out << M;
    out.close();
}

//функция чтения из текстового файла
Matrix Read_text_file(const string& file) {
    ifstream in(file);
    if (!in.is_open()) throw runtime_error("Файл не был открыт для чтения или его не существует.");
    else {
        int ROW = 0, COL = 0;
        while (!in.eof()) {
            char c;
            in.get(c);
            if (c == ' ' || c == '\n') COL++;
            if (c == '\n') ROW++;
        }
        --ROW;
        COL = (COL - 1) / ROW;
        Matrix result(ROW, COL);
        in.close();
        ifstream inf(file, ios::app);
        inf >> result;
        return result;
    }
}

//запись в бинарный файл
void Write_bin_file(const string& file, const Matrix& M) {
    ofstream out(file, ios::binary);
    if (!out.is_open()) throw runtime_error("Файл не был открыт для записи.");
    else {
        int ROW = M.Num_line();
        int COL = M.Num_column();
        out.write((char*)&ROW, sizeof ROW);
        out.write((char*)&COL, sizeof COL);
        for (int i = 0; i < ROW; i++) {
            for (int j = 0; j < COL; j++) {
                out.write((char*)&M.Element_ij(i, j), sizeof(double));
            }
        }
    }
}

//чтение из бинарного файла
Matrix Read_bin_file(const string& file) {
    ifstream in(file, ios::binary);
    if (!in.is_open()) throw runtime_error("Файл не был открыт для чтения или его не существует.");
    else {
        int ROW = 0, COL = 0;
        in.read((char*)&ROW, sizeof(int));
        in.read((char*)&COL, sizeof(int));
        Matrix result(ROW, COL);
        for (int i = 0; i < ROW; i++) {
            for (int j = 0; j < COL; j++) {
                in.read((char*)&result.Element_ij(i, j), sizeof(double));
            }
        }
        return result;
    }

}

vector<double> aver(vector<vector<double>> M, int col, int line) {
    vector<double> Average(col);
    int k = 0;
    for (auto element : M) {
        for (auto el : element) {
            Average[k] += el;
            ++k;
        }
        k = 0;
    }
    for (int i = 0; i < Average.size(); ++i) {
        Average[i] /= line;
    }
    return Average;
}

//Центрирование матрицы
void PCA::Centering() {
    int line = MATRIX.Num_line();
    int col = MATRIX.Num_column();
    vector<double> Average = aver(MATRIX.GetMatrix(), col, line);
    for (auto& element : MATRIX.GetMatrix()) {
        int i = 0;
        for (auto& el : element) {
            el -= Average[i];
            ++i;
        }
    }
};

//Шкалирование матрицы
void PCA::Scaling() {
    //this->Centering();
    int line = MATRIX.Num_line();
    int col = MATRIX.Num_column();
    vector<double> norm(col);
    int k = 0;
    for (auto element : MATRIX.GetMatrix()) {
        for (auto el : element) {
            norm[k] += el * el;
            ++k;
        }
        k = 0;
    }
    for (int i = 0; i < norm.size(); ++i) {
        norm[i] /= (line - 1);
        norm[i] = sqrt(norm[i]);
    }
    /*for(int i = 0; i < norm.size(); ++i){
      cout << norm[i] << endl;
    }*/
    for (auto& element : MATRIX.GetMatrix()) {
        int i = 0;
        for (auto& el : element) {
            el /= norm[i];
            ++i;
        }
    }
};

//aлгоритм NIPALS
tuple<Matrix, Matrix, Matrix> PCA::NIPALS(size_t A, double epsilon) {
    //this->Centering();
    //this->Scaling();
    Matrix EE = MATRIX;
    int line = MATRIX.Num_line();
    int col = MATRIX.Num_column();
    Matrix t(line, 1), p, t_old;
    vector<Matrix> vec_t, vec_p;
    if (A == -10) A = col;

    for (int j = 0; j < A; ++j) {
        for (int i = 0; i < line; ++i) {
            t.Element_ij(i, 0) = MATRIX.Element_ij(i, j);
        }
        do {
            //p = Transposition(Transposition(t) * EE);
            p = Transposition(t) * EE;
            double t1_t = (Transposition(t) * t).Element_ij(0, 0);
            for (auto& element : p.GetMatrix()) {
                for (auto& el : element) el /= t1_t;
            }
            p = Transposition(p);

            double norm_p = Euclidean_norm(p);

            for (auto& element : p.GetMatrix()) {
                for (auto& el : element) el /= norm_p;
            }

            t_old = t;
            t = EE * p;
            double p1_p = (Transposition(p) * p).Element_ij(0, 0);

            for (auto& element : t.GetMatrix()) {
                for (auto& el : element) el /= p1_p;
            }
        } while (Euclidean_norm(t_old - t) > epsilon);
        EE = EE - t * Transposition(p);
        vec_p.push_back(p);
        vec_t.push_back(t);
    }

    Matrix T1(line, vec_t.size()), P1(vec_p[0].Num_line(), vec_p.size());
    for (int i = 0; i < line; ++i) {
        for (int r = 0; r < vec_t.size(); ++r) {
            T1.Element_ij(i, r) = vec_t[r].Element_ij(i, 0);
            if (i < vec_p[0].Num_line())
                P1.Element_ij(i, r) = vec_p[r].Element_ij(i, 0);
        }
    }
    T = T1; P = P1; E = EE;
    nip = true;
    tuple<Matrix, Matrix, Matrix> h = { T, P, E };
    return h;
};

/*vector<Matrix> PCA::NIPALS(size_t A, double epsilon) {
    //this->Centering();
    //this->Scaling();
    Matrix EE = MATRIX;
    int line = MATRIX.Num_line();
    int col = MATRIX.Num_column();
    Matrix t(line, 1), p, t_old;
    vector<Matrix> vec_t, vec_p;
    if (A == -10) A = col;

    for (int j = 0; j < A; ++j) {
        for (int i = 0; i < line; ++i) {
            t.Element_ij(i, 0) = MATRIX.Element_ij(i, j);
        }
        do {
            //p = Transposition(Transposition(t) * EE);
            p = Transposition(t) * EE;
            double t1_t = (Transposition(t) * t).Element_ij(0, 0);
            for (auto& element : p.GetMatrix()) {
                for (auto& el : element) el /= t1_t;
            }
            p = Transposition(p);

            double norm_p = Euclidean_norm(p);

            for (auto& element : p.GetMatrix()) {
                for (auto& el : element) el /= norm_p;
            }

            t_old = t;
            t = EE * p;
            double p1_p = (Transposition(p) * p).Element_ij(0, 0);

            for (auto& element : t.GetMatrix()) {
                for (auto& el : element) el /= p1_p;
            }
        } while (Euclidean_norm(t_old - t) > epsilon);
        EE = EE - t * Transposition(p);
        vec_p.push_back(p);
        vec_t.push_back(t);
    }

    Matrix T1(line, vec_t.size()), P1(vec_p[0].Num_line(), vec_p.size());
    for (int i = 0; i < line; ++i) {
        for (int r = 0; r < vec_t.size(); ++r) {
            T1.Element_ij(i, r) = vec_t[r].Element_ij(i, 0);
            if (i < vec_p[0].Num_line())
                P1.Element_ij(i, r) = vec_p[r].Element_ij(i, 0);
        }
    }
    T = T1; P = P1; E = EE;
    nip = true;
    vector<Matrix> h = { T, P, E };
    return h;
};*/

//вывод всех матриц
void PCA::PRINT() {
    cout << MATRIX << '\n';
    cout << T << '\n';
    cout << P << '\n';
    cout << E << '\n';
}
//размах
vector<double> PCA::Scope() const {
    if (!nip) throw runtime_error("Не был применен алгоритм NIPALS!");
    vector<double> scop;
    for (int i = 0; i < T.Num_line(); ++i) {
        Matrix t(1, T.Num_column());
        for (int j = 0; j < T.Num_column(); ++j) {
            t.Element_ij(0, j) = T.Element_ij(i, j);
        }
        auto T1_T = Reverse(Transposition(T) * T);
        auto h = (t * T1_T) * Transposition(t);
        scop.push_back(h.Element_ij(0, 0));
    }
    return scop;
}
//отклонения
vector<double> PCA::Deviations() const {
    if (!nip) throw runtime_error("Не был применен алгоритм NIPALS!");
    vector<double> devs;
    for (int i = 0; i < E.Num_line(); ++i) {
        double v = 0;
        for (int j = 0; j < E.Num_column(); ++j) {
            v += E.Element_ij(i, j) * E.Element_ij(i, j);
        }
        devs.push_back(v);
    }
    return devs;
}
//вычисление среднего (для всех образцов) расстояния
double PCA::Average_distance() {
    if (!nip) throw runtime_error("Не был применен алгоритм NIPALS!");
    vector<double> vi = Deviations();
    double v0 = 0;
    for (auto i : vi) v0 += i;
    return v0 / E.Num_line();
}

//дисперсия для полной и объясненной
vector<double> PCA::Dispersion() {
    if (!nip) throw runtime_error("Не был применен алгоритм NIPALS!");
    vector<double> deviations = Deviations();
    for (auto& i : deviations) {
        i /= E.Num_column();
    }
    return deviations;
}
//полная дисперсия
double PCA::Full_dispersion() {
    if (!nip) throw runtime_error("Не был применен алгоритм NIPALS!");
    vector<double> dispersion = Dispersion();
    double full = 0;
    for (auto i : dispersion) full += i;
    return full / E.Num_line();
}
//объясненная дисперсия
double PCA::Explained_dispersion() {
    if (!nip) throw runtime_error("Не был применен алгоритм NIPALS!");
    double Exp = 0;
    for (size_t i = 0; i < MATRIX.Num_line(); ++i) {
        for (size_t j = 0; j < MATRIX.Num_column(); ++j) {
            Exp += MATRIX.Element_ij(i, j) * MATRIX.Element_ij(i, j);
        }
    }
    return (1 - E.Num_line() * Average_distance() / Exp);
}