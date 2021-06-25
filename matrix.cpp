#include "matrix.h"

#include "pch.h"

using std::cout;
using std::endl;
using std::istream;

const double EPS = 1E-10;

// 分配空间
void Matrix::initialize() {
    p = new double*[row_num];  // 分配row_num个指针
    for (int i = 0; i < row_num; ++i) {
        p[i] = new double[col_num];  // 为p[i]分配col_num个空间
    }
}

// 构造函数，采用0初始化
Matrix::Matrix(int r, int c): row_num(r), col_num(c) {
    initialize();
    for (int i = 0; i < r; i++){
        for (int j = 0; j < c; j++){
            p[i][j] = 0;
        }
    }  
}

// 声明一个值全为1的矩阵
Matrix::Matrix(int r, int c, double value): row_num(r), col_num(c) {
    initialize();
    for (int i=0; i<r; ++i){
        for (int j=0; j<c; ++j){
            p[i][j] = value;
        }
    }
}

Matrix::~Matrix(){
    for (int i=0; i<row_num; ++i){
        delete[] p[i];
    }
    delete[] p;
    p = NULL; // 不然会报问题pointer being freed was not allocated
}

// 矩阵赋值，进行拷贝，不共享内存
Matrix& Matrix::operator=(const Matrix& m){
    if (this == &m){ // 当前对象地址 与 m地址相同
        return *this; // 它们指向同一块内存，返回自己就行
    }
    if (row_num != m.row_num || col_num != m.col_num){
        // 先释放之前的空间
        for (int i=0; i<row_num; ++i){
            delete[] p[i];
        }
        delete[] p;
        // 重新设置形状，并重新分配空间
        row_num = m.row_num;
        col_num = m.col_num;
        initialize();
    }
    for (int i=0; i< row_num; ++i){
        for (int j=0; j<col_num; ++j){
            p[i][j] = m.p[i][j];
        }
    }
    return *this;
}

// 矩阵赋值，不进行拷贝，共享内存，（要求矩阵的大小被声明过了）
Matrix& Matrix::operator=(double* m){
    for (int i=0; i< row_num; ++i){
        for (int j=0; j< col_num; ++j){
            p[i][j] = *(m + i*row_num + j);
        }
    }
    return *this;
}

Matrix& Matrix::operator+=(const Matrix& m1){
    // TODO 补充assert形状相同
    for (int i=0; i<row_num; ++i){
        for (int j=0; j<col_num; ++j){
            p[i][j] += m1.p[i][j];
        }
    }
    return *this;
}
Matrix& Matrix::operator-=(const Matrix& m1){
    // TODO 补充assert形状相同
    for (int i=0; i<row_num; ++i){
        for (int j=0; j<col_num; ++j){
            p[i][j] -= m1.p[i][j];
        }
    }
    return *this;
}
Matrix& Matrix::operator*=(const Matrix& m1){
    Matrix temp = Matrix(row_num, m1.col_num);
    for (int i=0; i<row_num; i++){
        for (int j=0; j<m1.col_num; j++){
            for (int k=0; k<col_num; k++){
                temp.p[i][j] += (p[i][k] * m1.p[k][j]);
            }
        }
    }
    *this = temp;
    return *this;
}
Matrix Matrix::operator*(const Matrix& m1) const{
    Matrix temp = Matrix(row_num, m1.col_num);
    for (int i=0; i<row_num; i++){
        for (int j=0; j<m1.col_num; j++){
            for (int k=0; k<col_num; k++){
                temp.p[i][j] += (p[i][k] * m1.p[k][j]);
            }
        }
    }
    return temp;
}

Matrix Matrix::Hadamard(const Matrix& m1, const Matrix& m2){
    Matrix temp = Matrix(m1.row_num, m1.col_num);
    for (int i=0; i<m1.row_num; i++){
        for (int j=0; j<m1.col_num; j++){
            temp.p[i][j] = (m1.p[i][j] * m2.p[i][j]);
        }
    }
    return temp;
}

int Matrix::row() const{
    return row_num;
}
int Matrix::col() const{
    return col_num;
}
double Matrix::Point(int i, int j) const{
    return this->p[i][j];
}
void Matrix::Show() const {
    cout << "row=" << row_num <<", col="<<col_num<< endl;//显示矩阵的行数和列数
	for (int i = 0; i < row_num; i++) {
		for (int j = 0; j < col_num; j++) {
			cout << p[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

}
Matrix Matrix::eye(int n){
    Matrix E = Matrix(n, n);
    for (int i=0; i<n; i++){
        E.p[i][i] = 1;
    }
    return E;
}
Matrix Matrix::T(const Matrix& m){
    Matrix mT(m.col(), m.row());
    for (int i=0; i<m.row_num; i++){
        for (int j=0; j<m.col_num; j++){
            mT.p[j][i] = m.p[i][j];
        }
    }
    return mT;
}

std::ifstream& operator>>(std::ifstream& is, Matrix& A){
    int n;
    for (int i=0; i<A.row(); i++){
        for (int j=0; j<A.col(); j++){
            is >> n;
            A.p[i][j] = n;
        }
    }
    return is;
}

// Matrix Matrix::Solve(const Matrix& A, const Matrix& b){

// }

void Matrix::gaussianEliminate(const Matrix& A, const Matrix& B){
    int ar = A.row(), ac = A.col();
    int br = B.row(), bc = B.col();
    // TODO assert ar == br
    for (int i=0; i<ar; i++){
        int t = A.p[i][i];
        for (int j=i; j<ac; j++){
            A.p[i][j] /= t;
        }
        for (int j=0; j<bc; j++){
            B.p[i][j] /= t;
        }
        for (int k=i+1; k<ar; k++){ // i后面的每行，消去i列的位置
            for (int j=i; j<ac; j++){ // r_k - (r_k[i] * r_i)
                A.p[k][j] -= A.p[k][i]*A.p[i][j]; 
            }
            for (int j=0; j<bc; j++){
                B.p[i][j] -= A.p[k][i]*B.p[i][j];
            }
        }
    }
    A.Show();
    B.Show();
    // return B;
}

// Matrix Matrix::inv(Matrix& A){
//     int n = A.row();
//     Matrix E = eye(n);
//     for (int i=0; i< n; i++){

//     }

// }
