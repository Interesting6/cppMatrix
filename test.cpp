#include "iostream"
#include "matrix.h"
using namespace std;

int main(){
    Matrix A = Matrix(3, 3, 1);
    Matrix B = A;
    ifstream fin("inp.txt");
    fin >> A;
    A.Show();
    B.Show();
    Matrix C = Matrix(3, 1);
    fin >> C;
    C.Show();
    fin.close();
    Matrix D = Matrix(3, 1);
    Matrix::gaussianEliminate(A, C);
    C.Show();
}