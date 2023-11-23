#include "matrix.h"
#include "cmath"

Matrix::Matrix(int nSize, int mSize)
{
    data = QVector<QVector<double>>(nSize, QVector<double>(mSize));
    NSize = nSize;
    MSize = mSize;
}

void Matrix::setRotateMatrix(double xAngle, double yAngle, double zAngle){
    data = {{1, 0, 0, 0},
            {0, 1, 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1}};
    NSize = 4;
    MSize = 4;

    rotateX(xAngle);
    rotateY(yAngle);
    rotateZ(zAngle);
}

void Matrix::rotateX(double angle){
    Matrix angleMatrix = Matrix(4,4);
    angleMatrix.data = {{1, 0, 0, 0},
                        {0, cos(angle), -sin(angle), 0},
                        {0, sin(angle), cos(angle), 0},
                        {0, 0, 0, 1}};

    data = multiplyMatrixes(angleMatrix, *this).data;
}

void Matrix::rotateY(double angle){
    Matrix angleMatrix = Matrix(4,4);
    angleMatrix.data = {{cos(angle), 0, -sin(angle), 0},
                        {0, 1, 0, 0},
                        {sin(angle), 0,  cos(angle), 0},
                        {0, 0, 0, 1}};

    data = multiplyMatrixes(angleMatrix, *this).data;
}

void Matrix::rotateZ(double angle){
    Matrix angleMatrix = Matrix(4,4);
    angleMatrix.data = {{cos(angle), sin(angle), 0, 0},
                        {-sin(angle), cos(angle), 0, 0},
                        {0, 0, 1, 0},
                        {0, 0, 0, 1}};

    data = multiplyMatrixes(angleMatrix, *this).data;
}

Matrix Matrix::multiplyMatrixes(Matrix &m1, Matrix &m2){
    int a = m1.NSize;
    int b = m2.NSize;
    int p = m2.MSize;

    Matrix result = Matrix(a, p);

    for (int i = 0; i < a; ++i)
        for (int j = 0; j < p ; ++j)
            for (int k = 0; k < b; ++k)
                result.data[i][j] += m1.data[i][k] * m2.data[k][j];

    return result;
}

void Matrix::multiplyOnValue(double value){
    Matrix scaleM = Matrix(4,4);
    scaleM.data = {{value, 0, 0, 0},
                   {0, value, 0, 0},
                   {0, 0, value, 0},
                   {0, 0, 0, 1}};

    scale(scaleM);
}

void Matrix::scale(Matrix scaleM){
    data = multiplyMatrixes(scaleM, *this).data;
}

Matrix Matrix::transferAlong(Matrix &m1, Matrix &m2){
    Matrix transM(4,4);
    transM.data = {{1, 0, 0, m2.data[0][0]},
                   {0, 1, 0, m2.data[1][0]},
                   {0, 0, 1, m2.data[2][0]},
                   {0, 0, 0, 1}};

    Matrix result = multiplyMatrixes(transM, m1);

    return result;
}
