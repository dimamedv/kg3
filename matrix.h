#ifndef MATRIX_H
#define MATRIX_H
#include "QMainWindow"

class Matrix
{
public:
    void rotateX(double angle);
    void rotateY(double angle);
    void rotateZ(double angle);
    static Matrix multiplyMatrixes(Matrix &m1, Matrix &m2);
    static Matrix transferAlong(Matrix &m1, Matrix &m2);

    Matrix(int nSize, int mSize);
    QVector<QVector<double>>data;
    int NSize, MSize;
    void multiplyOnValue(double value);
    void setRotateMatrix(double xAngle, double yAngle, double zAngle);
    void scale(Matrix scaleM);
};

#endif // MATRIX_H
