#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <matrix.h>
#include <QPainter>
#include "QTimer"
#include "QComboBox"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public slots:
    void repaint();

public:

    struct Polygon;
    struct Point;

    MainWindow(QWidget *parent = nullptr);
    Ui::MainWindow *ui;
    QTimer* timer;
    ~MainWindow();
    void paintEvent(QPaintEvent *event);
    QVector<Polygon> polygons;
    QComboBox *projectionBox;

    bool intersectionShowMode = false;
    bool centroidIntersectionShowMode = false;
    bool centroidContainShowMode = false;

    QVector<double> cameraPos = {0,0,0};
    QVector<double> cameraAngles = {0,0,0};
    Matrix rotateMatrix = Matrix(4,4);
    double cameraScale = 20;
    double a;
    bool mousePressed = false;
    bool lab5Mod = false;
    QPointF lastMousePos;
    QVector<double> projectionCenter = {0,0,10};
    Matrix scaleAndTransferMatrix = Matrix(4,4);
    QPointF intPoint;

    QVector<QPointF> intersectionPoints;
    QVector<QPointF> intersectCentroidPoints;
    QVector<QPointF> containCentroidPoints;
    QVector<QPointF> doublerPoints;

    void setPolygons();

    QVector<MainWindow::Polygon> scaleAndTransferFigure();
    void drawPolygon(Polygon &poly, Matrix &currentRotateMatrix, QPainter &painter);
    QPolygonF scaleAndTransferPolygon(QPolygonF &poly);
    void draw(QVector<Polygon> &polygons, Matrix &currentRotateMatrix, QPainter &painter);
    void keyPressEvent(QKeyEvent *keyEvent);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void drawCenterProjectedPolygon(Polygon &poly, Matrix &currentRotateMatrix,  QPainter &painter);
    void drawCenterProjectedPolygons(QVector<Polygon> &polygons, Matrix &currentRotateMatrix, QPainter &painter);

    MainWindow::Point getCenterProjectPoint(MainWindow::Point &point);
    void drawFreeObliquePolygonProjection(Polygon &poly, Matrix &currentRotateMatrix,  double mod, QPainter &painter);
    MainWindow::Point getObliquePointProjection(MainWindow::Point &point, double mod);
    void drawFreeObliqueProjection(QVector<Polygon> &polygons, Matrix &currentRotateMatrix,  QPainter &painter);
    void drawCabinetObliqueProjection(QVector<Polygon> &polygons, Matrix &currentRotateMatrix,  QPainter &painter);
    void rotatePoint(Point &point, Matrix &angles);
    void sortByZ(QVector<Polygon> &polygons, Matrix &currentRotateMatrix,  QPainter &painter);
    bool areVectorsParallels(Point &first, Point &second);
    int arePolygonsIntersect(Polygon &p1, Polygon &p2, Polygon &realP1, Polygon &realP2);
    int arePolygonInPolygon(Polygon &p1, Polygon &p2, Polygon &realP1, Polygon &realP2);
    void drawSortedPolygons(QVector<Polygon> &polygons, Matrix &currentRotateMatrix,  QPainter &painter);
    double getZOfIntersection(Point &p1, Point &p2, Point &intersectionP);
    MainWindow::Point getMiddlePoint(Point &p1, Point &p2);
    void drawOrderedPolygon(Matrix &orderMatrix, QVector<bool> &drawedMatrix, int drawNomber, QVector<Polygon> &polygons, QPainter &painter);
    double getZOfIntersectionWithPolygon(Point &startPoint, Polygon &poly);
    void setScaleAndTransferMatrix(double scale, QPointF &transfer);
    Matrix getPointMatrix(Point &p);
    Matrix getQPointMatrix(QPointF &p);
    Matrix getCameraPosMatrix();

    bool arePointsCross(Point &firstP, Point &secondP, Point &thirdP, Point &fourthP);
    bool arePointsSame(Point &p1, Point &p2);
    bool arePointInPolygon(Point &point, Polygon &polgygon);
    bool isPolygonContainPoint(Point &point, Polygon &poly);
    MainWindow::Point getEquationOfLine(Point fPoint, Point sPoint);
    void move(double val);
private:

};
#endif // MAINWINDOW_H
