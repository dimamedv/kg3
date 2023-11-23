#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "matrix.h"
#include "QPainter"
#include "QKeyEvent"
#include "QTimer"
#include "QMouseEvent"
#include "QDebug"
#include "cmath"
#include "QComboBox"
#include "QFocusEvent"

void MainWindow::repaint()
{
    QWidget::repaint();
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    setFocusPolicy(Qt::StrongFocus);
    ui->setupUi(this);
    setPolygons();

    rotateMatrix.data = {{1,0,0, 0},
                         {0,1,0, 0},
                         {0,0,1, 0},
                         {0,0,0,1}};

    projectionBox = new QComboBox(this);
    projectionBox->addItem("центральная");
    projectionBox->addItem("косоугольная кабинетная");
    projectionBox->addItem("косоугольная свободная");
    projectionBox->addItem("параллельная");
    projectionBox->addItem("ортографическая ZOX");
    projectionBox->addItem("ортографическая ZOY");
    projectionBox->addItem("ортографическая XOY");
    projectionBox->hidePopup();

    timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(repaint()));
    timer->start(16);
}



struct MainWindow::Point{
public:
    Point(double x, double y, double z){
        X = x;
        Y = y;
        Z = z;
    }

    double getX(){
        return X;
    }

    double getY(){
        return Y;
    }

    double getZ(){
        return Z;
    }

    void muliply(double val){
        X *= val;
        Y *= val;
        Z *= val;
    }

    bool operator == (const Point &other){
        return X == other.X && Y == other.Y && Z == other.Z;
    }

    bool operator != (const Point &other){
        return X != other.X || Y != other.Y || Z != other.Z;
    }

    const Point operator + (const Point &p2){
        return Point(X + p2.X, Y + p2.Y, Z + p2.Z);
    }

    const Point operator - (const Point &p2){
        return Point(X - p2.X, Y - p2.Y, Z - p2.Z);
    }

private:
    double X, Y, Z;
};

struct MainWindow::Polygon{
    Polygon(){};
    int size(){
        return points.size();
    }

    void addPoint(double x, double y, double z){
        points.push_back(Point(x,y,z));
    }

    void addPoint(Point p){
        points.push_back(p);
    }

    Point *getPoint(size_t i){
        return &points[i];
    }

    void calculateCentroid(){
        double x = 0;
        double y = 0;
        double z = 0;
        foreach (Point p, points){
            x += p.getX();
            y += p.getY();
            z += p.getZ();
        }

        centroid = Point(x / points.size(), y / points.size(), z / points.size());
    }

    Point getCentroid(){
        return centroid;
    }

private:
    Point centroid = Point(0,0,0);
    QVector<Point> points;
};

Matrix MainWindow::getPointMatrix(Point &p){
    Matrix M(4,1);
    M.data = {{p.getX()},
              {p.getY()},
              {p.getZ()},
              {1}};

    return M;
}

Matrix MainWindow::getQPointMatrix(QPointF &p){
    Matrix M(4,1);

    M.data[0][0] = p.rx();
    M.data[1][0] = p.ry();
    M.data[3][0] = 1;

    return M;
}

Matrix MainWindow::getCameraPosMatrix(){
    Matrix CameraM(4,1);
    CameraM.data = {{cameraPos[0]},
                    {cameraPos[1]},
                    {cameraPos[2]},
                    {1}};

    return CameraM;
}

void MainWindow::setScaleAndTransferMatrix(double scale, QPointF &transfer){
    scaleAndTransferMatrix = Matrix(4,4);
    scaleAndTransferMatrix.data = {{scale, 0, 0, 0},
                                   {0, scale, 0, 0},
                                   {0, 0, scale, 0},
                                   {0, 0, 0, 1}};

    Matrix transM = Matrix(4,4);
    transM.data = {{1, 0, 0, transfer.rx()},
                   {0, 1, 0, transfer.ry()},
                   {0, 0, 1, 0},
                   {0, 0, 0, 1}};

    scaleAndTransferMatrix = Matrix::multiplyMatrixes(transM, scaleAndTransferMatrix);
}

QPolygonF MainWindow::scaleAndTransferPolygon(QPolygonF &poly){
    QPolygonF transferedPoly;

    for (int j = 0; j < poly.size();j++) {
        QPointF p = poly.at(j);
        Matrix point = getQPointMatrix(p);

        point = Matrix::multiplyMatrixes(scaleAndTransferMatrix, point);

        QPointF newP(point.data[0][0], point.data[1][0]);
        transferedPoly.push_back(newP);
    }

    return transferedPoly;
}

void MainWindow::drawPolygon(Polygon &poly, Matrix &currentRotateMatrix, QPainter &painter){
    QPolygonF newPoly;
    Matrix CameraM = getCameraPosMatrix();

    CameraM.multiplyOnValue(10);

    for (int i = 0 ; i < poly.size(); i++){
        Point point = *poly.getPoint(i);
        Matrix M = getPointMatrix(point);

        M = Matrix::multiplyMatrixes(currentRotateMatrix, M);

        point = Point(M.data[0][0], M.data[1][0], M.data[2][0]);
        newPoly.push_back(QPointF(point.getX(), point.getY()));
    }

    newPoly = scaleAndTransferPolygon(newPoly);

    painter.drawPolygon(newPoly);
}

void MainWindow::draw(QVector<Polygon> &polygons, Matrix &currentRotateMatrix , QPainter &painter){

    foreach(Polygon a, polygons){
        drawPolygon(a, currentRotateMatrix, painter);
    }
}

void MainWindow::setPolygons(){
                double boxSize = 2;
                int countOfSquares = 10;
                QVector<Point> sidePoints;
                Point p(- boxSize / 2, -boxSize / 2, -boxSize / 2);
                for (int i = 0; i < countOfSquares; i++) {
                    for (int j = 0; j < countOfSquares; j++) {
                        sidePoints.push_back(p);
                        p = Point(p.getX() + boxSize / countOfSquares, p.getY(), p.getZ());
                    }
                    p = Point(-boxSize / 2, p.getY() + boxSize / countOfSquares, p.getZ());
                }

                for (int i = 0; i < countOfSquares - 1; i++) {
                    for(int j = 0; j < countOfSquares - 1; j++){
                        Polygon newPoly;
                        newPoly.addPoint(sidePoints[i * countOfSquares + j]);
                        newPoly.addPoint(sidePoints[i * countOfSquares + j + 1]);
                        newPoly.addPoint(sidePoints[(i + 1) * countOfSquares + j + 1]);
                        newPoly.addPoint(sidePoints[(i + 1) * countOfSquares + j]);

                        polygons.push_back(newPoly);
                    }
                }

    QVector<Point> baseInnerFigure;

    int pointCount = 10;

    baseInnerFigure.push_back(Point(-4,0,0));
    baseInnerFigure.push_back(Point(-2,1,0));
    baseInnerFigure.push_back(Point(-1,2,0));
    baseInnerFigure.push_back(Point(1,2,0));
    baseInnerFigure.push_back(Point(2,1,0));
    baseInnerFigure.push_back(Point(4,0,0));
    baseInnerFigure.push_back(Point(2,-1,0));
    baseInnerFigure.push_back(Point(1,-2,0));
    baseInnerFigure.push_back(Point(-1,-2,0));
    baseInnerFigure.push_back(Point(-2,-1,0));

    if(pointCount > 10){

        pointCount -= 10;
        while(pointCount > 0){
            QVector<Point> newInnerFigure;
            newInnerFigure.push_back(baseInnerFigure[0]);
            for (int i = 1; i < baseInnerFigure.size(); i ++) {
                if(pointCount > 0){
                    Point middlePoint = Point((newInnerFigure.back().getX() + baseInnerFigure[i].getX()) / 2,
                                              (newInnerFigure.back().getY() + baseInnerFigure[i].getY()) / 2,
                                              (newInnerFigure.back().getZ() + baseInnerFigure[i].getZ()) / 2);
                    newInnerFigure.push_back(middlePoint);
                    pointCount --;
                }
                newInnerFigure.push_back(baseInnerFigure[i]);
            }

            if(pointCount > 0){
                Point middlePoint = Point((newInnerFigure.back().getX() + baseInnerFigure[0].getX()) / 2,
                        (newInnerFigure.back().getY() + baseInnerFigure[0].getY()) / 2,
                        (newInnerFigure.back().getZ() + baseInnerFigure[0].getZ()) / 2);
                newInnerFigure.push_back(middlePoint);
                pointCount --;
            }

            baseInnerFigure = newInnerFigure;
        }
    }
    Point apperPoint(0,0,4);
    Point lowerPoint(0,0,-4);

    for (int i = 0; i < baseInnerFigure.size() - 1; i++) {
        Polygon poly;

        poly.addPoint(baseInnerFigure[i]);
        poly.addPoint(baseInnerFigure[i + 1]);
        poly.addPoint(apperPoint);
        polygons.push_back(poly);

        poly = Polygon();

        poly.addPoint(baseInnerFigure[i]);
        poly.addPoint(baseInnerFigure[i + 1]);
        poly.addPoint(lowerPoint);
        polygons.push_back(poly);
    }

    Polygon poly;

    poly.addPoint(baseInnerFigure[0]);
    poly.addPoint(baseInnerFigure[baseInnerFigure.size() - 1]);
    poly.addPoint(apperPoint);
    polygons.push_back(poly);

    poly = Polygon();

    poly.addPoint(baseInnerFigure[0]);
    poly.addPoint(baseInnerFigure[baseInnerFigure.size() - 1]);
    poly.addPoint(lowerPoint);
    polygons.push_back(poly);


}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::mousePressEvent(QMouseEvent *event){
    if(event->button() == Qt::MiddleButton){
        qDebug() << "Mouse pressed";
        QPointF globalCursorPos = QCursor::pos();
        if(mousePressed){
            double deltaX = globalCursorPos.rx() - lastMousePos.rx();
            double deltaY = globalCursorPos.ry() - lastMousePos.ry();
            rotateMatrix.rotateX(deltaX);
            rotateMatrix.rotateY(deltaY);

            lastMousePos = globalCursorPos;
        }
        else{
            mousePressed = true;
            lastMousePos = globalCursorPos;
        }
    }
}



void MainWindow::mouseMoveEvent(QMouseEvent *event){
    QPointF globalCursorPos = QCursor::pos();
    if(mousePressed){
        double deltaX = 0.01 * (globalCursorPos.rx() - lastMousePos.rx());
        double deltaY = 0.01 *(globalCursorPos.ry() - lastMousePos.ry());
        rotateMatrix.rotateX(-deltaY);
        rotateMatrix.rotateY(-deltaX);

        lastMousePos = globalCursorPos;
    }
}

void MainWindow::mouseReleaseEvent(QMouseEvent *event){
    if(event->button() == Qt::MiddleButton){
        mousePressed = false;
    }
}

void MainWindow::move(double val){
//    foreach(Polygon poly,  polygons){
//        foreach(Point p, poly)
//            p = Point (p.getX(), p.getY(), p.getZ() + val);
//    }
}

void MainWindow::keyPressEvent(QKeyEvent *keyEvent){
    qDebug() << "mod changed";
    int key = keyEvent->key();
    if(key == Qt::Key_Plus)
        cameraScale += 0.4;
    else if(key == Qt::Key_Minus)
        cameraScale -= 0.4;

    else if(key == Qt::Key_Q)
        rotateMatrix.rotateX(-0.1);
    else if(key == Qt::Key_W)
        rotateMatrix.rotateX(0.1);
    else if(key == Qt::Key_A)
        rotateMatrix.rotateY(-0.1);
    else if(key == Qt::Key_S)
        rotateMatrix.rotateY(0.1);
    else if(key == Qt::Key_Z)
        rotateMatrix.rotateZ(-0.1);
    else if(key == Qt::Key_X)
        rotateMatrix.rotateZ(0.1);

    else if(key == Qt::Key_P)
        lab5Mod = !lab5Mod;

    else if(key == Qt::Key_Up)
        projectionCenter[2] += 1;
    else if(key == Qt::Key_Down)
        projectionCenter[2] -= 1;

    else if(key == Qt::Key_K)
        intersectionShowMode = !intersectionShowMode;
    else if(key == Qt::Key_J)
        centroidIntersectionShowMode = !centroidIntersectionShowMode;
    else if(key == Qt::Key_L)
        centroidContainShowMode = !centroidContainShowMode;

//    else if(key == Qt::Key_M)
//        move(0.5);
//    else if(key == Qt::Ke)
}

void MainWindow::rotatePoint(Point &point, Matrix &angles){
    Matrix M = getPointMatrix(point);

    M = Matrix::multiplyMatrixes(angles, M);

    point = Point(M.data[0][0], M.data[1][0], M.data[2][0]);
}


MainWindow::Point MainWindow::getCenterProjectPoint(MainWindow::Point &point){
    Point projectPoint = Point(projectionCenter[0],projectionCenter[1],projectionCenter[2]);

    Point projectedPoint = point;
    Matrix pointSection(4, 1);
    pointSection.data = {{ projectedPoint.getX() - projectPoint.getX()},
                         { projectedPoint.getY() - projectPoint.getY()},
                         { projectedPoint.getZ() - projectPoint.getZ()},
                         {1}};

    double mod = (projectedPoint.getZ() - 4) / pointSection.data[2][0];
    pointSection.multiplyOnValue(mod);

    Matrix projectedPointMatrix = getPointMatrix(point);

    projectedPointMatrix = Matrix::transferAlong(projectedPointMatrix, pointSection);
    projectedPoint = Point(projectedPointMatrix.data[0][0], projectedPointMatrix.data[1][0], 0);

    return projectedPoint;
}

void MainWindow::sortByZ(QVector<Polygon> &polygons, Matrix &currentRotateMatrix, QPainter &painter){
    painter.setBrush(QColor(255, 0, 0, 255));
    QVector<Polygon> transferedPolygons;

    foreach(Polygon poly, polygons){
        Polygon newPoly;

        for (int i = 0 ; i < poly.size(); i++){
            Point point = *poly.getPoint(i);

            rotatePoint(point, currentRotateMatrix);

            newPoly.addPoint(point);
        }

        newPoly.calculateCentroid();
        transferedPolygons.push_back(newPoly);
    }

    for (int i = 0; i < transferedPolygons.size(); i++) {
        Polygon currentPoly = transferedPolygons.at(i);
        int j = i - 1;
        while( j >= 0 && currentPoly.getCentroid().getZ() < transferedPolygons[j].getCentroid().getZ()){
            transferedPolygons[j + 1] = transferedPolygons[j];
            j --;
        }

        transferedPolygons[j + 1] = currentPoly;
    }

    foreach(Polygon poly, transferedPolygons){
        QPolygonF newPoly;
        for (int i = 0; i < poly.size(); i++){
            Point newPoint = getCenterProjectPoint(*poly.getPoint(i));
            newPoly.push_back(QPointF(newPoint.getX(), newPoint.getY()));
        }
        painter.drawPolygon(scaleAndTransferPolygon(newPoly));
    }
}

void MainWindow::drawCenterProjectedPolygon(Polygon &poly, Matrix &currentRotateMatrix, QPainter &painter){
    QPolygonF newPoly;
    Matrix CameraM = getCameraPosMatrix();

    CameraM.multiplyOnValue(10);

    for (int i = 0 ; i < poly.size(); i++){
        Point point = *poly.getPoint(i);

        rotatePoint(point, currentRotateMatrix);

        point = getCenterProjectPoint(point);

        newPoly.push_back(QPointF(point.getX(), point.getY()));
    }

    painter.drawPolygon(scaleAndTransferPolygon(newPoly));
}

void MainWindow::drawCenterProjectedPolygons(QVector<Polygon> &polygons, Matrix &currentRotateMatrix, QPainter &painter){
    foreach(Polygon poly, polygons)
        drawCenterProjectedPolygon(poly, currentRotateMatrix, painter);
}

MainWindow::Point MainWindow::getObliquePointProjection(MainWindow::Point &point, double mod){
    Point projectedPoint = point;
    Matrix projectedPointMatrix = getPointMatrix(projectedPoint);

    Matrix projectMatrix(4,4);
    projectMatrix.data = {{1, 0, mod * cos(0.79), 0},
                          {0, 1, mod * cos(0.79), 0},
                          {0, 0, 1, 0},
                          {0, 0, 0, 1}};

    projectedPointMatrix = Matrix::multiplyMatrixes( projectMatrix, projectedPointMatrix);

    projectedPoint = Point(projectedPointMatrix.data[0][0], projectedPointMatrix.data[1][0], projectedPointMatrix.data[2][0]);
    return projectedPoint;
}

void MainWindow::drawFreeObliquePolygonProjection(Polygon &poly, Matrix &currentRotateMatrix, double mod , QPainter &painter){
    QPolygonF newPoly;
    Matrix CameraM = getCameraPosMatrix();
    CameraM.multiplyOnValue(10);


    for (int i = 0 ; i < poly.size(); i++){
        Point point = *poly.getPoint(i);

        rotatePoint(point, currentRotateMatrix);

        point = getObliquePointProjection(point, mod);

        newPoly.push_back(QPointF(point.getX(), point.getY()));
    }
    painter.drawPolygon(scaleAndTransferPolygon(newPoly));
}


void MainWindow::drawFreeObliqueProjection(QVector<Polygon> &polygons, Matrix &currentRotateMatrix, QPainter &painter){
    foreach(Polygon poly, polygons)
        drawFreeObliquePolygonProjection(poly, currentRotateMatrix, 1, painter);
}


void MainWindow::drawCabinetObliqueProjection(QVector<Polygon> &polygons, Matrix &currentRotateMatrix, QPainter &painter){
    foreach(Polygon poly, polygons)
        drawFreeObliquePolygonProjection(poly, currentRotateMatrix, 0.5, painter);
}

bool MainWindow::areVectorsParallels(Point &first, Point &second){
    double k1, k2;
    if(first.getX() == 0 && first.getY() != 0)
        k1 = INFINITY;
    else if(first.getY() == 0)
        k1 = 0;
    else
        k1 = first.getY() / first.getX();

    if(second.getX() == 0 && second.getY() != 0)
        k2 = INFINITY;
    else if(second.getY() == 0)
        k2 = 0;
    else
        k2 = second.getY() / second.getX();

    return k1 == k2;
}

double MainWindow::getZOfIntersection(Point &p1, Point &p2, Point &intersectionP){
    double res;
    if(abs(p2.getX() - p1.getX()) > abs(p2.getY() - p1.getY()))
        res = (intersectionP.getX() - p1.getX()) * (p2.getZ() - p1.getZ()) / (p2.getX() - p1.getX()) +p1.getZ();
    else
        res = (intersectionP.getY() - p1.getY()) * (p2.getZ() - p1.getZ()) / (p2.getY() - p1.getY()) +p1.getZ();

    return res;
}

MainWindow::Point MainWindow::getMiddlePoint(Point &p1, Point &p2){
    return Point((p1.getX() + p2.getX())/2, (p1.getY() + p2.getY())/2, (p1.getZ() + p2.getZ()) / 2);
}

double MainWindow::getZOfIntersectionWithPolygon(Point &startPoint, Polygon &poly){
    Point A = *poly.getPoint(1);
    Point AB = *poly.getPoint(0) - *poly.getPoint(1);
    Point AC = *poly.getPoint(2) - *poly.getPoint(1);

    double nXMod = AB.getY() * AC.getZ() - AB.getZ() * AC.getY();
    double nYMod = AB.getZ() * AC.getX() - AB.getX() * AC.getZ();
    double nZMod = AB.getX() * AC.getY() - AB.getY() * AC.getX();

    double D = -A.getX() * nXMod - A.getY() * nYMod - A.getZ() * nZMod;

    return (-D - nXMod * startPoint.getX() - nYMod * startPoint.getY()) / nZMod;
}

MainWindow::Point MainWindow::getEquationOfLine(Point fPoint, Point sPoint){
    return Point(sPoint.getY() - fPoint.getY(), fPoint.getX() - sPoint.getX(), fPoint.getY() * (sPoint.getX() - fPoint.getX()) - fPoint.getX() * (sPoint.getY() - fPoint.getY()));
}

bool MainWindow::arePointsCross(Point &firstP, Point &secondP, Point &thirdP, Point &fourthP) {
    double k1 = (secondP.getY() - firstP.getY()) / (secondP.getX() - firstP.getX());
    double k2 = (fourthP.getY() - thirdP.getY()) / (fourthP.getX() - thirdP.getX());
    if(abs(k1 - k2) < 0.01)
        return false;

    Point fEquation = getEquationOfLine(firstP, secondP);
    Point sEquation = getEquationOfLine(thirdP, fourthP);

    Point pseudoX = fEquation;
    pseudoX.muliply(-1 / fEquation.getX());
    Point pseudoY = Point(0, sEquation.getY() + pseudoX.getY() * sEquation.getX(), sEquation.getZ() + pseudoX.getZ() * sEquation.getX());

    double y = -pseudoY.getZ() / pseudoY.getY();

    double x = pseudoX.getZ() + pseudoX.getY() * y;

    intPoint = QPointF(x, y);
    return 1;
}

bool MainWindow::arePointsSame(Point &p1, Point &p2){
    double delta = abs(p1.getX() - p2.getX());
    delta = abs(p1.getY() - p2.getY()) > delta ? abs(p1.getY() - p2.getY()): delta;
    delta = abs(p1.getZ() - p2.getZ()) > delta ? abs(p1.getZ() - p2.getZ()): delta;

    return delta < 0.1;
}

int MainWindow::arePolygonsIntersect(Polygon &p1, Polygon &p2, Polygon &realP1, Polygon &realP2){
    int innerOrder = 0;
    for (int i = 0; i <  p1.size(); i++) {
        Point firstP = *p1.getPoint(i);
        Point secondP = *p1.getPoint((i + 1) % p1.size());
        if(secondP.getX() < firstP.getX())
            qSwap(secondP, firstP);

        for (int j = 0; j < p2.size(); j ++) {
            Point thirdP = *p2.getPoint(j);
            Point fourthP = *p2.getPoint((j + 1) % p2.size());
            if(fourthP.getX() < thirdP.getX())
                qSwap(fourthP, thirdP);

            int countOfSamePoints = 0;
            countOfSamePoints += arePointsSame(firstP, thirdP);
            countOfSamePoints += arePointsSame(firstP, fourthP);
            countOfSamePoints += arePointsSame(secondP, thirdP);
            countOfSamePoints += arePointsSame(secondP, fourthP);

            if(countOfSamePoints == 0){
                if(arePointsCross(firstP, secondP, thirdP, fourthP))
                    if(intPoint.rx() > firstP.getX() && intPoint.rx() < secondP.getX() && intPoint.rx() > thirdP.getX() && intPoint.rx() < fourthP.getX()){
                        Point intersectionPoint = Point(intPoint.rx(),intPoint.ry(), 0);
                        intersectionPoints.push_back(QPointF(intersectionPoint.getX(), intersectionPoint.getY()));
                        return getZOfIntersection(*realP1.getPoint(i), *realP1.getPoint((i + 1) % p1.size()), intersectionPoint) > getZOfIntersection(*realP2.getPoint(j), *realP2.getPoint((j + 1) % p2.size()), intersectionPoint) ? 1 : -1;
                    }
            }

            else if(innerOrder == 0 && countOfSamePoints >= 2){
                Point p = realP1.getCentroid();
                if(arePointInPolygon(p, realP2)){
                    intersectCentroidPoints.push_back(QPointF(p1.getCentroid().getX(), p1.getCentroid().getY()));
                    innerOrder =  p.getZ() > getZOfIntersectionWithPolygon(p, realP2) ? 1 : -1;
                }
                else{
                    p = realP2.getCentroid();
                    if(arePointInPolygon(p, realP1)){
                        intersectCentroidPoints.push_back(QPointF(p2.getCentroid().getX(), p2.getCentroid().getY()));
                        innerOrder =  p.getZ() < getZOfIntersectionWithPolygon(p, realP1) ? 1 : -1;
                    }
                }
            }
        }
    }

    return innerOrder;
}


bool MainWindow::arePointInPolygon(Point &point, Polygon &polgygon){
    bool result = false;
    int j = polgygon.size() - 1;
    for (int i = 0; i < polgygon.size(); i++) {
        if ( ((polgygon.getPoint(i)->getY() < point.getY() && polgygon.getPoint(j)->getY() >= point.getY()) || (polgygon.getPoint(j)->getY() < point.getY() && polgygon.getPoint(i)->getY() >= point.getY())) &&
             (polgygon.getPoint(i)->getX() + (point.getY() - polgygon.getPoint(i)->getY()) / (polgygon.getPoint(j)->getY() - polgygon.getPoint(i)->getY()) * (polgygon.getPoint(j)->getX() - polgygon.getPoint(i)->getX()) < point.getX()))
            result = !result;
        j = i;
    }

    return result;

}

bool MainWindow::isPolygonContainPoint( Point &point, Polygon& poly){
    for (int i = 0; i < poly.size(); i++) {
        if(arePointsSame(*poly.getPoint(i),point))
            return true;
    }

    return false;
}

int MainWindow::arePolygonInPolygon(Polygon &p1, Polygon &p2, Polygon &realP1, Polygon &realP2){
    bool result = true;
    for (int k = 0; k < p1.size() && result; k++) {
        result = arePointInPolygon(*p1.getPoint(k), p2) || isPolygonContainPoint(*p1.getPoint(k), p2);
    }
    if(result){
        Point p = realP1.getCentroid();
        containCentroidPoints.push_back(QPointF(p1.getCentroid().getX(), p1.getCentroid().getY()));
        return p.getZ() > getZOfIntersectionWithPolygon(p, realP2) ? 1 : -1;
    }

    else {
        result = true;
        for (int k = 0; k < p2.size() && result; k++) {
            result = arePointInPolygon(*p2.getPoint(k), p1) || isPolygonContainPoint(*p2.getPoint(k), p1);
        }
        if(result){
            Point p = realP2.getCentroid();
            containCentroidPoints.push_back(QPointF(p2.getCentroid().getX(), p2.getCentroid().getY()));
            return p.getZ() < getZOfIntersectionWithPolygon(p, realP1) ? 1 : -1;
        }
    }
}

void MainWindow::drawOrderedPolygon(Matrix &orderMatrix, QVector<bool> &drawedMatrix, int drawNomber, QVector<Polygon> &polygons, QPainter &painter){
    if(!drawedMatrix[drawNomber]){
        drawedMatrix[drawNomber] = true;
        for (int i = 0; i < drawedMatrix.size(); i++) {
            if(!drawedMatrix[i] && orderMatrix.data[drawNomber][i] == 1)
                drawOrderedPolygon(orderMatrix, drawedMatrix, i, polygons, painter);
        }

        QPolygonF newPoly;
        for (int i = 0; i < polygons[drawNomber].size(); i++){
            QPointF newPoint(polygons[drawNomber].getPoint(i)->getX(), polygons[drawNomber].getPoint(i)->getY());
            newPoly.push_back(newPoint);
        }
        painter.drawPolygon(scaleAndTransferPolygon(newPoly));
    }
}

void MainWindow::drawSortedPolygons(QVector<Polygon> &polygons, Matrix &currentRotateMatrix, QPainter &painter){
    intersectionPoints.clear();
    intersectCentroidPoints.clear();
    containCentroidPoints.clear();
    doublerPoints.clear();
    Matrix drawOrder(polygons.size(), polygons.size());

    QVector<Polygon> transferedPolygons, projectedPolygons;
    foreach(Polygon poly, polygons){
        Polygon transferedPoly, projectedPoly;

        for (int i = 0 ; i < poly.size(); i++){
            Point point = *poly.getPoint(i);

            rotatePoint(point, currentRotateMatrix);

            transferedPoly.addPoint(point);
            projectedPoly.addPoint(getCenterProjectPoint(point));
        }

        transferedPoly.calculateCentroid();
        transferedPolygons.push_back(transferedPoly);
        projectedPoly.calculateCentroid();
        projectedPolygons.push_back(projectedPoly);
    }


    for (int i = 0; i < polygons.size(); i++) {
        for (int j = i + 1; j < polygons.size(); j++) {
            int order = arePolygonsIntersect(projectedPolygons[i], projectedPolygons[j], transferedPolygons[i], transferedPolygons[j]);
            if(order != 0){
                drawOrder.data[i][j] = order;
                drawOrder.data[j][i] = -order;
            }
            else{
                order = arePolygonInPolygon(projectedPolygons[i], projectedPolygons[j], transferedPolygons[i], transferedPolygons[j]);
                drawOrder.data[i][j] = order;
                drawOrder.data[j][i] = -order;
            }
        }
    }

    QVector<bool> drawedMatrix(projectedPolygons.size(), 0);

    for (int i = 0; i < projectedPolygons.size(); i++) {
        drawOrderedPolygon(drawOrder, drawedMatrix, i, projectedPolygons, painter);
    }

    if(intersectionShowMode)
        for (int i = 0; i < intersectionPoints.size(); i++) {
            Matrix pointMatrix = getQPointMatrix(intersectionPoints[i]);
            pointMatrix = Matrix::multiplyMatrixes(scaleAndTransferMatrix, pointMatrix);

            intersectionPoints[i] = QPointF(pointMatrix.data[0][0], pointMatrix.data[1][0]);
            if(i > 0 && intersectionPoints[i - 1] == intersectionPoints[i])
                doublerPoints.push_back(intersectionPoints[i]);

            painter.setBrush(QColor(0, 255, 0, 255));
            painter.drawEllipse(intersectionPoints[i].rx() - 2.5, intersectionPoints[i].ry() - 2.5, 5, 5);
        }
    if(centroidIntersectionShowMode)
        for (int i = 0; i < intersectCentroidPoints.size(); i++) {
            Matrix pointMatrix = getQPointMatrix(intersectCentroidPoints[i]);
            pointMatrix = Matrix::multiplyMatrixes(scaleAndTransferMatrix, pointMatrix);
            intersectCentroidPoints[i] = QPointF(pointMatrix.data[0][0], pointMatrix.data[1][0]);

            if(i > 0 && intersectCentroidPoints[i - 1] == intersectCentroidPoints[i])
                doublerPoints.push_back(intersectCentroidPoints[i]);

            painter.setBrush(QColor(0, 0, 255, 255));
            painter.drawEllipse(intersectCentroidPoints[i].rx() - 2.5, intersectCentroidPoints[i].ry() - 2.5, 5, 5);
        }

    if(centroidContainShowMode)
        for (int i = 0; i < containCentroidPoints.size(); i++) {
            Matrix pointMatrix = getQPointMatrix(containCentroidPoints[i]);
            pointMatrix = Matrix::multiplyMatrixes(scaleAndTransferMatrix, pointMatrix);
            containCentroidPoints[i] = QPointF(pointMatrix.data[0][0], pointMatrix.data[1][0]);

            if(i > 0 && containCentroidPoints[i - 1] == containCentroidPoints[i])
                doublerPoints.push_back(containCentroidPoints[i]);

            painter.setBrush(QColor(0, 150, 150, 255));
            painter.drawEllipse(containCentroidPoints[i].rx() - 2.5, containCentroidPoints[i].ry() - 2.5, 5, 5);
        }

    for (int i = 0; i < doublerPoints.size(); i++) {
        painter.setBrush(QColor(150, 0, 150, 255));
        painter.drawEllipse(doublerPoints[i].rx() - 3.5, doublerPoints[i].ry() - 3.5, 7, 7);
    }
}

void MainWindow::paintEvent(QPaintEvent* event)
{
    QPainter painter(this);
    if (width() > height())
        a = (height() - 30) / 2;
    else a = (width() - 30) / 2;

    if (width() < 30 || height() < 30)
        return;

    painter.setPen(QPen(Qt::black, 2, Qt::SolidLine, Qt::FlatCap));
    if(lab5Mod)
        painter.setBrush(QColor(255, 0, 0, 255));
    else
        painter.setBrush(QColor(0,0,0,0));

    QPointF pos;

    Matrix curentRotateMatrix(4,4);
    if(!lab5Mod){
        pos = QPointF(width() * 3 / 4, height() * 3 / 4);
        setScaleAndTransferMatrix(cameraScale, pos);
        switch (projectionBox->currentIndex()) {
        case(0):
            drawCenterProjectedPolygons(polygons, rotateMatrix,  painter);
            break;
        case(1):
            drawCabinetObliqueProjection(polygons, rotateMatrix,  painter);
            break;
        case(2):
            drawFreeObliqueProjection(polygons, rotateMatrix,  painter);
            break;
        case(3):
            draw(polygons, rotateMatrix,  painter);
            break;
        case(4):
            curentRotateMatrix.setRotateMatrix(1.54, 0, 0);
            draw(polygons, curentRotateMatrix,  painter);
            break;
        case(5):
            curentRotateMatrix.setRotateMatrix(0, 1.54, 0);
            draw(polygons, curentRotateMatrix,  painter);
            break;
        case(6):
            curentRotateMatrix.setRotateMatrix(0, 0, 0);
            draw(polygons, curentRotateMatrix,  painter);
            break;
        }


        QVector<double> baseCameraPos = cameraPos;
        QVector<double> baseCameraAngles = cameraAngles;

        pos = QPointF(width() / 4, height() / 4);
        curentRotateMatrix.setRotateMatrix(1.54, 0, 0);
        setScaleAndTransferMatrix(20, pos);
        draw(polygons, curentRotateMatrix, painter);


        pos = QPointF(width() * 3 / 4, height()/ 4);
        curentRotateMatrix.setRotateMatrix(0,1.54, 0);
        setScaleAndTransferMatrix(20, pos);
        draw(polygons, curentRotateMatrix, painter);

        pos = QPointF(width() / 4, height() * 3 / 4);
        curentRotateMatrix.setRotateMatrix( 0, 0, 0);
        setScaleAndTransferMatrix(20, pos);
        draw(polygons, curentRotateMatrix, painter);
    }

    else{
        projectionBox->hide();
        pos = QPointF(width() / 2, height() / 2);
        setScaleAndTransferMatrix(cameraScale, pos);
        pos = QPointF(width() / 2, height() / 2);
//            sortByZ(polygons, rotateMatrix,  painter);
        drawSortedPolygons(polygons, rotateMatrix, painter);
    }
}
