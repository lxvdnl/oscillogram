#include "mainwindow.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

#include "ui_mainwindow.h"

typedef double (*funcF)(double);
typedef double (*funcG)(double);
typedef double (*funcSurface)(double, double, double, double);
#define MAXX 50

static int N = 2, tmpN = N;
const double step = 0.01;
static double P = 1.5, tmpP = P, R = 0.95, tmpR = R, M = 0.15,
              xEnd = MAXX / P + step;
const double xBegin = 0, yBegin = 1, zBegin = 1.5;
static double E1 = 0, E2 = 0, E3 = 0;
static double U1 = 0, U2 = 2.5, U3 = 5;
static double Y1 = 1, Y2 = 1, Y3 = 3;

double f(double z) { return z; }

double g(double x) { return M * cos(x) - P; }

double surfaceFoo(double x, double E, double Y, double U) {
    return E - M * Y * cos(x - U);
}

double findSecond(const std::vector<std::pair<double, double>>& points,
                  double first) {
    for (const auto& point : points) {
        if (std::abs(point.first - first) < 1e-9) {
            return point.second;
        }
    }
    return 0.0;
}

double calculate_derivative(
    const std::vector<std::pair<double, double>>& points, double x) {
    return (findSecond(points, x + step) - findSecond(points, x - step)) /
           (2 * step);
}

std::vector<std::pair<double, double>> rungeKutta(
    double x0, double y0, double z0, funcF _f, funcG _g,
    std::vector<std::pair<double, double>> _surface) {
    auto findYOnSurface = [&](double x) {
        double closestY = _surface[0].second;
        double minDist = std::fabs(_surface[0].first - x);
        for (const auto& point : _surface) {
            double dist = std::fabs(point.first - x);
            if (dist < minDist) {
                minDist = dist;
                closestY = point.second;
            }
        }
        return closestY;
    };
    std::vector<std::pair<double, double>> points;
    while (x0 <= xEnd) {
        points.push_back(std::make_pair(x0, y0));
        double k1 = step * _f(z0);
        double l1 = step * _g(x0);
        double k2 = step * _f(z0 + l1 / 2);
        double l2 = step * _g(x0 + step / 2);
        double k3 = step * _f(z0 + l2 / 2);
        double l3 = step * _g(x0 + step / 2);
        double k4 = step * _f(z0 + l3);
        double l4 = step * _g(x0 + step);

        double nextY = y0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        double nextZ = z0 + (l1 + 2 * l2 + 2 * l3 + l4) / 6;

        double _surraceIt = findYOnSurface(x0 + step);
        double derivative = calculate_derivative(_surface, x0 + step);
        if (nextY <= _surraceIt && nextZ - derivative < 0)
            nextZ = -R * nextZ + (1 + R) * derivative;

        x0 += step;
        if (nextY < _surraceIt) {
            y0 = _surraceIt;
        } else {
            y0 = nextY;
        }
        z0 = nextZ;
    }

    return points;
}

std::vector<std::pair<double, double>> buildFoo(funcSurface _f) {
    std::vector<std::pair<double, double>> points;
    if (N == 1) {
        for (double x = 0; x <= MAXX; x += step) {
            points.push_back(std::make_pair(x, _f(x, E1, Y1, U1)));
        }
    } else if (N == 2) {
        for (double x = 0; x <= MAXX; x += step) {
            points.push_back(std::make_pair(
                x, std::max(_f(x, E1, Y1, U1), _f(x, E2, Y2, U2))));
        }
    } else if (N == 3) {
        for (double x = 0; x <= MAXX; x += step) {
            points.push_back(std::make_pair(
                x, std::max(_f(x, E1, Y1, U1),
                            std::max(_f(x, E2, Y2, U2), _f(x, E3, Y3, U3)))));
        }
    }

    return points;
}

void MainWindow::buildSurface() { surface = buildFoo(surfaceFoo); }

void MainWindow::draw() {
    QVector<double> x, y, xSurf, ySurf;
    for (const auto& point : surface) {
        xSurf.append(point.first);
        ySurf.append(point.second);
    }
    double maxY = std::numeric_limits<double>::min();
    for (const auto& point : points) {
        x.append(point.first);
        y.append(point.second);
        if (maxY < point.second) maxY = point.second;
    }

    ui->widget->clearGraphs();

    ui->widget->xAxis->setRange(0, xEnd);
    ui->widget->yAxis->setRange(-maxY / 3, maxY + maxY / 3);

    ui->widget->addGraph();
    ui->widget->graph(0)->setData(x, y);
    ui->widget->graph(0)->setPen(QColor(Qt::blue));

    ui->widget->addGraph();
    ui->widget->graph(1)->setData(xSurf, ySurf);
    ui->widget->graph(1)->setPen(QColor(Qt::red));

    ui->widget->replot();
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::initializeData() {
    points = rungeKutta(xBegin, yBegin, zBegin, f, g, surface);
}

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
    buildSurface();
    initializeData();
    draw();
    ui->pValSlider->setValue(int(P * 20.0));
    ui->nValSlider->setValue(N);
    ui->rValSlider->setValue(int(R * 20.0));
    ui->rValLine->setText(QString::number(R));
    ui->pValLine->setText(QString::number(P));
    ui->nValLine->setText(QString::number(N));
}

void MainWindow::on_pValSlider_sliderMoved(int position) {
    P = position / 20.0;
    if (MAXX / P <= MAXX)
        xEnd = MAXX / P + step;
    else
        xEnd = MAXX + step;
    ui->pValLine->setText(QString::number(P));
    initializeData();
    draw();
}

void MainWindow::on_nValSlider_sliderMoved(int position) {
    N = position;
    ui->nValLine->setText(QString::number(N));
    buildSurface();
    initializeData();
    draw();
}

void MainWindow::on_pValLine_textEdited(const QString& arg1) {
    tmpP = arg1.toDouble(nullptr);
}

void MainWindow::on_pValLine_editingFinished() {
    if (tmpP > 0 && tmpP < 2) {
        P = tmpP;
        if (MAXX / P <= MAXX)
            xEnd = MAXX / P + step;
        else
            xEnd = MAXX + step;
        ui->pValSlider->setValue(int(P * 20.0));
        initializeData();
        draw();
    }
}

void MainWindow::on_nValLine_textEdited(const QString& arg1) {
    tmpN = arg1.toInt(nullptr);
}

void MainWindow::on_nValLine_editingFinished() {
    if (tmpN > 0 && tmpN < 0.2) {
        N = tmpN;
        ui->nValSlider->setValue(N);
        buildSurface();
        initializeData();
        draw();
    }
}

void MainWindow::on_rValSlider_sliderMoved(int position) {
    R = position / 20.0;
    ui->rValLine->setText(QString::number(R));
    initializeData();
    draw();
}

void MainWindow::on_rValLine_textEdited(const QString& arg1) {
    tmpR = arg1.toDouble(nullptr);
}

void MainWindow::on_rValLine_editingFinished() {
    if (tmpR > 0 && tmpR < 1) {
        R = tmpR;
        ui->rValSlider->setValue(int(R * 20.0));
        initializeData();
        draw();
    }
}
