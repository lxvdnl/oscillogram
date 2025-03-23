//!!!!
//!!!!!

#include "mainwindow.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

#include "ui_mainwindow.h"

typedef double (*funcF)(double);
typedef double (*funcG)();
typedef double (*funcSurface)(double, double, double, double);
#define MAXX 50

static int N = 1, tmpN = N;
const double step = 0.001;
static double P = 1.5, tmpP = P, R = 0.4, tmpR = R, M = 0.1,
              xEnd = MAXX / P + step;
const double xBegin = 0, yBegin = 1, zBegin = 1.5;
static double E1 = 0.02;
static double U1 = 0.3;
static double Y1 = 4;

double f(double z) { return z; }

double g() { return -P; }

double surfaceFoo(double x, double E, double Y, double U) {
    return E - M * Y * cos(x - U);
}

double derivativeSurfaceFoo(double x, double E, double Y, double U) {
    return M * Y * sin(x - U);
}

double interpolate(double x, const std::vector<std::pair<double, double>>& points) {
    for (size_t i = 0; i < points.size() - 1; ++i) {
        if (points[i].first <= x && x <= points[i + 1].first) {
            double x1 = points[i].first, y1 = points[i].second;
            double x2 = points[i + 1].first, y2 = points[i + 1].second;
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
        }
    }
    return points.back().second;
}

double calculate_interpolated_derivative(double x, const std::vector<std::pair<double, double>>& surface) {
    double y1 = interpolate(x - step, surface);
    double y2 = interpolate(x + step, surface);
    return (y2 - y1) / (2 * step);
}

std::vector<std::pair<double, double>> rungeKutta(
    double x0, double y0, double z0, funcF _f, funcG _g,
    double step, double tolerance, double minStep, double maxStep) {

    std::vector<std::pair<double, double>> points;

    while (x0 <= xEnd) {
        points.push_back(std::make_pair(x0, y0));

        double currentStep = step;
        double nextY, nextZ;

        while (true) {
            // Рунге-Кутта 4-го порядка
            double k1_y = currentStep * _f(z0);
            double k1_z = currentStep * _g();
            double k2_y = currentStep * _f(z0 + k1_z / 2);
            double k2_z = currentStep * _g();
            double k3_y = currentStep * _f(z0 + k2_z / 2);
            double k3_z = currentStep * _g();
            double k4_y = currentStep * _f(z0 + k3_z);
            double k4_z = currentStep * _g();

            double tempY = y0 + (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;
            double tempZ = z0 + (k1_z + 2 * k2_z + 2 * k3_z + k4_z) / 6;

            // Оценка погрешности
            double errorY = std::abs(tempY - y0);
            double errorZ = std::abs(tempZ - z0);

            if (errorY <= tolerance && errorZ <= tolerance) {
                // Условие точности выполнено
                nextY = tempY;
                nextZ = tempZ;
                break;
            } else if (currentStep > minStep) {
                // Уменьшаем шаг
                currentStep /= 2.0;
            } else {
                // Минимальный шаг достигнут
                nextY = tempY;
                nextZ = tempZ;
                break;
            }
        }

        // Проверка пересечения поверхности
        double surfaceValue = surfaceFoo(x0, E1, Y1, U1);
        double df_dtau = derivativeSurfaceFoo(x0, E1, Y1, U1);

        if (nextY <= surfaceValue && nextZ - df_dtau < 0) {
            // Поиск точки пересечения методом линейной интерполяции
            double t = (surfaceValue - y0) / (nextY - y0);
            double xImpact = x0 + t * currentStep;
            double yImpact = surfaceValue;
            double zImpact = z0 + t * (nextZ - z0);

            // Обновляем значения на момент удара
            nextY = yImpact;
            nextZ = -R * zImpact + (1 + R) * df_dtau;

            // Добавляем точку удара
            points.push_back(std::make_pair(xImpact, yImpact));
        }

        // Переходим к следующему шагу
        x0 += currentStep;
        y0 = nextY;
        z0 = nextZ;

        // Увеличиваем шаг, если можно
        if (currentStep < maxStep) {
            currentStep *= 2.0;
        }
    }

    return points;
}


std::vector<double> rungeKuttaWithImpacts(double x0, double y0, double z0, funcF _f, funcG _g,
                                          double E, double Y, double U, double curP,
                                          double step, double tolerance, double minStep, double maxStep) {
    std::vector<double> postImpactSpeeds;
    int impacts = 0;
    P = curP;

    while (impacts < 1200) {
        double currentStep = step;
        bool impactOccurred = false;
        double nextY, nextZ;

        while (true) {
            // Рунге-Кутта 4-го порядка
            double k1_y = currentStep * _f(z0);
            double k1_z = currentStep * _g();
            double k2_y = currentStep * _f(z0 + k1_z / 2);
            double k2_z = currentStep * _g();
            double k3_y = currentStep * _f(z0 + k2_z / 2);
            double k3_z = currentStep * _g();
            double k4_y = currentStep * _f(z0 + k3_z);
            double k4_z = currentStep * _g();

            double tempY = y0 + (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;
            double tempZ = z0 + (k1_z + 2 * k2_z + 2 * k3_z + k4_z) / 6;

            // Точность: оценка ошибки (простая модель)
            double errorY = std::abs(tempY - y0);
            double errorZ = std::abs(tempZ - z0);

            if (errorY <= tolerance && errorZ <= tolerance) {
                // Условие выполнено: шаг подходит
                nextY = tempY;
                nextZ = tempZ;
                break;
            } else if (currentStep > minStep) {
                // Уменьшаем шаг
                currentStep /= 2.0;
            } else {
                // Минимальный шаг достигнут, выходим
                nextY = tempY;
                nextZ = tempZ;
                break;
            }
        }

        // Проверка на пересечение поверхности
        double surfaceValue = surfaceFoo(x0, E, Y, U);
        double df_dtau = derivativeSurfaceFoo(x0, E, Y, U);

        if (nextY <= surfaceValue && nextZ - df_dtau < 0) {
            // Поиск точки пересечения методом линейной интерполяции
            double t = (surfaceValue - y0) / (nextY - y0);
            double yImpact = surfaceValue;
            double zImpact = z0 + t * (nextZ - z0);

            // Обновление значений на точке удара
            nextY = yImpact;
            nextZ = -R * zImpact + (1 + R) * df_dtau;
            impactOccurred = true;
        }

        // Обновляем параметры
        x0 += currentStep;
        y0 = nextY;
        z0 = nextZ;

        // Добавляем скорость удара, если он произошёл
        if (impactOccurred) {
            postImpactSpeeds.push_back(nextZ);
            impacts++;
        }

        // Ограничиваем размер массива скоростей
        if (postImpactSpeeds.size() > 200) {
            postImpactSpeeds.erase(postImpactSpeeds.begin());
        }

        // Увеличиваем шаг, если можно
        if (!impactOccurred && currentStep < maxStep) {
            currentStep *= 2.0;
        }
    }

    return postImpactSpeeds;
}


std::vector<std::pair<double, double>> buildFoo(funcSurface _f) {
    std::vector<std::pair<double, double>> points;
    for (double x = 0; x <= MAXX; x += step) {
        points.push_back(std::make_pair(x, _f(x, E1, Y1, U1)));
    }
    return points;
}

std::vector<std::pair<double, double>> bifurcationDiagram(std::function<void(int)> progressCallback, double pStart, double pEnd, double pStep) {
    std::vector<std::pair<double, double>> bifurcationData;

    int totalSteps = (pEnd - pStart) / pStep;
    int currentStep = 0;

    for (double currentP = pStart; currentP <= pEnd; currentP += pStep) {
        double x0 = xBegin;
        double y0 = yBegin;
        double z0 = zBegin;

        std::vector<double> impacts = rungeKuttaWithImpacts(x0, y0, z0, f, g, E1, Y1, U1, currentP, 0.5, 1e-3, 1e-6, 1.0);
        std::cout << "Size = " << impacts.size() << "  first val = " << impacts[0] << " , second val = " << impacts[1] << "\n";

        for (double impactValue : impacts) {
            bifurcationData.push_back({currentP, impactValue});
        }

        currentStep++;
        progressCallback(static_cast<int>(100.0 * currentStep / totalSteps));
    }

    return bifurcationData;
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

void MainWindow::plotBifurcationDiagram() {
    ui->progressBar->setRange(0, 100);
    ui->progressBar->setValue(0);
    double pStart = 0.196, pEnd = 0.218, pStep = 0.00002;

    auto data = bifurcationDiagram([this](int progress) {
        ui->progressBar->setValue(progress);
        QCoreApplication::processEvents();
    }, pStart, pEnd, pStep);

    auto minmax = std::minmax_element(data.begin(), data.end(),
                                          [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                                              return a.second < b.second;
                                          });

    QVector<double> x, y;
    for (const auto& point : data) {
        x.append(point.first);
        y.append(point.second);
    }

    ui->widget->clearGraphs();

    ui->widget->xAxis->setRange(pStart - 0.001, pEnd + 0.001);
    ui->widget->yAxis->setRange(minmax.first->second - 0.01, minmax.second->second + 0.01);

    ui->widget->addGraph();

    ui->widget->graph(0)->setData(x, y);

    ui->widget->graph(0)->setPen(QColor(Qt::black));
    ui->widget->graph(0)->setLineStyle(QCPGraph::lsNone);

    QCPScatterStyle scatterStyle(QCPScatterStyle::ssCircle);
    scatterStyle.setSize(2);
    scatterStyle.setPen(QPen(Qt::black));
    scatterStyle.setBrush(QBrush(Qt::black));
    ui->widget->graph(0)->setScatterStyle(scatterStyle);

    ui->widget->replot();
    ui->progressBar->setValue(100);

}



MainWindow::~MainWindow() { delete ui; }

void MainWindow::initializeData() {
    points = rungeKutta(xBegin, yBegin, zBegin, f, g, 0.1, 1e-6, 0.00001, 1.0);
}

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
    std::cout << "UI Loaded" << std::endl;
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

void MainWindow::on_bifurcationButton_clicked() {
    plotBifurcationDiagram();
}
