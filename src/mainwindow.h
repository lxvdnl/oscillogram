#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QVector>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT

   public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

   private slots:
    void on_pValSlider_sliderMoved(int position);

    void on_nValSlider_sliderMoved(int position);

    void on_pValLine_textEdited(const QString &arg1);

    void on_pValLine_editingFinished();

    void on_nValLine_textEdited(const QString &arg1);

    void on_nValLine_editingFinished();

    void on_rValSlider_sliderMoved(int position);

    void on_rValLine_textEdited(const QString &arg1);

    void on_rValLine_editingFinished();

   private:
    Ui::MainWindow *ui;

    void buildSurface();
    void initializeData();
    void draw();

    std::vector<std::pair<double, double>> points, surface;
};

#endif  // MAINWINDOW_H
