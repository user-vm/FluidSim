#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QResizeEvent>

namespace Ui {
  class MainWindow;
}

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  explicit MainWindow(QWidget *parent = 0);
  ~MainWindow();

public slots:
  //void resetCamera();
  void showControls();
  void showBouyancyParamsHelp();
  void updateNumFrames();
  void updateFrameDuration();
  void updateTimestep();
  void updateFramerate();
  void bake();
  void freeBake();
  void updateVoxelSize();
  void updateSimSizeVoxel();

private:
  Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
