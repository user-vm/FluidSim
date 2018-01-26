#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include <sstream>
#include <iomanip>
#include <unordered_map>

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow)
{
  ui->setupUi(this);
}

MainWindow::~MainWindow()
{
  delete ui;
}

void MainWindow::showControls()
{
  QMessageBox controlBox;
  QFont courierFont("Courier", 10, QFont::Normal, false);

  controlBox.setFont(courierFont);

  std::stringstream controls;
  std::vector<std::vector<std::string>> controlAssoc =
  {{"Left Click + Drag","Rotate around lookAt point"},
   {"Middle Click + Drag","Move lookAt point"},
   {"Right Click + Drag","Zoom in/out"},
   {"Scroll","Zoom in/out"},
   {"Left key","Decrement frame"},
   {"Right key","Increment frame"},
   {"Space","Pause/Play"}
  };

  for(uint i=0;i<controlAssoc.size();i++)
    controls<<std::setfill(' ')<<std::setw(26)<<std::left<<controlAssoc[i][0]<<std::setfill(' ')<<std::setw(27)<<std::right<<controlAssoc[i][1]<<'\n';

  std::string controlsString;
  controlsString = controls.str();

  controlBox.setText(controlsString.c_str());

  controlBox.exec();
}

void MainWindow::showBouyancyParamsHelp(){
  QMessageBox bouyancyParamsHelpBox;

  bouyancyParamsHelpBox.setText("The bouyancy force is the vector (0, -as + b(T - T_ambient), 0), where \'s\' is the smoke concentration.");

  bouyancyParamsHelpBox.exec();
}

void MainWindow::updateNumFrames(){

  bool ok;

  uint numFrames = ui->framesText->text().toUInt(&ok);

  if(ok){
      if(numFrames==0){
        numFrames = 1;
        ui->framesText->setText("1");
        }
      ui->openGLWidget->totalSimTime = (numFrames-1) * ui->openGLWidget->frameDuration;
    }
  else
    ui->framesText->setText(QString(std::to_string(uint(ui->openGLWidget->totalSimTime/ui->openGLWidget->frameDuration)+1).c_str()));
}

void MainWindow::updateFrameDuration(){

  bool ok;

  float frameDuration = ui->frameLengthText->text().toFloat(&ok);

  if(ok && frameDuration > 0){
      if(frameDuration < 0.005){
        frameDuration = 0.005;
        ui->frameLengthText->setText("0.005");
}
      ui->openGLWidget->totalSimTime *= frameDuration / ui->openGLWidget->frameDuration; //new frame duration over old
      ui->openGLWidget->frameDuration = frameDuration;
      ui->framerateText->setText(QString(std::to_string(1.0/frameDuration).c_str()));
    }
  else
    ui->framesText->setText(QString(std::to_string(ui->openGLWidget->frameDuration).c_str()));
}

void MainWindow::updateFramerate(){

  bool ok;

  float framerate = ui->framerateText->text().toFloat(&ok);

  if(ok){
      ui->frameLengthText->setText(QString(std::to_string(1.0/framerate).c_str()));
      updateFrameDuration();
    }
}

void MainWindow::updateTimestep(){

  bool ok;

  float timestep = ui->timestepText->text().toFloat(&ok);

  if(ok && timestep>0.0){
      ui->openGLWidget->dt = timestep;
    }
  else
    ui->timestepText->setText(QString(std::to_string(ui->openGLWidget->dt).c_str()));
}

void MainWindow::bake(){

  ui->tab_1->setEnabled(false);
  ui->tab_2->setEnabled(false);
  ui->bakeButton->setEnabled(false);
  ui->frameLengthText->setEnabled(false);
  ui->framerateText->setEnabled(false);
  ui->framesText->setEnabled(false);
  ui->timestepText->setEnabled(false);

  ui->openGLWidget->reset();
}

void MainWindow::freeBake(){

  ui->tab_1->setEnabled(true);
  ui->tab_1->setEnabled(true);
  ui->tab_2->setEnabled(true);
  ui->bakeButton->setEnabled(true);
  ui->frameLengthText->setEnabled(true);
  ui->framerateText->setEnabled(true);
  ui->framesText->setEnabled(true);
  ui->timestepText->setEnabled(true);

  ui->openGLWidget->freeBake();
}

void MainWindow::updateVoxelSize(){

  bool ok;

  float voxelSize = ui->timestepText->text().toFloat(&ok);

  if(ok && voxelSize>0.0){
      ui->openGLWidget->dx = voxelSize;
    }
  else
    ui->timestepText->setText(QString(std::to_string(ui->openGLWidget->dx).c_str()));
}

void MainWindow::updateSimSizeVoxel(){

  bool okX, okY, okZ;

  uint voxelSizeX = ui->simSizeVoxXText->text().toUInt(&okX);
  uint voxelSizeY = ui->simSizeVoxYText->text().toUInt(&okY);
  uint voxelSizeZ = ui->simSizeVoxZText->text().toUInt(&okZ);

  if(okX && okY && okZ && voxelSizeX>0 && voxelSizeY>0 && voxelSizeZ>0){

      //0 -> Scale, 1 -> Center, 2 -> Do not move
      ui->openGLWidget->resizeSimulation(voxelSizeX, voxelSizeY, voxelSizeZ, ui->resizeMethodBox->currentIndex());

      ui->openGLWidget->xSimSize = voxelSizeX;
      ui->openGLWidget->ySimSize = voxelSizeY;
      ui->openGLWidget->zSimSize = voxelSizeZ;
    }
  else{
      ui->simSizeVoxXText->setText(QString(std::to_string(ui->openGLWidget->xSimSize).c_str()));
      ui->simSizeVoxYText->setText(QString(std::to_string(ui->openGLWidget->ySimSize).c_str()));
      ui->simSizeVoxZText->setText(QString(std::to_string(ui->openGLWidget->zSimSize).c_str()));
    }
}
