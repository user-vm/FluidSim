/*
 * Basic GL Window modified from the example here
 * http://qt-project.org/doc/qt-5.0/qtgui/openglwindow.html
 * adapted to use NGL
 */
#include "OpenGLWindow.h"
#include <QKeyEvent>
#include <QApplication>
#include <memory>
#include <iostream>
#include <ngl/Vec3.h>
#include <QElapsedTimer>
#include <cmath>

constexpr float cubeSize=0.05;

OpenGLWindow::OpenGLWindow()
{
  setTitle("Qt5 compat profile OpenGL 3.2");

}

OpenGLWindow::~OpenGLWindow()
{
  // now we have finished clear the device
  std::cout<<"deleting buffer\n";
  glDeleteBuffers(1,&m_vboPointer);
}

void OpenGLWindow::initializeGL()
{

  glewInit();

  glClearColor(0.5f, 0.5f, 0.5f, 1.0f);			   // Grey Background

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_MULTISAMPLE);

  makeCubes(cubeSize);
  bake();
}

void  OpenGLWindow::makeCubes( GLfloat _size)
{
  // allocate enough space for our verts
  // as we are doing lines it will be 2 verts per line
  // and we need to add 1 to each of them for the <= loop
  // and finally muliply by 12 as we have 12 values per line pair
  qint64 k=-1;
  m_cubeSubVBOSize=3*3*2*6*3;//3 coordinates per vertex, 3 vertices per tri, 2 tris per quad, 6 quads per cube, 3 data types (vertex, normal, colour)
  m_vboSize= m_cubeSubVBOSize*xSimSize*ySimSize*zSimSize*3;

  std::unique_ptr<GLfloat []>vertexData( new GLfloat[m_vboSize]);

  std::vector<ngl::Vec3> verts=
  {
    //12 triangles, two for each face
    //face z=-0.5
    ngl::Vec3(0.5,0.5,-0.5),
    ngl::Vec3(-0.5,0.5,-0.5),
    ngl::Vec3(-0.5,-0.5,-0.5), //+1
    ngl::Vec3(0.5,0.5,-0.5),
    ngl::Vec3(0.5,-0.5,-0.5),
    ngl::Vec3(-0.5,-0.5,-0.5), //-1
    //face z=0.5
    ngl::Vec3(0.5,0.5,0.5),
    ngl::Vec3(-0.5,0.5,0.5),
    ngl::Vec3(-0.5,-0.5,0.5), //
    ngl::Vec3(0.5,0.5,0.5),
    ngl::Vec3(0.5,-0.5,0.5),
    ngl::Vec3(-0.5,-0.5,0.5),
    //face x=-0.5
    ngl::Vec3(-0.5,0.5,0.5),
    ngl::Vec3(-0.5,-0.5,0.5),
    ngl::Vec3(-0.5,-0.5,-0.5),
    ngl::Vec3(-0.5,0.5,0.5),
    ngl::Vec3(-0.5,0.5,-0.5),
    ngl::Vec3(-0.5,-0.5,-0.5),
    //face x=0.5
    ngl::Vec3(0.5,0.5,0.5),
    ngl::Vec3(0.5,-0.5,0.5),
    ngl::Vec3(0.5,-0.5,-0.5),
    ngl::Vec3(0.5,0.5,0.5),
    ngl::Vec3(0.5,0.5,-0.5),
    ngl::Vec3(0.5,-0.5,-0.5),
    //face y=-0.5
    ngl::Vec3(0.5,-0.5,0.5),
    ngl::Vec3(-0.5,-0.5,0.5),
    ngl::Vec3(-0.5,-0.5,-0.5),
    ngl::Vec3(0.5,-0.5,0.5),
    ngl::Vec3(0.5,-0.5,-0.5),
    ngl::Vec3(-0.5,-0.5,-0.5),
    //face y=0.5
    ngl::Vec3(0.5,0.5,0.5),
    ngl::Vec3(-0.5,0.5,0.5),
    ngl::Vec3(-0.5,0.5,-0.5),
    ngl::Vec3(0.5,0.5,0.5),
    ngl::Vec3(0.5,0.5,-0.5),
    ngl::Vec3(-0.5,0.5,-0.5),
  };

  ngl::Vec3 normalHolder;
  bool xCoeff=false, yCoeff=false, zCoeff=false;

  //we will interleave the vertex, normal and colour data

  for(size_t d1=0;d1<xSimSize;d1++)
    for(size_t d2=0;d2<ySimSize;d2++)
      for(size_t d3=0;d3<zSimSize;d3++)
        for(size_t i=0;i<verts.size();i+=3){
            //figure out what the normals will be

            xCoeff=false;
            yCoeff=false;
            zCoeff=false;

            if(verts[i].m_x==verts[i+1].m_x&&verts[i+1].m_x==verts[i+2].m_x)
              xCoeff = true;
            else if(verts[i].m_y==verts[i+1].m_y&&verts[i+1].m_y==verts[i+2].m_y)
              yCoeff = true;
            else if(verts[i].m_z==verts[i+1].m_z&&verts[i+1].m_z==verts[i+2].m_z)
              zCoeff = true;

            //add vertices in three by three, because each vertex in a tri has the same normal vector

            for(int t=0;t<3;t++){
                //vertex position
                vertexData[++k]=(verts[i+t].m_x+d1*1 - xSimSize/2.0)*_size;
                vertexData[++k]=(verts[i+t].m_y+d2*1 - ySimSize/2.0)*_size;
                vertexData[++k]=(verts[i+t].m_z+d3*1 - zSimSize/2.0)*_size;

                //vertex normal; factor of two because each quad center is 0.5 away from cube origin
                vertexData[++k]=verts[i+t].m_x * 2 * xCoeff;
                vertexData[++k]=verts[i+t].m_y * 2 * yCoeff;
                vertexData[++k]=verts[i+t].m_z * 2 * zCoeff;

                //vertex color
                vertexData[++k]=0.5 + d1*0.5/xSimSize;// + verts[i+t].m_x;
                vertexData[++k]=0.5 + d2*0.5/ySimSize;// + verts[i+t].m_y;
                vertexData[++k]=0.5 + d3*0.5/zSimSize;// + verts[i+t].m_z;

                std::cout<<vertexData[k]<<"\n";

              }

          }

  // now we will create our VBO first we need to ask GL for an Object ID

  for(size_t i=0;i<verts.size();i++){
      std::cout<<verts[i].m_x<<" "<<verts[i].m_y<<" "<<verts[i].m_z<<"\n";
    }

  glGenBuffers(1, &m_vboPointer);
  // now we bind this ID to an Array buffer
  glBindBuffer(GL_ARRAY_BUFFER, m_vboPointer);
  // finally we stuff our data into the array object
  // First we tell GL it's an array buffer
  // then the number of bytes we are storing (need to tell it's a sizeof(FLOAT)
  // then the pointer to the actual data
  // Then how we are going to draw it (in this case Statically as the data will not change)
  glBufferData(GL_ARRAY_BUFFER, m_vboSize*sizeof(GL_FLOAT) , vertexData.get(), GL_DYNAMIC_DRAW);

  startTimer(100);
  timer.start();

}

void OpenGLWindow::paintGL()
{

  // set the viewport
  glViewport(m_xOffset,m_yOffset,m_width,m_height);
  // clear the colour and depth buffers ready to draw.
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // enable  vertex array drawing
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  // bind our VBO data to be the currently active one

  //glUseProgram(shaderProgramID);

  glVertexPointer(3,GL_FLOAT,9*sizeof(GL_FLOAT),(void*)0); //SET THE STRIDE
  glNormalPointer(GL_FLOAT,9*sizeof(GL_FLOAT),(void*)(3*sizeof(GL_FLOAT)));
  glColorPointer(3,GL_FLOAT,9*sizeof(GL_FLOAT),(void*)(6*sizeof(GL_FLOAT)));

  glPushMatrix();
  glRotatef(360.0* (m_spin?(timer.elapsed()/5000.0):1),0.0,1.0,0.0);  //* timer.elapsed()/5000
  glRotatef(45.0,-0.4,1.0,0.0);
  glDrawArrays(GL_TRIANGLES, 0, xSimSize*ySimSize*zSimSize*m_cubeSubVBOSize*10);
  glPopMatrix();
  // now turn off the VBO client state as we have finished with it
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
}

void OpenGLWindow::timerEvent(QTimerEvent *)
{
  update();
}

void OpenGLWindow::keyPressEvent(QKeyEvent *_event)
{
  switch (_event->key())
    {
    case Qt::Key_Escape : QApplication::exit(EXIT_SUCCESS); break;
    }
}

void OpenGLWindow::resizeGL(int _w, int _h)
{

  m_width  = static_cast<int>( ((_w<_h)?_w:_h) * devicePixelRatio() );
  m_height = static_cast<int>( ((_w<_h)?_w:_h) * devicePixelRatio() );
  m_xOffset = (_w - m_width)/2;
  m_yOffset = (_h - m_height) / 2;
}
/*
void OpenGLWindow::addPressureFrameData(GridsHolder grids){

  pressureFrameData.resize(pressureFrameData.size() + x_Size*y_Size*z_Size);

}

void OpenGLWindow::addVelocityFrameData(){

}
*/
void OpenGLWindow::initializePressure(GridsHolder *gridsHolder){

  // here we will populate the pressure grid with its initial values

  TwoStepMatrix3D* p = gridsHolder->getTwoStepMatrix3DByName("p");

  size_t x_Size = p->xSize();
  size_t y_Size = p->ySize();
  size_t z_Size = p->zSize();

  for(size_t i=0;i<x_Size;i++)
    for(size_t j=0;j<y_Size;j++)
      for(size_t k=0;k<z_Size;k++)

        p->setOld(i,j,k,std::fabs(((x_Size*1.0)/2.0-i)/x_Size));
}

void OpenGLWindow::initializeVelocity(GridsHolder *gridsHolder){

  // here we will populate the velocity grids with their initial values
  TwoStepMatrix3D* u = gridsHolder->getTwoStepMatrix3DByName("u");
  TwoStepMatrix3D* v = gridsHolder->getTwoStepMatrix3DByName("v");
  TwoStepMatrix3D* w = gridsHolder->getTwoStepMatrix3DByName("w");

  size_t x_Size = u->xSize();
  size_t y_Size = u->ySize();
  size_t z_Size = u->zSize();

  for(size_t i=0;i<x_Size;i++)
    for(size_t j=0;j<y_Size;j++)
      for(size_t k=0;k<z_Size;k++)

        u->setOld(i,j,k,std::fabs(((x_Size*1.0)/2.0-i)/x_Size));

  x_Size = v->xSize();
  y_Size = v->ySize();
  z_Size = v->zSize();

  for(size_t i=0;i<x_Size;i++)
    for(size_t j=0;j<y_Size;j++)
      for(size_t k=0;k<z_Size;k++)

        v->setOld(i,j,k,std::fabs(((y_Size*1.0)/2.0-j)/y_Size));

  x_Size = u->xSize();
  y_Size = u->ySize();
  z_Size = u->zSize();

  for(size_t i=0;i<x_Size;i++)
    for(size_t j=0;j<y_Size;j++)
      for(size_t k=0;k<z_Size;k++)

        w->setOld(i,j,k,std::fabs(((z_Size*1.0)/2.0-k)/z_Size));

}

void OpenGLWindow::bake(){

  float dx = cubeSize;

  float simTime = 0.0;

  // need to initialize matrices

  std::vector<GridTuple*> gridsToMake = {new GridTuple("u",GRID_3D_TWOSTEP,xSimSize+1,ySimSize,zSimSize),
                                         new GridTuple("v",GRID_3D_TWOSTEP,xSimSize,ySimSize+1,zSimSize),
                                         new GridTuple("w",GRID_3D_TWOSTEP,xSimSize,ySimSize,zSimSize+1),
                                         new GridTuple("p",GRID_3D_TWOSTEP,xSimSize,ySimSize,zSimSize),
                                         new GridTuple("r",GRID_3D,xSimSize,ySimSize,zSimSize),
                                         new GridTuple("z",GRID_3D,xSimSize,ySimSize,zSimSize),
                                         new GridTuple("s",GRID_3D,xSimSize,ySimSize,zSimSize),
                                         new GridTuple("d",GRID_3D,xSimSize,ySimSize,zSimSize),
                                         new GridTuple("precon",GRID_3D,xSimSize,ySimSize,zSimSize),
                                         new GridTuple("q",GRID_3D,xSimSize,ySimSize,zSimSize)};


  std::unique_ptr<GridsHolder> grids(new GridsHolder(gridsToMake,
                                                     dx, dt, tol, maxIterations, rho, g));

  //addPressureFrameData();
  //addVelocityFrameData();
  size_t currentFrame = 0;

  initializePressure(grids.get());
  initializeVelocity(grids.get());

  pressureFrameData = std::unique_ptr<FrameData>(new FrameData(xSimSize,ySimSize,zSimSize));
  uFrameData = std::unique_ptr<FrameData>(new FrameData(xSimSize+1,ySimSize,zSimSize));
  vFrameData = std::unique_ptr<FrameData>(new FrameData(xSimSize,ySimSize+1,zSimSize));
  wFrameData = std::unique_ptr<FrameData>(new FrameData(xSimSize,ySimSize,zSimSize+1));

  pressureFrameData->addFrame(grids.get(),"p");
  uFrameData->addFrame(grids.get(),"u");
  vFrameData->addFrame(grids.get(),"v");
  wFrameData->addFrame(grids.get(),"w");

  float tempDt = dt;

  // prepare A matrix; since walls are not implemented yet, all diag values will be 6, and all other values will be -1, except at the (xSimSize-1, ySimSize-1, zSimSize-1) corner

  while(simTime<=totalSimTime-dt){

      bool makeFrame = false;

      if(size_t((simTime + dt)/frameDuration) > currentFrame){
          makeFrame = true;
          tempDt = (currentFrame + 1) * frameDuration - simTime;
        }

      // only need information for two consectutive time steps
      // ADVECT RUNS WITH [0] AS OLD, BODY REWRITES [1], PROJECT RUNS WITH [1] AS OLD, REPEAT
      grids.get()->advectDummy({"u","v","w","p"}, tempDt);

      // body function is just updating the velocities to account for gravity
      grids.get()->bodyDummy(tempDt);

      grids.get()->projectDummy(tempDt);
      simTime += tempDt;

      if(makeFrame){
          pressureFrameData->addFrame(grids.get(),"p");
          uFrameData->addFrame(grids.get(),"u");
          vFrameData->addFrame(grids.get(),"v");
          wFrameData->addFrame(grids.get(),"w");

          currentFrame++;
          std::cout<<currentFrame<<"/"<<totalSimTime/frameDuration<<"\n";
        }

      tempDt = dt;
    }

  std::cout<<"Finished bake.\n";
}

size_t OpenGLWindow::FrameData::xSize(){

  return x_Size;
}

size_t OpenGLWindow::FrameData::ySize(){

  return y_Size;
}

size_t OpenGLWindow::FrameData::zSize(){

  return z_Size;
}

bool OpenGLWindow::FrameData::addFrame(GridsHolder* gridsHolder, std::string gridName){

  return addFrame(gridsHolder, gridName, num_Frames);
}

bool OpenGLWindow::FrameData::addFrame(GridsHolder* gridsHolder, std::string gridName, size_t atFrame){

  if(atFrame > num_Frames){
      std::cout<<"Index of frame to add is not in scope of FrameData object \""<<_name<<"\"";
    return false;}

  // grid data saves onyl two timesteps; that's why you need to copy them
  // consider file caching in the future? (whatever that means)

  bool doSizesMatch = true;

  GridType gridType = gridsHolder->getTypeByName(gridName);

  //gridsHolder.getAnyByName();

  if(gridType == GRID_3D){
      Matrix3D* grid = gridsHolder->getMatrix3DByName(gridName);

      if(grid->xSize() != x_Size){
          std::cout<<"X-dimensions of FrameData object \""<<_name<<"\" and grid object \""<<gridName<<"\" do not match.\n";
          doSizesMatch = false;
        }

      if(grid->ySize() != y_Size){
          std::cout<<"Y-dimensions of FrameData object \""<<_name<<"\" and grid object \""<<gridName<<"\" do not match.\n";
          doSizesMatch = false;
        }

      if(grid->zSize() != z_Size){
          std::cout<<"Z-dimensions of FrameData object \""<<_name<<"\" and grid object \""<<gridName<<"\" do not match.\n";
          doSizesMatch = false;
        }

      if(!doSizesMatch)
        return false;

      data.resize(data.size() + x_Size*y_Size*z_Size);

      size_t pos = atFrame * x_Size * y_Size * z_Size;

      for(size_t i=0;i<x_Size;i++)
        for(size_t j=0;j<y_Size;j++)
          for(size_t k=0;k<z_Size;k++,pos++)

            data[pos] = grid->get(i,j,k);

      num_Frames++;
      return true;
    }

  if(gridType == GRID_3D_TWOSTEP){
      TwoStepMatrix3D* grid = gridsHolder->getTwoStepMatrix3DByName(gridName);

      if(grid->xSize() != x_Size){
          std::cout<<"X-dimensions of FrameData object \""<<_name<<"\" and grid object \""<<gridName<<"\" do not match.\n";
          doSizesMatch = false;
        }

      if(grid->ySize() != y_Size){
          std::cout<<"Y-dimensions of FrameData object \""<<_name<<"\" and grid object \""<<gridName<<"\" do not match.\n";
          doSizesMatch = false;
        }

      if(grid->zSize() != z_Size){
          std::cout<<"Z-dimensions of FrameData object \""<<_name<<"\" and grid object \""<<gridName<<"\" do not match.\n";
          doSizesMatch = false;
        }

      if(!doSizesMatch)
        return false;

      data.resize(data.size() + x_Size*y_Size*z_Size);

      size_t pos = atFrame * x_Size * y_Size * z_Size;

      for(size_t i=0;i<x_Size;i++)
        for(size_t j=0;j<y_Size;j++)
          for(size_t k=0;k<z_Size;k++,pos++)

            data[pos] = grid->getOld(i,j,k);

      num_Frames++;

      return true;
    }

  std::cout<<"Grid \""<<gridName<<"\" is neither of type Matrix3D nor TwoStepMatrix3D.\n";
  return false;
}

OpenGLWindow::FrameData::FrameData(size_t xSize, size_t ySize, size_t zSize){

  x_Size = xSize;
  y_Size = ySize;
  z_Size = zSize;
  num_Frames = 0;
}
