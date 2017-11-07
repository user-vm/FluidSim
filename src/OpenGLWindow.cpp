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

void OpenGLWindow::addPressureFrameData(){

}

void OpenGLWindow::addVelocityFrameData(){

}

void OpenGLWindow::initializePressure(std::vector<std::vector<std::vector<float>>> pInitial){

  // here we will populate the pressure grid with its initial values
  // nothing here at the moment, so grid will be all zeroes
}

void OpenGLWindow::initializeVelocity(std::vector<std::vector<std::vector<float>>> uInitial, std::vector<std::vector<std::vector<float>>> vInitial, std::vector<std::vector<std::vector<float>>> wInitial){

  // here we will populate the velocity grids with their initial values
  // nothing here at the moment, so grids will be all zeroes
}

void OpenGLWindow::bake(){

  dx = cubeSize;

  simTime = 0.0;
  size_t k=0;
  std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> args = {w,v,u,p};
  std::vector<bool> isCentered = {true, true, true, false};
  std::vector<float> outsideValues = {0,0,0,0}; // take outside velocity to be 0 for now (smoke-like)
  std::vector<std::vector<std::vector<SevenPointLagrangianMatrixElement>>> A;

  std::vector<std::vector<std::vector<float>>> r, z, s, d, precon, q;

  addPressureFrameData();
  addVelocityFrameData();
  size_t currentFrame = 0;

  // build all the grid matrices

  p.resize(2);
  u.resize(2);
  v.resize(2);
  w.resize(2);

  std::vector<std::vector<std::vector<std::vector<float>>>> centeredVarVect = {r,s,z,d,precon,q,p[0]};

  u.resize(xSimSize+1);
  v.resize(xSimSize);
  w.resize(xSimSize);

  for(int i=0;i<xSimSize;i++){
      for(auto q1: centeredVarVect)
        q1.resize(ySimSize);
      u[0][i].resize(ySimSize);
      v[0][i].resize(ySimSize+1);
      w[0][i].resize(ySimSize);
      for(int j=0;j<ySimSize;j++){
          for(auto q : centeredVarVect)
            q[i][j].resize(zSimSize);

          u[0][i][j].resize(zSimSize);
          v[0][i][j].resize(zSimSize);
          w[0][i][j].resize(zSimSize+1);

        }
      v[ySimSize].resize(zSimSize);
    }

  for(int j=0;j<ySimSize;j++)
    u[0][xSimSize][j].resize(zSimSize);

  u[1] = u[0];
  v[1] = v[0];
  w[1] = w[0];

  p[1] = p[0];

  initializePressure(p[0]);
  initializeVelocity(u[0],v[0],w[0]);

  // prepare A matrix; since walls are not implemented yet, all diag values will be 6, and all other values will be -1, except at the (xSimSize-1, ySimSize-1, zSimSize-1) corner

  for(int ai=0;ai<xSimSize;ai++)
    for(int aj=0;aj<ySimSize;aj++)
      for(int ak=0;ak<zSimSize;ak++){

          A[ai][aj][ak].diag = 6;

          if(ai<xSimSize-1)
            A[ai][aj][ak].iUp = -1;
          else
            A[ai][aj][ak].iUp = 0;

          if(aj<ySimSize-1)
            A[ai][aj][ak].jUp = -1;
          else
            A[ai][aj][ak].jUp = 0;

          if(ak<zSimSize-1)
            A[ai][aj][ak].kUp = -1;
          else
            A[ai][aj][ak].kUp = 0;

        }

  while(simTime<=totalSimTime-dt){

      // only need information for two consectutive time steps
      // ADVECT RUNS WITH [0] AS OLD, BODY REWRITES [1], PROJECT RUNS WITH [1] AS OLD, REPEAT
      advect(args,isCentered,outsideValues);

      // body function is just updating the velocities to account for gravity
      for(int i=0;i<=xSimSize;i++)
        for(int j=0;j<=ySimSize;j++)
          for(int k=0;k<=zSimSize;k++){

            if(i==xSimSize){
              u[1][i][j][k] += g.m_x * dt;
              continue;}
            else
              if(j==ySimSize){
                v[1][i][j][k] += g.m_y * dt;
                continue;}
              else
                if(k==zSimSize){
                    w[1][i][j][k] += g.m_z * dt;
                    continue;}

            u[1][i][j][k] += g.m_x * dt;
            v[1][i][j][k] += g.m_y * dt;
            w[1][i][j][k] == g.m_z * dt;
            }

      project(A,z,d,r,s,precon,q);
      simTime += dt;

      if(fmod(simTime,dt) >= currentFrame){
          addPressureFrameData();
          addVelocityFrameData();
          currentFrame++;
        }
    }
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
