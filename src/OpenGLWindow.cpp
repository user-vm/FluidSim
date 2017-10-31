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
}

bool OpenGLWindow::advect(){
  std::vector args;
  std::vector<bool> isCentered;
  return advect(args, isCentered);
}

bool OpenGLWindow::advect(std::vector args, std::vector<bool> isCentered, std::vector<float> outsideValues)
{
  // if args is empty, assume only u,v,w and p need to be updated
  if(args.size()==0){
      isCentered.clear();

      args.push_back(w);
      args.push_back(false);

      args.push_back(v);
      args.push_back(false);

      args.push_back(u);
      args.push_back(false);

      args.push_back(p);
      args.push_back(true);
    }
  else
    if(args.size() != isCentered.size()){
        std::cout<<"args and isCentered mismatch";
        return false;
      }
  int k=0;
  for(auto c:args){

      if(isCentered[k++]==0){

          ngl::Vec3 velocity;

          for(ix=0;ix<xSimSize;ix++)
            for(iy=0;iy<ySimSize;iy++)
              for(iz=0;iz<zSimSize;iz++){

                  // get velocity components at center of current grid cell through trilinear interpolation
                  velocity.m_x = (u[ix][iy][iz]+u[ix+1][iy][iz])/2;
                  velocity.m_y = (v[ix][iy][iz]+v[ix][iy+1][iz])/2;
                  velocity.m_z = (w[ix][iy][iz]+w[ix][iy][iz+1])/2;

                  // the current grid cell (of indices (ix,iy,iz)) is viewed as a particle
                  // (xn,yn,zn) is the corresponding position of this particle; (xp,yp,zp) is its projected past position a time dt ago
                  // we use Forward Euler (might update to RK2)
                  // (ixp,iyp,izp) are the indices of the grid point (xpFloor,ypFloor,zpFloor) such that xp is in [xpFloor, xpFloor+dx), yp is in [ypFloor, ypFloor+dx), zp is in [zpFloor, zpFloor+dx)
                  xn = dx * (ix + 0.5);
                  xp = xn - dt * velocity.m_x;
                  ixp = floor(xp/dx - 0.5);

                  yn = dx * (iy + 0.5);
                  yp = yn - dt * velocity.m_y;
                  iyp = floor(yp/dx - 0.5);

                  zn = dx * (iz + 0.5);
                  zp = zn - dt * velocity.m_z;
                  izp = floor(zp/dx - 0.5);

                  for(i=0;i<2;i++)
                    for(j=0;j<2;j++)
                      for(k=0;k<2;k++){
                          if(ixp+i>=xSimSize || iyp+j>=ySimSize || izp+k>=zSimSize) // outside the simulation volume, there is a constant value for the quantity
                            cNew[x][y][z] += outsideValues[k];
                          cNew[x][y][z]+=(i?ax:(1-ax))*(j?ay:(1-ay))*(k?(1-az):az)*c[ixp+i][iyp+j][izp+k];
                        }
                }
        }
      else{

          ngl::Vec3 uVelocity, vVelocity, wVelocity;
          size_t xSize,ySize,zSize;

          xSize = c.size();
          ySize = c[0].size();
          zSize = c[0][0].size();

          for(ix=0;ix<xSize;ix++)
            for(iy=0;iy<ySize;iy++)
              for(iz=0;iz<zSize;iz++){

                  // don't need interpolation for xp, yp, zp, since they are at grid cell edges, just like the velocities

                  // the current grid wall (of indices (ix,iy,iz)) is viewed as a particle; the wall plane depends on the current data grid
                  // (xn,yn,zn) is the corresponding position of this particle; (xp,yp,zp) is its projected past position a time dt ago
                  // we use Forward Euler (might update to RK2)
                  // (ixp,iyp,izp) are the indices of the grid wall (xpFloor,ypFloor,zpFloor) such that xp is in [xpFloor, xpFloor+dx), yp is in [ypFloor, ypFloor+dx), zp is in [zpFloor, zpFloor+dx)
                  xn = dx * ix;
                  xp = xn - dt * u[x][y][z];
                  ixp = floor(xp/dx);

                  yn = dx * iy;
                  yp = yn - dt * v[x][y][z];
                  iyp = floor(yp/dx);

                  zn = dx * iz;
                  zp = zn - dt * w[x][y][z];
                  izp = floor(zp/dx);

                  for(i=0;i<2;i++)
                    for(j=0;j<2;j++)
                      for(k=0;k<2;k++){
                          if(ixp+i>=xSize || iyp+j>=ySize || izp+k>=zSize) // outside the simulation volume, there is a constant value for the quantity
                            cNew[x][y][z] += outsideValues[k];
                          cNew[x][y][z]+=(i?ax:(1-ax))*(j?ay:(1-ay))*(k?(1-az):az)*c[ixp+i][iyp+j][izp+k];
                        }

                }


        }
    }
}

void OpenGLWindow::project()
{

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

  //the simulation uses a staggered MAC grid
  //the p values are at the center of the grid cells (p[i][j][k] = p_i,j,k)
  //p = new GLfloat[xSimSize][ySimSize][zSimSize];

  p.resize(xSimSize);
  tm.resize(xSimSize);
  u.resize(xSimSize+1);
  v.resize(xSimSize);
  w.resize(xSimSize);

  for(int i=0;i<xSimSize;i++){
      p[i].resize(ySimSize);
      tm[i].resize(ySimSize);
      u[i].resize(ySimSize);
      v[i].resize(ySimSize+1);
      w[i].resize(ySimSize);
      for(int j=0;j<ySimSize;j++){
          p[i][j].resize(zSimSize);
          tm[i][j].resize(zSimSize);
          u[i][j].resize(zSimSize);
          v[i][j].resize(zSimSize);
          w[i][j].resize(zSimSize+1);}
    }

  for(int i=0;i<xSimSize;i++){
      v[i][ySimSize].resize(zSimSize);
    }

  u[xSimSize].resize(zSimSize);
  for(int j=0;j<ySimSize;j++)
    u[xSimSize][j].resize(zSimSize);

  //the u,v,w values are at the centers of the faces that delimit cells (u[i][j][k] = u_i-1/2,j,k; v[i][j][k] = v_i,j-1/2,k; w[i][j][k] = w_i,j,k-1/2)
  //u = new float[xSimSize+1][ySimSize][zSimSize];
  //v = new float[xSimSize][ySimSize+1][zSimSize];
  //w = new float[xSimSize][ySimSize][zSimSize+1];

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

  if(m_spin == true){
      startTimer(100);
      timer.start();}

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
