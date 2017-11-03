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

bool OpenGLWindow::project(std::vector<std::vector<std::vector<SevenPointLagrangianMatrixElement>>> A, std::vector<std::vector<std::vector<float>>> z,
                           std::vector<std::vector<std::vector<float>>> d, std::vector<std::vector<std::vector<float>>> r,
                           std::vector<std::vector<std::vector<float>>> s, std::vector<std::vector<std::vector<float>>> precon,
                           std::vector<std::vector<std::vector<float>>> q)
{
  int oldIndex = 1;
  int newIndex = 0;

  float sigma;

  bool breakIteration = true; // flag for whether the result is within tolerance before looping

  for(size_t i=0;i<xSimSize;i++)
    for(size_t j=0;j<ySimSize;j++)
      for(size_t k=0;k<zSimSize;k++){

        d[i][j][k] = (u[oldIndex][i][j][k] - u[oldIndex][i+1][j][k] + v[oldIndex][i][j][k] - v[oldIndex][i][j+1][k] + w[oldIndex][i][j][k] - w[oldIndex][i][j][k+1])/2;
        r[i][j][k] = d[i][j][k];
        if(abs(r[i][j][k])>tol)
          breakIteration = false;
        p[0][i][j][k] = 0;
  }

  if(breakIteration)
    it = maxIterations + 2; // give it a value that will skip the loop, but also lets us know that maxIterations was not truly exceeded

  //first we apply the preconditioner
  applyPreconditioner(sigma, A, z, d, r, s, precon, q);

  // now loop

  float maxAbsR, a, b;
  int it;

  for(it=0;it<maxIterations;it++){

      maxAbsR = 0;

      applyA(s,z,A);
      a = rho / dotProduct(z,s);

      for(int i=0;i<xSimSize;i++)
        for(int j=0;j<ySimSize;j++)
          for(int k=0;k<zSimSize;k++){

              p[newIndex][i][j][k] += a * s[i][j][k];
              r[i][j][k] -= a * z[i][j][k];

              if(abs(r[i][j][k]) > maxAbsR)
                maxAbsR = abs(r[i][j][k]);
            }

      if(maxAbsR <= tol)
        break;

      applyPreconditioner(sigma, A, z, d, r, s, precon, q);

      b = sigma / rho;

      for(int i=0;i<xSimSize;i++)
        for(int j=0;j<ySimSize;j++)
          for(int k=0;k<zSimSize;k++)

            s[i][j][k] = z[i][j][k] + b * s[i][j][k];

    }

  // now compute the new velocities

  for(i=0;i<xSimSize;i++)
    for(j=0;j<ySimSize;j++)
      for(k=0;k<zSimSize;k++){
        u[newIndex][i][j][k] = -dt/rho * (p[newIndex][i+1][j][k] - p[newIndex][i][j][k]) / dx + u[oldIndex][i][j][k];
        v[newIndex][i][j][k] = -dt/rho * (p[newIndex][i][j+1][k] - p[newIndex][i][j][k]) / dx + v[oldIndex][i][j][k];
        w[newIndex][i][j][k] = -dt/rho * (p[newIndex][i][j][k+1] - p[newIndex][i][j][k]) / dx + v[oldIndex][i][j][k];
      }

  if(it < maxIterations || it == maxIterations + 2)
    return true;

  return false;
}

// applies the preconditioner, also does the dotproduct for sigma, so we don't loop the whole grid again
// returns true on success, false on failure
bool OpenGLWindow::applyPreconditioner(float& sigma, std::vector<std::vector<std::vector<SevenPointLagrangianMatrixElement>>> A,
                                       std::vector<std::vector<std::vector<float>>> z, std::vector<std::vector<std::vector<float>>> d,
                                       std::vector<std::vector<std::vector<float>>> r, std::vector<std::vector<std::vector<float>>> s,
                                       std::vector<std::vector<std::vector<float>>> precon, std::vector<std::vector<std::vector<float>>> q){

  float e, t;

  for(size_t i=0;i<xSimSize;i++)
    for(size_t j=0;j<ySimSize;j++)
      for(size_t k=0;k<zSimSize;k++){

          //THIS ISN'T GOING TO WORK BECAUSE YOU SET P TO ZERO
          //if(p[i][j][k] == 0) // if there is no fluid in this cell, skip it
          //  continue; // need to set stuff to zero? (probably not)

          // apply preconditioner (is the i,j or k = 0 limit behaviour correct?)

          e = A[i][j][k].diag;
          if(i>0)
            e-= pow((A[i-1][j][k].iUp * r[i-1][j][k]),2) + tau * (A[i-1][j][k].iUp * (A[i-1][j][k].jUp + A[i-1][j][k].kUp)) * pow(precon[i-1][j][k],2);
          if(j>0)
            e-= pow((A[i][j-1][k].jUp * r[i][j-1][k]),2) + tau * (A[i][j-1][k].jUp * (A[i][j-1][k].iUp + A[i][j-1][k].kUp)) * pow(precon[i][j-1][k],2);
          if(k>0)
            e-= pow((A[i][j][k-1].iUp * r[i][j][k-1]),2) + tau * (A[i][j][k-1].iUp * (A[i][j][k-1].jUp + A[i][j][k-1].kUp)) * pow(precon[i][j][k-1],2);

          precon[i][j][k] = 1.0/sqrt(e+1E-30);

          t = r[i][j][k];
          if(i>0)
            t-= A[i-1][j][k].iUp * precon[i-1][j][k] * q[i-1][j][k];
          if(j>0)
            t-= A[i][j-1][k].iUp * precon[i][j-1][k] * q[i][j-1][k];
          if(k>0)
            t-= A[i][j][k-1].iUp * precon[i][j][k-1] * q[i][j][k-1];

          q[i][j][k] = t * precon[i][j][k];
        }

  sigma = 0;

  for(int i=xSimSize-1;i>=0;i--)
    for(int j=ySimSize-1;j>=0;j--)
      for(int k=zSimSize-1;k>=0;k--){

          t = q[i][j][k];
          if(i<xSimSize)
            t-= A[i][j][k].iUp * precon[i][j][k] * z[i+1][j][k];
          if(j<ySimSize)
            t-= A[i][j][k].jUp * precon[i][j][k] * z[i][j+1][k];
          if(k<zSimSize)
            t-= A[i][j][k].kUp * precon[i][j][k] * z[i][j][k+1];

          z[i][j][k] = t * precon[i][j][k];

          // s is the search vector
          s[i][j][k] = z[i][j][k];

          // set sigma as the dot product of z and r
          sigma += z[i][j][k] * r[i][j][k];
        }
}

float OpenGLWindow::dotProduct(std::vector<std::vector<std::vector<float>>> aMatrix, std::vector<std::vector<std::vector<float>>> bMatrix){

  size_t xSize = aMatrix.size();
  size_t ySize = aMatrix[0].size();
  size_t zSize = aMatrix[0][0].size();

  if(bMatrix.size() != aMatrix.size() || bMatrix[0].size() != aMatrix[0].size() || bMatrix[0][0].size() != aMatrix[0][0].size())
    return NAN;

  float d = 0;

  for(int i=0;i<xSize;i++)
    for(int j=0;j<ySize;j++)
      for(int k=0;k<zSize;k++)

        d += aMatrix[i][j][k] * bMatrix[i][j][k];

  return d;

}

// returns true on success, false on failure
bool OpenGLWindow::applyA(std::vector<std::vector<std::vector<float>>> aMatrix, std::vector<std::vector<std::vector<float>>> targetMatrix,
                          std::vector<std::vector<std::vector<SevenPointLagrangianMatrixElement>>> A){

  size_t xSize = aMatrix.size();
  size_t ySize = aMatrix[0].size();
  size_t zSize = aMatrix[0][0].size();

  if(targetMatrix.size() != A.size() || targetMatrix[0].size() != A[0].size() || targetMatrix[0][0].size() != A[0][0].size())
    return false;

  float d = 0;

  for(int i=0;i<xSize;i++)
    for(int j=0;j<ySize;j++)
      for(int k=0;k<zSize;k++){

          targetMatrix[i][j][k] = A[i][j][k].diag * aMatrix[i][j][k];

          if(i<xSimSize)
            targetMatrix[i][j][k]+= A[i][j][k].iUp * aMatrix[i+1][j][k];
          if(j<ySimSize)
            targetMatrix[i][j][k]+= A[i][j][k].jUp * aMatrix[i][j+1][k];
          if(k<zSimSize)
            targetMatrix[i][j][k]+= A[i][j][k].kUp * aMatrix[i][j][k+1];

          // A[i][j][k][i-1][j][k] = A[i-1][j][k][i][j][k] (symmetry)

          if(i<xSimSize)
            targetMatrix[i][j][k]+= A[i-1][j][k].iUp * aMatrix[i-1][j][k];
          if(j<ySimSize)
            targetMatrix[i][j][k]+= A[i][j-1][k].jUp * aMatrix[i][j-1][k];
          if(k<zSimSize)
            targetMatrix[i][j][k]+= A[i][j][k-1].kUp * aMatrix[i][j][k-1];
        }
}

bool OpenGLWindow::advect(std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> args,
                          std::vector<bool> isCentered, std::vector<float> outsideValues)
{
  int oldIndex = 0;
  int newIndex = 1;


  // if args is empty, assume only u,v,w and p need to be updated

  if(args.size() != isCentered.size()){
      std::cout<<"args and isCentered mismatch";
      return false;
  }
  if(args.size() != outsideValues.size()){
      std::cout<<"args and outsideValues mismatch";
      return false;
  }

  int k=0;
  for(auto c:args){

      if(isCentered[k]==0){

          ngl::Vec3 velocity;
          float xn, xp, yn, yp, zn, zp, ax, ay, az;
          int ixp, iyp, izp;

          for(int ix=0;ix<xSimSize;ix++)
            for(int iy=0;iy<ySimSize;iy++)
              for(int iz=0;iz<zSimSize;iz++){

                  // get velocity components at center of current grid cell through trilinear interpolation
                  velocity.m_x = (u[oldIndex][ix][iy][iz]+u[oldIndex][ix+1][iy][iz])/2;
                  velocity.m_y = (v[oldIndex][ix][iy][iz]+v[oldIndex][ix][iy+1][iz])/2;
                  velocity.m_z = (w[oldIndex][ix][iy][iz]+w[oldIndex][ix][iy][iz+1])/2;

                  // the current grid cell (of indices (ix,iy,iz)) is viewed as a particle
                  // (xn,yn,zn) is the corresponding position of this particle; (xp,yp,zp) is its projected past position a time dt ago
                  // we use Forward Euler (might update to RK2)
                  // (ixp,iyp,izp) are the indices of the grid point (xpFloor,ypFloor,zpFloor) such that xp is in [xpFloor, xpFloor+dx), yp is in [ypFloor, ypFloor+dx), zp is in [zpFloor, zpFloor+dx)
                  xn = dx * (ix + 0.5);
                  xp = xn - dt * velocity.m_x;
                  ixp = floor(xp/dx - 0.5);
                  ax = dx * (ixp + 0.5);

                  yn = dx * (iy + 0.5);
                  yp = yn - dt * velocity.m_y;
                  iyp = floor(yp/dx - 0.5);
                  ay = dx * (iyp + 0.5);

                  zn = dx * (iz + 0.5);
                  zp = zn - dt * velocity.m_z;
                  izp = floor(zp/dx - 0.5);
                  az = dx * (izp + 0.5);

                  for(int i=0;i<2;i++)
                    for(int j=0;j<2;j++)
                      for(int k=0;k<2;k++){
                          if(ixp+i>=xSimSize || iyp+j>=ySimSize || izp+k>=zSimSize){ // outside the simulation volume, there is a constant value for the quantity
                            c[newIndex][ix][iy][iz] += outsideValues[k];
                            continue;}
                          c[newIndex][ix][iy][iz]+=(i?ax:(1-ax))*(j?ay:(1-ay))*(k?(1-az):az)*c[newIndex][ixp+i][iyp+j][izp+k];
                        }
                }
        }
      else{

          size_t xSize,ySize,zSize;
          float xn, xp, yn, yp, zn, zp, ax, ay, az;
          int ixp, iyp, izp;

          xSize = c.size();
          ySize = c[0].size();
          zSize = c[0][0].size();

          for(size_t ix=0;ix<xSize;ix++)
            for(size_t iy=0;iy<ySize;iy++)
              for(size_t iz=0;iz<zSize;iz++){

                  // don't need interpolation for xp, yp, zp, since they are at grid cell edges, just like the velocities

                  // the current grid wall (of indices (ix,iy,iz)) is viewed as a particle; the wall plane depends on the current data grid
                  // (xn,yn,zn) is the corresponding position of this particle; (xp,yp,zp) is its projected past position a time dt ago
                  // we use Forward Euler (might update to RK2)
                  // (ixp,iyp,izp) are the indices of the grid wall (xpFloor,ypFloor,zpFloor) such that xp is in [xpFloor, xpFloor+dx), yp is in [ypFloor, ypFloor+dx), zp is in [zpFloor, zpFloor+dx)
                  xn = dx * ix;
                  xp = xn - dt * u[oldIndex][ix][iy][iz];
                  ixp = floor(xp/dx);
                  ax = dx * ixp;

                  yn = dx * iy;
                  yp = yn - dt * v[oldIndex][ix][iy][iz];
                  iyp = floor(yp/dx);
                  ay = dx * iyp;

                  zn = dx * iz;
                  zp = zn - dt * w[oldIndex][ix][iy][iz];
                  izp = floor(zp/dx);
                  az = dx * izp;

                  for(int i=0;i<2;i++)
                    for(int j=0;j<2;j++)
                      for(int k=0;k<2;k++){
                          if(ixp+i>=xSize || iyp+j>=ySize || izp+k>=zSize){ // outside the simulation volume, there is a constant value for the quantity
                            c[newIndex][ix][iy][iz] += outsideValues[k];
                            continue;}
                          c[newIndex][ix][iy][iz]+=(i?ax:(1-ax))*(j?ay:(1-ay))*(k?(1-az):az)*c[newIndex][ixp+i][iyp+j][izp+k];
                        }

                }

        }
    }
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
