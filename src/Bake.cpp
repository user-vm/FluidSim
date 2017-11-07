#include "Bake.h"

// Matrix3D functions

Matrix3D::Matrix3D(size_t xSize, size_t ySize, size_t zSize) : _xSize(xSize), _ySize(ySize), _zSize(zSize){

  _xSize = xSize;
  _ySize = ySize;
  _zSize = zSize;

  data = std::vector<float>();
  data.resize(_xSize*_ySize*_zSize);
}

Matrix3D::~Matrix3D(){

  data = std::vector<float>();
}

bool Matrix3D::set(size_t x,size_t y,size_t z,float value){

  if (x >= _xSize || y >= _ySize || z >= _zSize)
    return false;

  data[_zSize * ((_ySize * x) + y) + z] = value;

  return true;
}

float Matrix3D::get(size_t x,size_t y,size_t z){

  if (x >= _xSize || y >= _ySize || z >= _zSize)
    return NAN;

  return data[_zSize * ((_ySize * x) + y) + z];
}

size_t Matrix3D::xSize(){
  return _xSize;
}

size_t Matrix3D::ySize(){
  return _ySize;
}

size_t Matrix3D::zSize(){
  return _zSize;
}

float Matrix3D::dotproduct(Matrix3D bMatrix){

  float result = 0.0;

  if (this->_xSize != bMatrix._xSize || this->_ySize != bMatrix._ySize || this->_zSize != bMatrix._zSize)
    return NAN;

  // this dotproduct treats matrices like vectors, so we can just parse the whole std::vector objects and multiply element by element
  for(size_t i;i<= _xSize;i++)
    for(size_t j;j<=_ySize;j++)
      for(size_t k;k<=_zSize;k++)
        result += this->get(i,j,k) * bMatrix.get(i,j,k);
  return result;
}

bool Matrix3D::apply7PLMatrix(SevenPointLagrangianMatrix A, Matrix3D targetMatrix){

  if(targetMatrix.xSize() != A.xSize() || targetMatrix.ySize() != A.ySize() || targetMatrix.zSize() != A.zSize())
    return false;

  if(xSize() != A.xSize() || ySize() != A.ySize() || zSize() != A.zSize())
    return false;

  for(size_t i=0;i<_xSize;i++)
    for(size_t j=0;j<_ySize;j++)
      for(size_t k=0;k<_zSize;k++){

          targetMatrix.set(i,j,k,A.get(i,j,k).diag * this->get(i,j,k));

          if(i<_xSize-1)
            targetMatrix.set(i,j,k,A.get(i,j,k).iUp * this->get(i+1,j,k));
          if(j<_ySize-1)
            targetMatrix.set(i,j,k,A.get(i,j,k).jUp * this->get(i,j+1,k));
          if(k<_zSize-1)
            targetMatrix.set(i,j,k,A.get(i,j,k).kUp * this->get(i,j,k+1));

          // A[i][j][k][i-1][j][k] = A[i-1][j][k][i][j][k] (symmetry)

          if(i>0)
            targetMatrix.set(i,j,k,A.get(i-1,j,k).iUp * this->get(i-1,j,k));
          if(j>0)
            targetMatrix.set(i,j,k,A.get(i,j-1,k).jUp * this->get(i,j-1,k));
          if(k>0)
            targetMatrix.set(i,j,k,A.get(i,j,k-1).kUp * this->get(i,j,k-1));
        }

  return true;
}

bool Matrix3D::isCentered(){
  return _isCentered;
}

//bool Matrix3D::apply7PLMatrix(SevenPointLagrangianMatrix A, Matrix3D targetMatrix);

// SevenPointLagrangianMatrix functions

SevenPointLagrangianMatrix::SevenPointLagrangianMatrix(size_t xSize, size_t ySize, size_t zSize){

  _xSize = xSize;
  _ySize = ySize;
  _zSize = zSize;

  data = std::vector<SevenPointLagrangianMatrixElement>();
  data.resize(_xSize*_ySize*_zSize);
}

SevenPointLagrangianMatrix::~SevenPointLagrangianMatrix(){

  data = std::vector<SevenPointLagrangianMatrixElement>();
}

size_t SevenPointLagrangianMatrix::xSize(){
  return _xSize;
}

size_t SevenPointLagrangianMatrix::ySize(){
  return _ySize;
}

size_t SevenPointLagrangianMatrix::zSize(){
  return _zSize;
}

SevenPointLagrangianMatrix::SevenPointLagrangianMatrixElement SevenPointLagrangianMatrix::get(size_t x,size_t  y, size_t z){

  if (x >= _xSize || y >= _ySize || z >= _zSize)
    return NAN_7PLM_ELEMENT;

  return data[_zSize * ((_ySize * x) + y) + z];
}

// TwoStepMatrix3D functions

TwoStepMatrix3D::TwoStepMatrix3D(size_t x_Size, size_t y_Size, size_t z_Size){

  newM = new Matrix3D(x_Size, y_Size, z_Size);
  oldM = new Matrix3D(x_Size, y_Size, z_Size);
}

TwoStepMatrix3D::~TwoStepMatrix3D(){

  delete newM;
  delete oldM;
}

float TwoStepMatrix3D::getNew(size_t x, size_t y, size_t z){

  return newM->get(x,y,z);
}

float TwoStepMatrix3D::getOld(size_t x, size_t y, size_t z){

  return oldM->get(x,y,z);
}

bool TwoStepMatrix3D::setNew(size_t x, size_t y, size_t z, float value){

  return newM->set(x,y,z,value);
}

bool TwoStepMatrix3D::setOld(size_t x, size_t y, size_t z, float value){

  return oldM->set(x,y,z,value);
}

size_t TwoStepMatrix3D::xSize(){

  return newM->xSize();
}

size_t TwoStepMatrix3D::ySize(){

  return newM->ySize();
}

size_t TwoStepMatrix3D::zSize(){

  return newM->zSize();
}

void TwoStepMatrix3D::swap(){

  std::swap(newM, oldM);
}

bool TwoStepMatrix3D::isCentered(){

  return newM->isCentered();
}

bool TwoStepMatrix3D::setOutsideValue(float value){

  bool newWorked = newM->setOutsideValue(value);
  bool oldWorked = oldM->setOutsideValue(value);

  return (newWorked && oldWorked);
}

float TwoStepMatrix3D::getOutsideValue(){

  return newM->getOutsideValue();
}

// GridTuple functions

GridTuple::GridTuple(std::string gridName, GridType gridType, size_t x, size_t y, size_t z){

  _gridName = gridName;
  _gridType = gridType;
  _gridName = gridName;
  _x = x;
  _y = y;
  _z = z;
}

GridTuple::~GridTuple(){}

// GridsHolder functions

GridsHolder::GridsHolder(std::vector<GridTuple> listOfGrids, float gridCellSize){

  // custom lambda compare function; list gets sorted alphabetically
  std::sort(listOfGrids.begin(),listOfGrids.end(), [](GridTuple a, GridTuple b){
    return a._gridName < b._gridName;} );

  for(auto i: listOfGrids){

      // if a duplicate element name is hit, skip it
      if(!_gridNames.empty() && i._gridName == _gridNames.back()){
          std::cout<<"Skipping grid \""<<i._gridName<<"\" of type "<<(i._gridType?((i._gridType==1)\
                     ?"TwoStepMatrix3D":"SevenPointLagrangianMatrix"):"Matrix3D")<<" with duplicate name\n";
          continue;}

      switch(i._gridType){

        case GRID_3D:         _grids.push_back(GridElement(new Matrix3D(i._x,i._y,i._z)));
        case GRID_3D_7PL:     _grids.push_back(GridElement(new SevenPointLagrangianMatrix(i._x,i._y,i._z)));
        case GRID_3D_TWOSTEP: _grids.push_back(GridElement(new TwoStepMatrix3D(i._x,i._y,i._z)));
        default:              continue;
        }

      _gridNames.push_back(i._gridName);
      _gridTypes.push_back(i._gridType);

    }

  dx = gridCellSize;

}

size_t GridsHolder::size(){
  return _gridNames.size();
}

Matrix3D* GridsHolder::getMatrix3DByName(std::string name){
  for(size_t i=0;i<_gridTypes.size();i++)
    if(_gridTypes[i] == GRID_3D && _gridNames[i]==name)
      return boost::get<Matrix3D*>(_grids[i]);

  return NULL;
}

SevenPointLagrangianMatrix* GridsHolder::getSevenPointLagrangianMatrixByName(std::string name){
  for(size_t i=0;i<_gridTypes.size();i++)
    if(_gridTypes[i] == GRID_3D_7PL && _gridNames[i]==name)
      return boost::get<SevenPointLagrangianMatrix*>(_grids[i]);

  return NULL;
}

TwoStepMatrix3D* GridsHolder::getTwoStepMatrix3DByName(std::string name){
  for(size_t i=0;i<_gridTypes.size();i++)
    if(_gridTypes[i] == GRID_3D_TWOSTEP && _gridNames[i]==name)
      return boost::get<TwoStepMatrix3D*>(_grids[i]);

  return NULL;
}

template<typename T>
T* GridsHolder::getAnyByName(std::string name){
  for(size_t i=0;i<_gridTypes.size();i++)
    if(_gridNames[i]==name)
      return boost::get<T*>(_grids[i]);

  return NULL;
}

bool GridsHolder::advect(std::vector<std::string> gridsToAdvectNames){

  return advect(gridsToAdvectNames, default_dt);
}

bool GridsHolder::advect(std::vector<std::string> gridsToAdvectNames, float dt)
                          //std::vector<float> outsideValues)
{

  TwoStepMatrix3D* u = getTwoStepMatrix3DByName("u");
  TwoStepMatrix3D* v = getTwoStepMatrix3DByName("v");
  TwoStepMatrix3D* w = getTwoStepMatrix3DByName("w");

  TwoStepMatrix3D* c;

  size_t x_Size,y_Size,z_Size;

  bool advectedSomething = false;

  for(auto cName:gridsToAdvectNames){

      c = getTwoStepMatrix3DByName(cName);

      if(c == NULL){
          std::cout<<"Grid named "<<cName<<" is not in grid list; skipping...\n";
          continue;
        }

      advectedSomething = true;

      x_Size = c->xSize();
      y_Size = c->ySize();
      z_Size = c->zSize();

      if(!c->isCentered()){

          ngl::Vec3 velocity;
          float xn, xp, yn, yp, zn, zp, ax, ay, az;
          int ixp, iyp, izp;

          x_Size = c->xSize();
          y_Size = c->ySize();
          z_Size = c->zSize();

          for(size_t ix=0;ix<x_Size;ix++)
            for(size_t iy=0;iy<y_Size;iy++)
              for(size_t iz=0;iz<z_Size;iz++){

                  // get velocity components at center of current grid cell through trilinear interpolation
                  velocity.m_x = (u->getOld(ix,iy,iz) + u->getOld(ix+1,iy,iz))/2;
                  velocity.m_y = (v->getOld(ix,iy,iz) + v->getOld(ix,iy+1,iz))/2;
                  velocity.m_z = (w->getOld(ix,iy,iz) + w->getOld(ix,iy,iz+1))/2;

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

                  for(size_t i=0;i<2;i++)
                    for(size_t j=0;j<2;j++)
                      for(size_t k=0;k<2;k++){
                          if(ixp+i>=x_Size || iyp+j>=y_Size || izp+k>=z_Size){ // outside the simulation volume, there is a constant value for the quantity
                            c->setNew(ix,iy,iz,c->getOutsideValue());
                            continue;}
                          c->setNew(ix,iy,iz,(i?ax:(1-ax))*(j?ay:(1-ay))*(k?(1-az):az) * c->getNew(ixp+i,iyp+j,izp+k));
                        }
                }
        }
      else{

          float xn, xp, yn, yp, zn, zp, ax, ay, az;
          int ixp, iyp, izp;

          for(size_t ix=0;ix<x_Size;ix++)
            for(size_t iy=0;iy<y_Size;iy++)
              for(size_t iz=0;iz<z_Size;iz++){

                  // don't need interpolation for xp, yp, zp, since they are at grid cell edges, just like the velocities

                  // the current grid wall (of indices (ix,iy,iz)) is viewed as a particle; the wall plane depends on the current data grid
                  // (xn,yn,zn) is the corresponding position of this particle; (xp,yp,zp) is its projected past position a time dt ago
                  // we use Forward Euler (might update to RK2)
                  // (ixp,iyp,izp) are the indices of the grid wall (xpFloor,ypFloor,zpFloor) such that xp is in [xpFloor, xpFloor+dx), yp is in [ypFloor, ypFloor+dx), zp is in [zpFloor, zpFloor+dx)
                  xn = dx * ix;
                  xp = xn - dt * u->getOld(ix,iy,iz);
                  ixp = floor(xp/dx);
                  ax = dx * ixp;

                  yn = dx * iy;
                  yp = yn - dt * v->getOld(ix,iy,iz);
                  iyp = floor(yp/dx);
                  ay = dx * iyp;

                  zn = dx * iz;
                  zp = zn - dt * w->getOld(ix,iy,iz);
                  izp = floor(zp/dx);
                  az = dx * izp;

                  for(size_t i=0;i<2;i++)
                    for(size_t j=0;j<2;j++)
                      for(size_t k=0;k<2;k++){
                          if(ixp+i>=x_Size || iyp+j>=y_Size || izp+k>=z_Size){ // outside the simulation volume, there is a constant value for the quantity
                            c->setNew(ix,iy,iz,c->getOutsideValue());
                            continue;}
                          c->setNew(ix,iy,iz,(i?ax:(1-ax))*(j?ay:(1-ay))*(k?(1-az):az)*c->getNew(ixp+i,iyp+j,izp+k));
                        }

                }

        }
    }

  return advectedSomething;
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

bool GridsHolder::applyPreconditioner(std::string targetName, float& sigma){

  applyPreconditioner(targetName, sigma, NULL);
}

// applies the preconditioner, also does the dotproduct for sigma, so we don't loop the whole grid again
// returns true on success, false on failure
bool GridsHolder::applyPreconditioner(std::string targetName, float& sigma, std::vector<std::array<std::string,2>> renamedVariables){

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
