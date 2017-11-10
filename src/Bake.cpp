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

float Matrix3D::dotProduct(Matrix3D *bMatrix){

  float result = 0.0;

  if (this->_xSize != bMatrix->_xSize || this->_ySize != bMatrix->_ySize || this->_zSize != bMatrix->_zSize)
    return NAN;

  // this dotproduct treats matrices like vectors, so we can just parse the whole std::vector objects and multiply element by element
  for(size_t i;i<= _xSize;i++)
    for(size_t j;j<=_ySize;j++)
      for(size_t k;k<=_zSize;k++)
        result += this->get(i,j,k) * bMatrix->get(i,j,k);
  return result;
}

bool Matrix3D::apply7PLMatrix(SevenPointLagrangianMatrix *A, Matrix3D *targetMatrix){

  if(targetMatrix->xSize() != A->xSize() || targetMatrix->ySize() != A->ySize() || targetMatrix->zSize() != A->zSize())
    return false;

  if(xSize() != A->xSize() || ySize() != A->ySize() || zSize() != A->zSize())
    return false;

  for(size_t i=0;i<_xSize;i++)
    for(size_t j=0;j<_ySize;j++)
      for(size_t k=0;k<_zSize;k++){

          targetMatrix->set(i,j,k,A->get(i,j,k).diag * this->get(i,j,k));

          if(i<_xSize-1)
            targetMatrix->set(i,j,k,A->get(i,j,k).iUp * this->get(i+1,j,k));
          if(j<_ySize-1)
            targetMatrix->set(i,j,k,A->get(i,j,k).jUp * this->get(i,j+1,k));
          if(k<_zSize-1)
            targetMatrix->set(i,j,k,A->get(i,j,k).kUp * this->get(i,j,k+1));

          // A[i][j][k][i-1][j][k] = A[i-1][j][k][i][j][k] (symmetry)

          if(i>0)
            targetMatrix->set(i,j,k,A->get(i-1,j,k).iUp * this->get(i-1,j,k));
          if(j>0)
            targetMatrix->set(i,j,k,A->get(i,j-1,k).jUp * this->get(i,j-1,k));
          if(k>0)
            targetMatrix->set(i,j,k,A->get(i,j,k-1).kUp * this->get(i,j,k-1));
        }

  return true;
}

bool Matrix3D::isCentered(){
  return _isCentered;
}

float Matrix3D::getOutsideValue(){

  return outsideValue;
}

bool Matrix3D::setOutsideValue(float value){

  outsideValue = value;
  return true;
}

//bool Matrix3D::apply7PLMatrix(SevenPointLagrangianMatrix A, Matrix3D targetMatrix);

// SevenPointLagrangianMatrix functions

SevenPointLagrangianMatrix::SevenPointLagrangianMatrix(size_t xSize, size_t ySize, size_t zSize, bool defaultInitialization){

  _xSize = xSize;
  _ySize = ySize;
  _zSize = zSize;

  data = std::vector<SevenPointLagrangianMatrixElement>();
  data.resize(_xSize*_ySize*_zSize);

  if(defaultInitialization){

      size_t pos = 0;

      for(size_t ai=0;ai<xSize;ai++)
        for(size_t aj=0;aj<ySize;aj++)
          for(size_t ak=0;ak<zSize;ak++,pos++){

              data[pos].diag = 6;

              if(ai<xSize-1)
                data[pos].iUp = -1;
              else
                data[pos].iUp = 0;

              if(aj<ySize-1)
                data[pos].jUp = -1;
              else
                data[pos].jUp = 0;

              if(ak<zSize-1)
                data[pos].kUp = -1;
              else
                data[pos].kUp = 0;

            }
    }

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

GridTuple::GridTuple(std::string gridName, GridType gridType, size_t x, size_t y, size_t z) : GridTuple(gridName, gridType, x, y, z, 0.0){}

GridTuple::GridTuple(std::string gridName, GridType gridType, size_t x, size_t y, size_t z, float outsideValue){

  _gridName = gridName;
  _gridType = gridType;
  _x = x;
  _y = y;
  _z = z;
  if(outsideValue != NAN)
    _outsideValue = outsideValue;
  else
    _outsideValue = 0.0;
}

//GridTuple::~GridTuple(){}

// GridsHolder functions

GridsHolder::GridsHolder(std::vector<GridTuple*> listOfGrids, float gridCellSize, float timeStep,
                         float projectionTolerance, size_t maxIterations, float density, ngl::Vec3 g){

  // custom lambda compare function; list gets sorted alphabetically
  std::sort(listOfGrids.begin(),listOfGrids.end(), [](GridTuple* a, GridTuple* b){
    return a->_gridName < b->_gridName;} );

  for(auto i: listOfGrids){

      // if a duplicate element name is hit, skip it
      if(!_gridNames.empty() && i->_gridName == _gridNames.back()){
          std::cout<<"Skipping grid \""<<i->_gridName<<"\" of type "<<(i->_gridType?((i->_gridType==1)\
                     ?"TwoStepMatrix3D":"SevenPointLagrangianMatrix"):"Matrix3D")<<" with duplicate name\n";
          continue;}

      switch(i->_gridType){

        case GRID_3D:         {_grids.push_back(GridElement(new Matrix3D(i->_x,i->_y,i->_z))); break;}
        case GRID_3D_7PL:     {_grids.push_back(GridElement(new SevenPointLagrangianMatrix(i->_x,i->_y,i->_z, true))); break;}
        case GRID_3D_TWOSTEP: {_grids.push_back(GridElement(new TwoStepMatrix3D(i->_x,i->_y,i->_z))); break;}
        default:              continue;
        }

      _gridNames.push_back(i->_gridName);
      _gridTypes.push_back(i->_gridType);
      _outsideValues.push_back(i->_outsideValue);

    }

  dx = gridCellSize;
  default_dt = timeStep;
  _density = density;
  default_tol = projectionTolerance;
  _g = g;
  default_maxIterations = maxIterations;

}

GridsHolder::~GridsHolder(){}

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
/*
template<typename T>
T* GridsHolder::getAnyByName(std::string name){
  for(size_t i=0;i<_gridTypes.size();i++)
    if(_gridNames[i]==name)
      return boost::get<T*>(_grids[i]);

  return NULL;
}
*/
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

void GridsHolder::setDefaultTimestep(float value){

  default_dt = value;
}

bool GridsHolder::project(){

  return project(default_dt, default_tol, default_maxIterations);
}

bool GridsHolder::project(float dt){

  return project(dt, default_tol, default_maxIterations);
}

bool GridsHolder::project(float dt, float tol, size_t maxIterations)
{
  // need A - 7PLM
  // need z, d, s, r, precon, q - M3D
  // need u, v, w, p - 2SM3D

  // old is 1 now, and new is 0

  if(dt == NAN || dt<=0.0)
    dt = default_dt;

  if(tol ==NAN || tol<=0.0)
    tol = default_tol;

  if(maxIterations == NAN)
    maxIterations = default_maxIterations;

  TwoStepMatrix3D* u = getTwoStepMatrix3DByName("u");
  TwoStepMatrix3D* v = getTwoStepMatrix3DByName("v");
  TwoStepMatrix3D* w = getTwoStepMatrix3DByName("w");
  TwoStepMatrix3D* p = getTwoStepMatrix3DByName("p");

  bool objectMissing = false;

  if(u==NULL){
    std::cout<<"Matrix \"u\" is missing from grid holder\n";
    objectMissing = true;}

  if(v==NULL){
    std::cout<<"Matrix \"v\" is missing from grid holder\n";
    objectMissing = true;}

  if(w==NULL){
    std::cout<<"Matrix \"w\" is missing from grid holder\n";
    objectMissing = true;}

  if(p==NULL){
    std::cout<<"Matrix \"p\" is missing from grid holder\n";
    objectMissing = true;}

  if(objectMissing)
    return false;

  Matrix3D* z = getMatrix3DByName("z");
  Matrix3D* d = getMatrix3DByName("d");
  Matrix3D* s = getMatrix3DByName("s");
  Matrix3D* r = getMatrix3DByName("r");
  Matrix3D* precon = getMatrix3DByName("precon");
  Matrix3D* q = getMatrix3DByName("q");

  SevenPointLagrangianMatrix* A = getSevenPointLagrangianMatrixByName("A");

  std::vector<std::pair<Matrix3D*,std::string>> utilityMatrixList = {std::make_pair(z,"z"),
                                                                     std::make_pair(d,"d"),
                                                                     std::make_pair(s,"s"),
                                                                     std::make_pair(r,"r"),
                                                                     std::make_pair(precon,"precon"),
                                                                     std::make_pair(q,"q")};

  for(auto i: utilityMatrixList)
      if(i.first == NULL){
          std::cout<<"Matrix3D "<<i.second<<" not in grid holder; building it...";
          append(std::unique_ptr<GridTuple>(new GridTuple(i.second,GRID_3D,p->xSize(),p->ySize(),p->zSize())));
        }

  if(A == NULL){
      std::cout<<"SevenPointLagrangianMatrix A not in grid holder; building it...";
    }

  // advect's new (rewritten by body) is project's old, so we need to swap the old and new components of the
  // TwoStepMatrix3D grids
  u->swap();
  v->swap();
  w->swap();
  p->swap();

  float sigma;

  bool breakIteration = true; // flag for whether the result is within tolerance before looping

  size_t xSize;
  size_t ySize;
  size_t zSize;

  xSize = p->xSize();
  ySize = p->ySize();
  zSize = p->zSize();

  for(size_t i=0;i<xSize;i++)
    for(size_t j=0;j<ySize;j++)
      for(size_t k=0;k<zSize;k++){

        d->set(i, j, k, (u->getOld(i,j,k) - u->getOld(i+1,j,k) + v->getOld(i,j,k) - v->getOld(i,j+1,k) + w->getOld(i,j,k) - w->getOld(i,j,k+1)/2));
        r->set(i,j,k, d->get(i,j,k));
        if(abs(r->get(i,j,k))>tol)
          breakIteration = false;
        p->setNew(i,j,k, 0);
  }

  size_t it;

  if(breakIteration)
    it = maxIterations + 2; // give it a value that will skip the loop, but also lets us know that maxIterations was not truly exceeded

  //first we apply the preconditioner
  applyPreconditioner("r",sigma);

  // now loop

  float maxAbsR, a, b;

  for(it=0;it<maxIterations;it++){

      maxAbsR = 0;

      s->apply7PLMatrix(A,z);
      a = _density / s->dotProduct(z);

      for(size_t i=0;i<xSize;i++)
        for(size_t j=0;j<ySize;j++)
          for(size_t k=0;k<zSize;k++){

              p->setNew(i,j,k, p->getNew(i,j,k) + a * s->get(i,j,k));
              r->set(i,j,k,r->get(i,j,k) - a * z->get(i,j,k));

              if(abs(r->get(i,j,k)) > maxAbsR)
                maxAbsR = abs(r->get(i,j,k));
            }

      if(maxAbsR <= tol)
        break;

      applyPreconditioner("r",sigma);

      b = sigma / _density;

      for(size_t i=0;i<xSize;i++)
        for(size_t j=0;j<ySize;j++)
          for(size_t k=0;k<zSize;k++)

            s->set(i,j,k, z->get(i,j,k) + b * s->get(i,j,k));

    }

  // now compute the new velocities

  for(size_t i=0;i<xSize;i++)
    for(size_t j=0;j<ySize;j++)
      for(size_t k=0;k<zSize;k++){
        u->setNew(i,j,k, -dt/_density * (p->getNew(i+1,j,k) - p->getNew(i,j,k)) / dx + u->getOld(i,j,k));
        v->setNew(i,j,k, -dt/_density * (p->getNew(i,j+1,k) - p->getNew(i,j,k)) / dx + v->getOld(i,j,k));
        w->setNew(i,j,k, -dt/_density * (p->getNew(i,j,k+1) - p->getNew(i,j,k)) / dx + v->getOld(i,j,k));
      }

  // swap grids back to prepare for next iteration; this means frames are drawn using the old components
  u->swap();
  v->swap();
  w->swap();
  p->swap();

  if(it < maxIterations || it == maxIterations + 2)
    return true;

  return false;
}

bool GridsHolder::applyPreconditioner(std::string targetName, float& sigma){

  std::vector<std::array<std::string,2>> renamedVariables;
  return applyPreconditioner(targetName, sigma, renamedVariables);
}

// applies the preconditioner, also does the dotproduct for sigma, so we don't loop the whole grid again
// returns true on success, false on failure
bool GridsHolder::applyPreconditioner(std::string targetName, float& sigma, std::vector<std::array<std::string,2>> renamedVariables){

  float e, t;

  float tau = TAU_TUNING_CONSTANT;

  //these are useless in this function
  //TwoStepMatrix3D* u = getTwoStepMatrix3DByName("u");
  //TwoStepMatrix3D* v = getTwoStepMatrix3DByName("v");
  //TwoStepMatrix3D* w = getTwoStepMatrix3DByName("w");
  //TwoStepMatrix3D* p = getTwoStepMatrix3DByName("p");

  // this is for checking matrix availability; skip to save time
  /*
  bool objectMissing = false;

  if(u==NULL){
    std::cout<<"Matrix \"u\" is missing from grid holder\n";
    objectMissing = true;}

  if(v==NULL){
    std::cout<<"Matrix \"v\" is missing from grid holder\n";
    objectMissing = true;}

  if(w==NULL){
    std::cout<<"Matrix \"w\" is missing from grid holder\n";
    objectMissing = true;}

  if(p==NULL){
    std::cout<<"Matrix \"p\" is missing from grid holder\n";
    objectMissing = true;}

  if(objectMissing)
    return false;
  */

  Matrix3D* z = getMatrix3DByName("z");
  //Matrix3D* d = getMatrix3DByName("d");
  Matrix3D* s = getMatrix3DByName("s");
  Matrix3D* r = getMatrix3DByName("r");
  Matrix3D* precon = getMatrix3DByName("precon");
  Matrix3D* q = getMatrix3DByName("q");

  SevenPointLagrangianMatrix* A = getSevenPointLagrangianMatrixByName("A");

  // this is for checking matrix availability; skip to save time
  /*
  std::vector<std::pair<Matrix3D*,std::string>> utilityMatrixList = {std::make_pair(z,"z"),
                                                                     std::make_pair(d,"d"),
                                                                     std::make_pair(s,"s"),
                                                                     std::make_pair(r,"r"),
                                                                     std::make_pair(precon,"precon"),
                                                                     std::make_pair(q,"q")};

  for(auto i: utilityMatrixList)
      if(i.first == NULL){
          std::cout<<"Matrix3D "<<i.second<<" not in grid holder; building it...";
          append(std::unique_ptr<GridTuple>(new GridTuple(i.second,GRID_3D,p->xSize(),p->ySize(),p->zSize())));
        }

  if(A == NULL){
      std::cout<<"SevenPointLagrangianMatrix A not in grid holder; building it...";
    }
  */

  size_t xSize = A->xSize();
  size_t ySize = A->ySize();
  size_t zSize = A->zSize();

  for(size_t i=0;i<xSize;i++)
    for(size_t j=0;j<ySize;j++)
      for(size_t k=0;k<zSize;k++){

          //THIS ISN'T GOING TO WORK BECAUSE YOU SET P TO ZERO
          //if(p[i][j][k] == 0) // if there is no fluid in this cell, skip it
          //  continue; // need to set stuff to zero? (probably not)

          // apply preconditioner (is the i,j or k = 0 limit behaviour correct?)

          e = A->get(i,j,k).diag;
          if(i>0)
            e-= pow((A->get(i-1,j,k).iUp * r->get(i-1,j,k)),2) + tau * (A->get(i-1,j,k).iUp * (A->get(i-1,j,k).jUp + A->get(i-1,j,k).kUp)) * pow(precon->get(i-1,j,k),2);
          if(j>0)
            e-= pow((A->get(i,j-1,k).jUp * r->get(i,j-1,k)),2) + tau * (A->get(i,j-1,k).jUp * (A->get(i,j-1,k).iUp + A->get(i,j-1,k).kUp)) * pow(precon->get(i,j-1,k),2);
          if(k>0)
            e-= pow((A->get(i,j,k-1).iUp * r->get(i,j,k-1)),2) + tau * (A->get(i,j,k-1).iUp * (A->get(i,j,k-1).jUp + A->get(i,j,k-1).kUp)) * pow(precon->get(i,j,k-1),2);

          precon->set(i,j,k, 1.0/sqrt(e+1E-30)); // small offset to protect from divide by zero

          t = r->get(i,j,k);
          if(i>0)
            t-= A->get(i-1,j,k).iUp * precon->get(i-1,j,k) * q->get(i-1,j,k);
          if(j>0)
            t-= A->get(i,j-1,k).iUp * precon->get(i,j-1,k) * q->get(i,j-1,k);
          if(k>0)
            t-= A->get(i,j,k-1).iUp * precon->get(i,j,k-1) * q->get(i,j,k-1);

          q->set(i,j,k, t * precon->get(i,j,k));
        }

  sigma = 0;

  for(size_t i=xSize;i>0;--i)
    for(size_t j=ySize;j>0;--j)
      for(size_t k=zSize;k>0;--k){

          t = q->get(i,j,k);
          if(i<xSize)
            t-= A->get(i,j,k).iUp * precon->get(i,j,k) * z->get(i+1,j,k);
          if(j<ySize)
            t-= A->get(i,j,k).jUp * precon->get(i,j,k) * z->get(i,j+1,k);
          if(k<zSize)
            t-= A->get(i,j,k).kUp * precon->get(i,j,k) * z->get(i,j,k+1);

          z->set(i,j,k, t * precon->get(i,j,k));

          // s is the search vector
          s->set(i,j,k, z->get(i,j,k));

          // set sigma as the dot product of z and r
          sigma += z->get(i,j,k) * r->get(i,j,k);
        }

  return true;
}

bool GridsHolder::body(){

  return body(default_dt);
}

bool GridsHolder::body(float dt){

  TwoStepMatrix3D* u = getTwoStepMatrix3DByName("u");
  TwoStepMatrix3D* v = getTwoStepMatrix3DByName("v");
  TwoStepMatrix3D* w = getTwoStepMatrix3DByName("w");
  //TwoStepMatrix3D* p = getTwoStepMatrix3DByName("p");

  // this is for checking matrix availability; skip to dave time
  /*
  bool objectMissing = false;

  if(u==NULL){
    std::cout<<"Matrix \"u\" is missing from grid holder\n";
    objectMissing = true;}

  if(v==NULL){
    std::cout<<"Matrix \"v\" is missing from grid holder\n";
    objectMissing = true;}

  if(w==NULL){
    std::cout<<"Matrix \"w\" is missing from grid holder\n";
    objectMissing = true;}

  if(p==NULL){
    std::cout<<"Matrix \"p\" is missing from grid holder\n";
    objectMissing = true;}

  if(objectMissing)
    return false;
  */

  size_t xSize = u->xSize() - 1;
  size_t ySize = u->ySize();
  size_t zSize = u->zSize();

  for(size_t i=0;i<=xSize;i++)
    for(size_t j=0;j<=ySize;j++)
      for(size_t k=0;k<=zSize;k++){

          if(i==xSize){
              u->setNew(i,j,k, u->getNew(i,j,k) + _g.m_x * dt);
              continue;}
          else
            if(j==ySize){
                v->setNew(i,j,k, v->getNew(i,j,k) + _g.m_y * dt);
                continue;}
            else
              if(k==zSize){
                  w->setNew(i,j,k, w->getNew(i,j,k) + _g.m_z * dt);
                  continue;}

          u->setNew(i,j,k, u->getNew(i,j,k) + _g.m_x * dt);
          v->setNew(i,j,k, v->getNew(i,j,k) + _g.m_y * dt);
          w->setNew(i,j,k, w->getNew(i,j,k) + _g.m_z * dt);
        }

  return true;
}

bool GridsHolder::advectDummy(std::vector<std::string> gridsToAdvectNames){

  return advectDummy(gridsToAdvectNames, default_dt);
}

bool GridsHolder::advectDummy(std::vector<std::string> gridsToAdvectNames, float dt){

  // do some stuff that modifies the grid values

  TwoStepMatrix3D* theGrid;
  bool updatedAtLeastOne = false;

  for(auto gridName:gridsToAdvectNames){

      theGrid = getTwoStepMatrix3DByName(gridName);

      if(theGrid == NULL){
          std::cout<<"Grid \""<<gridName<<"\" of type TwoStepMatrix3D not in GridsHolder object.\n";
          continue;}
      else
        updatedAtLeastOne = true;

      for(size_t i=0;i<theGrid->xSize();i++)
        for(size_t j=0;j<theGrid->ySize();j++)
          for(size_t k=0;k<theGrid->zSize();k++)

            // the p and u distributions move in the x directions, v moves in the y direction, w moves in the z direction
            theGrid->setNew(i,j,k,theGrid->getOld(((gridName=="p"||gridName=="u")?((i+1) % theGrid->xSize()):i),
                                                  ((gridName=="v")?((j+1) % theGrid->ySize()):j),
                                                  ((gridName=="w")?((k+1) % theGrid->zSize()):k)));
    }

  return updatedAtLeastOne;
}

bool GridsHolder::bodyDummy(){

  return bodyDummy(default_dt);
}

bool GridsHolder::bodyDummy(float dt){

  return true;
}

bool GridsHolder::projectDummy(){

  return projectDummy(default_dt);
}

bool GridsHolder::projectDummy(float dt){

  TwoStepMatrix3D* grid;

  for(auto i:{"u","v","w","p"}){
      grid = getTwoStepMatrix3DByName(i);
      grid->swap();
    }
  return true;
}

GridType GridsHolder::getTypeByName(std::string name){

  for(size_t i=0;i<_gridTypes.size();i++)
    if(_gridNames[i]==name)
      return _gridTypes[i];

  return GRID_INVALID;
}

bool GridsHolder::append(std::vector<GridTuple> listOfGrids){

  // custom lambda compare function; list gets sorted alphabetically
  std::sort(listOfGrids.begin(),listOfGrids.end(), [](GridTuple a, GridTuple b){
    return a._gridName < b._gridName;} );

  bool addedAtLeastOne = false;

  for(auto i: listOfGrids){

      // if a duplicate element name is hit, skip it
      if(!_gridNames.empty() && i._gridName == _gridNames.back()){
          std::cout<<"Skipping grid \""<<i._gridName<<"\" of type "<<(i._gridType?((i._gridType==1)\
                     ?"TwoStepMatrix3D":"SevenPointLagrangianMatrix"):"Matrix3D")<<" with duplicate name\n";
          continue;}

      switch(i._gridType){

        case GRID_3D:         {_grids.push_back(GridElement(new Matrix3D(i._x,i._y,i._z))); break;}
        case GRID_3D_7PL:     {_grids.push_back(GridElement(new SevenPointLagrangianMatrix(i._x,i._y,i._z, true))); break;}
        case GRID_3D_TWOSTEP: {_grids.push_back(GridElement(new TwoStepMatrix3D(i._x,i._y,i._z))); break;}
        default:              continue;
        }

      _gridNames.push_back(i._gridName);
      _gridTypes.push_back(i._gridType);
      _outsideValues.push_back(i._outsideValue);
      addedAtLeastOne = true;

    }

  return addedAtLeastOne;
}

bool GridsHolder::append(GridTuple gridTuple){

  if(!_gridNames.empty() && gridTuple._gridName == _gridNames.back())
      std::cout<<"Skipping grid \""<<gridTuple._gridName<<"\" of type "<<(gridTuple._gridType?((gridTuple._gridType==1)\
                 ?"TwoStepMatrix3D":"SevenPointLagrangianMatrix"):"Matrix3D")<<" with duplicate name\n";

  switch(gridTuple._gridType){

    case GRID_3D:         {_grids.push_back(GridElement(new Matrix3D(gridTuple._x,gridTuple._y,gridTuple._z))); break;}
    case GRID_3D_7PL:     {_grids.push_back(GridElement(new SevenPointLagrangianMatrix(gridTuple._x,gridTuple._y,gridTuple._z, true))); break;}
    case GRID_3D_TWOSTEP: {_grids.push_back(GridElement(new TwoStepMatrix3D(gridTuple._x,gridTuple._y,gridTuple._z))); break;}
    default:              return false;
    }

  _gridNames.push_back(gridTuple._gridName);
  _gridTypes.push_back(gridTuple._gridType);
  _outsideValues.push_back(gridTuple._outsideValue);

}

bool GridsHolder::append(std::vector<std::unique_ptr<GridTuple>> listOfGrids){

  // custom lambda compare function; list gets sorted alphabetically
  std::sort(listOfGrids.begin()->get(),listOfGrids.end()->get(), [](GridTuple a, GridTuple b){
    return a._gridName < b._gridName;} );

  bool addedAtLeastOne = false;

  for(size_t i=0;i<listOfGrids.size();i++){

      // if a duplicate element name is hit, skip it
      if(!_gridNames.empty() && listOfGrids[i].get()->_gridName == _gridNames.back()){
          std::cout<<"Skipping grid \""<<listOfGrids[i].get()->_gridName<<"\" of type "<<(listOfGrids[i].get()->_gridType?((listOfGrids[i].get()->_gridType==1)\
                     ?"TwoStepMatrix3D":"SevenPointLagrangianMatrix"):"Matrix3D")<<" with duplicate name\n";
          continue;}

      switch(listOfGrids[i].get()->_gridType){

        case GRID_3D:         {_grids.push_back(GridElement(new Matrix3D(listOfGrids[i].get()->_x,listOfGrids[i].get()->_y,listOfGrids[i].get()->_z))); break;}
        case GRID_3D_7PL:     {_grids.push_back(GridElement(new SevenPointLagrangianMatrix(listOfGrids[i].get()->_x,listOfGrids[i].get()->_y,listOfGrids[i].get()->_z, true))); break;}
        case GRID_3D_TWOSTEP: {_grids.push_back(GridElement(new TwoStepMatrix3D(listOfGrids[i].get()->_x,listOfGrids[i].get()->_y,listOfGrids[i].get()->_z))); break;}
        default:              continue;
        }

      _gridNames.push_back(listOfGrids[i].get()->_gridName);
      _gridTypes.push_back(listOfGrids[i].get()->_gridType);
      _outsideValues.push_back(listOfGrids[i].get()->_outsideValue);
      addedAtLeastOne = true;

    }

  return addedAtLeastOne;
}

bool GridsHolder::append(std::unique_ptr<GridTuple> gridTuple){

  if(!_gridNames.empty() && gridTuple.get()->_gridName == _gridNames.back())
      std::cout<<"Skipping grid \""<<gridTuple.get()->_gridName<<"\" of type "<<(gridTuple.get()->_gridType?((gridTuple.get()->_gridType==1)\
                 ?"TwoStepMatrix3D":"SevenPointLagrangianMatrix"):"Matrix3D")<<" with duplicate name\n";

  switch(gridTuple.get()->_gridType){

    case GRID_3D:         {_grids.push_back(GridElement(new Matrix3D(gridTuple.get()->_x,gridTuple.get()->_y,gridTuple.get()->_z))); break;}
    case GRID_3D_7PL:     {_grids.push_back(GridElement(new SevenPointLagrangianMatrix(gridTuple.get()->_x,gridTuple.get()->_y,gridTuple.get()->_z, true))); break;}
    case GRID_3D_TWOSTEP: {_grids.push_back(GridElement(new TwoStepMatrix3D(gridTuple.get()->_x,gridTuple.get()->_y,gridTuple.get()->_z))); break;}
    default:              return false;
    }

  _gridNames.push_back(gridTuple.get()->_gridName);
  _gridTypes.push_back(gridTuple.get()->_gridType);
  _outsideValues.push_back(gridTuple.get()->_outsideValue);

  return true;
}
