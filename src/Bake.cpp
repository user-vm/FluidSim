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

bool Matrix3D::apply7PLMatrix(SevenPointLagrangianMatrix A, Matrix3D targetMatrix);

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

GridsHolder::GridsHolder(std::vector<GridTuple> listOfGrids){

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
