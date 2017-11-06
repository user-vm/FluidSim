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

TwoStepMatrix3D::TwoStepMatrix3D(size_t x_Size, size_t y_Size, size_t z_Size){

  std::unique_ptr<Matrix3D> newPointer(new Matrix3D(x_Size, y_Size, z_Size));
  newM.swap(newPointer);

  std::unique_ptr<Matrix3D> oldPointer(new Matrix3D(x_Size, y_Size, z_Size));
  oldM.swap(oldPointer);
}

TwoStepMatrix3D::~TwoStepMatrix3D(){

  newM.reset();
  oldM.reset();
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

  for(auto i: listOfGrids){

      switch(i._gridType){

        case GRID_3D:         _grids.push_back(Matrix3D(i._x,i._y,i._z));
        case GRID_3D_7PL:     _grids.push_back(SevenPointLagrangianMatrix(i._x,i._y,i._z));
        case GRID_3D_TWOSTEP: {
            TwoStepMatrix3D X = {};
            _grids.push_back(TwoStepMatrix3D(i._x,i._y,i._z));
        default:              continue;
        }

      _gridNames.push_back(i._gridName);
      _gridTypes.push_back(i._gridType);

    }
}

size_t GridsHolder::size(){
  return _gridNames.size();
}
