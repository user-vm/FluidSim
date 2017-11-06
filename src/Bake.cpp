#include "Bake.h"
#include <cmath>

Matrix3D::Matrix3D(size_t xSize, size_t ySize, size_t zSize) : _xSize(xSize), _ySize(ySize), _zSize(zSize){

  _xSize = xSize;
  _ySize = ySize;
  _zSize = zSize;

  data->resize(_xSize*_ySize*_zSize);
}

Matrix3D::~Matrix3D(){

  data.reset();
}

bool Matrix3D::set(size_t x,size_t y,size_t z,float value){

  if (x >= _xSize || y >= _ySize || z >= _zSize)
    return false;

  (*(data.get()))[_zSize * ((_ySize * x) + y) + z] = value;

  return true;
}

float Matrix3D::get(size_t x,size_t y,size_t z){

  if (x >= _xSize || y >= _ySize || z >= _zSize)
    return NAN;

  return (*(data.get()))[_zSize * ((_ySize * x) + y) + z];
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

  float d = 0;

  for(int i=0;i<xSize();i++)
    for(int j=0;j<ySize();j++)
      for(int k=0;k<zSize();k++){

          targetMatrix.get(i,j,k) = A(i,j,k).diag * this->(i,j,k);

          if(i<xSimSize)
            targetMatrix.get(i,j,k)+= A.get(i,j,k).iUp * this->get(i+1,j,k);
          if(j<ySimSize)
            targetMatrix.get(i,j,k)+= A.get(i,j,k).jUp * this->get(i,j+1,k);
          if(k<zSimSize)
            targetMatrix.get(i,j,k)+= A.get(i,j,k).kUp * this->get(i,j,k+1);

          // A[i][j][k][i-1][j][k] = A[i-1][j][k][i][j][k] (symmetry)

          if(i<xSimSize)
            targetMatrix.get(i,j,k)+= A.get(i-1,j,k).iUp * this->get(i-1,j,k);
          if(j<ySimSize)
            targetMatrix.get(i,j,k)+= A.get(i,j-1,k).jUp * this->get(i,j-1,k);
          if(k<zSimSize)
            targetMatrix.get(i,j,k)+= A.get(i,j,k-1).kUp * this->get(i,j,k-1);
        }
}

bool Matrix3D::apply7PLMatrix(SevenPointLagrangianMatrix A, Matrix3D targetMatrix);

// SevenPointLagrangianMatrix functions

SevenPointLagrangianMatrix::SevenPointLagrangianMatrix(){

  xSize = xSize;
  ySize = ySize;
  zSize = zSize;

  data.resize(xSize*ySize*zSize);
}

SevenPointLagrangianMatrix::~SevenPointLagrangianMatrix(){

  data.reset();
}

size_t SevenPointLagrangianMatrix::xSize(){
  return xSize;
}

size_t SevenPointLagrangianMatrix::ySize(){
  return ySize;
}

size_t SevenPointLagrangianMatrix::zSize(){
  return zSize;
}
