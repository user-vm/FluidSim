#include "Bake.h"
#include <cmath>

Matrix3D::Matrix3D(size_t xSize, size_t ySize, size_t zSize){

  _xSize = xSize;
  _ySize = ySize;
  _zSize = zSize;

  data = std::vector<float>();
  data.resize(_xSize*_ySize*_zSize);
}

Matrix3D::Matrix3D(){

  // finish this
}

bool Matrix3D::set(size_t x,size_t y,size_t z,float value){

  if (x >= xSize || y >= ySize || z >= zSize)
    return false;

  data[_zSize * ((_ySize * x) + y) + z] = value;

  return true;
}

float Matrix3D::get(size_t x,size_t y,size_t z){

  if (x >= xSize || y >= ySize || z >= zSize)
    return NAN;

  return data[_zSize * ((_ySize * x) + y) + z];
}

float Matrix3D::dotproduct(Matrix3D bMatrix){

  float result = 0.0;

  if (this->_xSize != bMatrix._xSize || this->_ySize != bMatrix._ySize || this->_zSize != bMatrix._zSize)
    return NAN;

  for(int i;i<= this->data.size();i++)
    result += this->data[i] * bMatrix.data[i];

  return result;
}

float Matrix3D::apply7PLMatrix(SevenPointLagrangianMatrix A, Matrix3D targetMatrix){


}

bool Matrix3D::apply7PLMatrix(SevenPointLagrangianMatrix A, Matrix3D targetMatrix);
