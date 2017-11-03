#ifndef BAKE_H_
#define BAKE_H_
#include <vector>
#include <cstdio>

//template<typename T>
//using Matrix3D = std::vector<std::vector<std::vector<T>>>;
//using Matrix4D = std::vector<std::vector<std::vector<std::vector<T>>>>;
//using TwoStepMatrix3D =

// 3D fixed-size matrix
template<typename T>
class Matrix3D<T>{
public:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief ctor for 3D matrix
  /// @param [in] xSize the x-dimension of the matrix (outermost)
  /// @param [in] ySize the y-dimension of the matrix
  /// @param [in] zSize the z-dimension of the matrix (innermost)
  //----------------------------------------------------------------------------------------------------------------------
  Matrix3D(size_t xSize, size_t ySize, size_t zSize);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief dtor for 3D matrix
  //----------------------------------------------------------------------------------------------------------------------
  ~Matrix3D();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief setter
  /// @param [in] x the x-coordinate of the element to set
  /// @param [in] y the y-coordinate of the element to set
  /// @param [in] z the z-coordinate of the element to set
  /// @param [in] value the value to set the element to
  //----------------------------------------------------------------------------------------------------------------------
  bool set(size_t x,size_t y,size_t z,T value);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief getter
  //----------------------------------------------------------------------------------------------------------------------
  T get(size_t x,size_t y,size_t z);

  T dotproduct(Matrix3D<T> bMatrix);

  bool apply7PLMatrix(SevenPointLagrangianMatrix A, Matrix3D<T> targetMatrix);

private:
  std::vector data;
  size_t _x, _y, _z;
  std::vector data;
};

template<typename T>
class TwoStepMatrix3D<T>{
public:
  Matrix3D<T> newM;
  Matrix3D<T> oldM;

}

class GridTuple{

}

class GridsHolder{
public:
  GridsHolder(...);

}
