#ifndef BAKE_H_
#define BAKE_H_
#include <vector>
#include <string>
#include <cstdio>
#include <memory>
#include <cmath>
#include <boost/variant.hpp>

enum GridType {GRID_3D = 1, GRID_3D_TWOSTEP = 2, GRID_3D_7PL = 3};

// a representation of a 6D seven point Lagrangian matrix, containing only four values due to symmetry
class SevenPointLagrangianMatrix{
public:
  SevenPointLagrangianMatrix(size_t xSize, size_t ySize, size_t zSize);
  ~SevenPointLagrangianMatrix();

  struct SevenPointLagrangianMatrixElement{
    float diag; // A(i,j,k)(i,j,k)
    float iUp;  // A(i,j,k)(i+1,j,k)
    float jUp;  // A(i,j,k)(i,j+1,k)
    float kUp;  // A(i,j,k)(i,j,k+1)
  };

  static constexpr SevenPointLagrangianMatrixElement NAN_7PLM_ELEMENT = {NAN, NAN, NAN, NAN};

  SevenPointLagrangianMatrixElement get(size_t x,size_t y,size_t z);

  size_t xSize();

  size_t ySize();

  size_t zSize();

private:

  std::vector<SevenPointLagrangianMatrixElement> data;
  size_t _xSize;
  size_t _ySize;
  size_t _zSize;
};

// 3D fixed-size matrix
class Matrix3D{
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
  bool set(size_t x,size_t y,size_t z,float value);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief getter
  //----------------------------------------------------------------------------------------------------------------------
  float get(size_t x,size_t y,size_t z);

  float dotproduct(Matrix3D bMatrix);

  bool apply7PLMatrix(SevenPointLagrangianMatrix A, Matrix3D targetMatrix);

  size_t xSize();

  size_t ySize();

  size_t zSize();

private:
  size_t _xSize, _ySize, _zSize;
  std::vector<float> data;
};

// container for two Matrix3D, "old" and "new"
struct TwoStepMatrix3D{
  Matrix3D newM;
  Matrix3D oldM;
}
/*
class TwoStepMatrix3D{
public:
  TwoStepMatrix3D(size_t x_Size, size_t y_Size, size_t z_Size);
  ~TwoStepMatrix3D();

private:
  std::unique_ptr<Matrix3D> newM;
  std::unique_ptr<Matrix3D> oldM;
};*/

// class for creating tuples with information to add a grid to the GridHolder singleton
class GridTuple{
public:
  GridTuple(std::string gridName, GridType gridType, size_t x, size_t y, size_t z);
  ~GridTuple();

  friend class GridsHolder;

private:
  GridType _gridType;
  std::string _gridName;
  size_t _x,_y,_z;
};

class GridsHolder{

public:
  GridsHolder(std::vector<GridTuple> listOfGrids);
  ~GridsHolder();

  size_t size();

private:
  // union doesn't work for these classes, at least not like this
  /*
  union Matrix{
    Matrix3D _3D;
    TwoStepMatrix3D _twoStep3D;
    SevenPointLagrangianMatrix _sevenPointLagrangian;
  };*/
  std::vector<std::string> _gridNames;
  std::vector<boost::variant<Matrix3D, SevenPointLagrangianMatrix, TwoStepMatrix3D>> _grids;
  std::vector<GridType> _gridTypes;
};

#endif
