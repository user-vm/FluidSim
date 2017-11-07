#ifndef BAKE_H_
#define BAKE_H_
#include <vector>
#include <string>
#include <cstdio>
#include <memory>
#include <cmath>
#include <boost/variant/get.hpp>
#include <boost/variant.hpp>
#include <utility>
#include <iostream>
#include <ngl/Vec3.h>

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

  friend class GridsHolder;

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

  bool isCentered();

  float getOutsideValue();

  bool setOutsideValue(float value);

  friend class GridsHolder;
  friend class TwoStepMatrix3D;

private:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief ctor for 3D matrix
  /// @param [in] xSize the x-dimension of the matrix (outermost)
  /// @param [in] ySize the y-dimension of the matrix
  /// @param [in] zSize the z-dimension of the matrix (innermost)
  //----------------------------------------------------------------------------------------------------------------------
  Matrix3D(size_t xSize, size_t ySize, size_t zSize);
  size_t _xSize, _ySize, _zSize;
  std::vector<float> data;
  bool _isCentered;
  float outsideValue;
};
/*
// container for two Matrix3D, "old" and "new"
struct TwoStepMatrix3D{
  Matrix3D newM;
  Matrix3D oldM;
}
*/
class TwoStepMatrix3D{
public:

  ~TwoStepMatrix3D();

  bool setOld(size_t x, size_t y, size_t z, float value);
  bool setNew(size_t x, size_t y, size_t z, float value);
  float getNew(size_t x, size_t y, size_t z);
  float getOld(size_t x, size_t y, size_t z);

  size_t xSize();

  size_t ySize();

  size_t zSize();

  void swap();

  bool isCentered();

  bool setOutsideValue(float value);

  float getOutsideValue();

  friend class GridsHolder;

private:
  TwoStepMatrix3D(size_t x_Size, size_t y_Size, size_t z_Size);

  Matrix3D* newM;
  Matrix3D* oldM;
};

// class for creating tuples with information to add a grid to the GridHolder singleton
class GridTuple{
public:
  //GridTuple()
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
  GridsHolder(std::vector<GridTuple> listOfGrids, float gridCellSize);
  ~GridsHolder();

  size_t size();

  typedef boost::variant<Matrix3D*, SevenPointLagrangianMatrix*, TwoStepMatrix3D*> GridElement;

  Matrix3D* getMatrix3DByName(std::string name);

  SevenPointLagrangianMatrix* getSevenPointLagrangianMatrixByName(std::string name);

  TwoStepMatrix3D* getTwoStepMatrix3DByName(std::string name);

  bool advect(std::vector<std::string> gridsToAdvectNames);

  bool advect(std::vector<std::string> gridsToAdvectNames, float dt);

  bool applyPreconditioner(std::string targetName, float& sigma);

  bool applyPreconditioner(std::string targetName, float& sigma, std::vector<std::array<std::string,2>> renamedVariables);

private:
  // union doesn't work for these classes, at least not like this
  /*
  union Matrix{
    Matrix3D _3D;
    TwoStepMatrix3D _twoStep3D;
    SevenPointLagrangianMatrix _sevenPointLagrangian;
  };*/
  //typedef std::array<std::unique_ptr<Matrix3D>,2> TwoStepMatrix3D;

  // this is probably stupid, so it'll be private
  template<class T>
  T* getAnyByName(std::string name);

  GridType getTypeByName(std::string name);

  std::vector<std::string> _gridNames;
  std::vector<GridElement> _grids;
  std::vector<GridType> _gridTypes;

  // grid cell size
  float dx;

  // timestep; for now changes only when next frame is less than default_dt away
  float default_dt;
};

#endif
