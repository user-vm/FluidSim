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

#define TAU_TUNING_CONSTANT 0.97

enum GridType {GRID_3D = 1, GRID_3D_TWOSTEP = 2, GRID_3D_7PL = 3, GRID_INVALID = 0};
/*
class DisplayableMatrix{

public:
  size_t xSize() = 0;
  size_t ySize() = 0;
  size_t zSize() = 0;

  float get(size_t x, size_t y, size_t z);
};
*/
// a representation of a 6D seven point Lagrangian matrix, containing only four values due to symmetry
class SevenPointLagrangianMatrix{
public:
  SevenPointLagrangianMatrix(size_t xSize, size_t ySize, size_t zSize, bool defaultInitialization);
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

  float dotProduct(Matrix3D* bMatrix);

  bool apply7PLMatrix(SevenPointLagrangianMatrix* A, Matrix3D* targetMatrix);

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

  // default get, so the interface can be used with this class
  // gets from old
  float get(size_t x, size_t y, size_t z);

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
  GridTuple(std::string gridName, GridType gridType, size_t x, size_t y, size_t z, float outsideValue);
  //~GridTuple();

  friend class GridsHolder;

private:
  GridType _gridType;
  std::string _gridName;
  size_t _x,_y,_z;
  float _outsideValue;
};

class GridsHolder{

public:
  GridsHolder(std::vector<std::unique_ptr<GridTuple>> listOfGrids, float gridCellSize, float timeStep,
              float projectionTolerance, size_t maxIterations,
              float density, ngl::Vec3 g);
  //~GridsHolder();

  bool append(std::vector<GridTuple> listOfGrids);

  bool append(GridTuple gridTuple);

  bool append(std::unique_ptr<GridTuple>);

  bool append(std::vector<std::unique_ptr<GridTuple>> listOfGrids);

  size_t size();

  //typedef boost::variant<Matrix3D*, SevenPointLagrangianMatrix*, TwoStepMatrix3D*> GridElement;

  Matrix3D* getMatrix3DByName(std::string name);

  SevenPointLagrangianMatrix* getSevenPointLagrangianMatrixByName(std::string name);

  TwoStepMatrix3D* getTwoStepMatrix3DByName(std::string name);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief advection function
  //----------------------------------------------------------------------------------------------------------------------
  bool advect(std::vector<std::string> gridsToAdvectNames);
  bool advect(std::vector<std::string> gridsToAdvectNames, float dt);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief function to apply preconditioner  to residual "vector" (3d matrix) and set sigma to dotproduct of z and r
  //----------------------------------------------------------------------------------------------------------------------
  bool applyPreconditioner(std::string targetName, float& sigma);

  // todo: the renamed thing doesn't do anything
  bool applyPreconditioner(std::string targetName, float& sigma, std::vector<std::array<std::string,2>> renamedVariables);

  void setDefaultTimestep(float value);

  void setDefaultDx(float value);

  void setDefaultMaxIterations(size_t value);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief projection function
  //----------------------------------------------------------------------------------------------------------------------
  bool project();
  bool project(float dt);

  // just provide NAN, zero or negative values (the latter only for the first two) to use defaults
  bool project(float dt, float tol, size_t maxIterations);

  bool body();
  bool body(float dt);

  // this might not be good
  // template<class T>
  // T* getAnyByName(std::string name);

  GridType getTypeByName(std::string name);

  // placeholder functions for testing FrameData and displaying the frames
  bool advectDummy(std::vector<std::string> gridsToAdvectNames);
  bool advectDummy(std::vector<std::string> gridsToAdvectNames, float dt);

  bool bodyDummy();
  bool bodyDummy(float dt);

  bool projectDummy();
  bool projectDummy(float dt);

private:
  // union doesn't work for these classes, at least not like this
  /*
  union Matrix{
    Matrix3D _3D;
    TwoStepMatrix3D _twoStep3D;
    SevenPointLagrangianMatrix _sevenPointLagrangian;
  };*/
  //typedef std::array<std::unique_ptr<Matrix3D>,2> TwoStepMatrix3D;

  std::vector<std::string> _gridNames;
  //std::vector<GridElement> _grids;

  std::vector<Matrix3D*> _grids_M3D;
  std::vector<TwoStepMatrix3D*> _grids_2SM3D;
  std::vector<SevenPointLagrangianMatrix*> _grids_7PL;

  std::vector<GridType> _gridTypes;
  std::vector<float> _outsideValues;

  // grid cell size, default 0.01
  float dx;

  // timestep; for now changes only when next frame is less than default_dt away; default 0.01
  float default_dt;

  // default tolerance
  float default_tol;

  // default maximum iterations of projection function
  size_t default_maxIterations;

  // density (might be variable over volume or time, so this might not be enough)
  float _density;

  // gravitational acceleration
  ngl::Vec3 _g;
};

#endif
