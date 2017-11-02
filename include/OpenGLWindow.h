#ifndef OPENGLWINDOW_H_
#define OPENGLWINDOW_H_
#include <GL/glew.h>
#include <QOpenGLWindow>
#include <QElapsedTimer>
#include <ngl/Vec3.h>

class OpenGLWindow : public QOpenGLWindow
{
    // need to tell Qt to run the MOC
    Q_OBJECT
  public:
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief ctor for OpenGL window must set the surface type to OpenGL
    /// @param [in] parent the parent window to the class
    //----------------------------------------------------------------------------------------------------------------------
    explicit OpenGLWindow();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief dtor, remember to remove the device once finished
    //----------------------------------------------------------------------------------------------------------------------

    ~OpenGLWindow();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief render method called every update
    //----------------------------------------------------------------------------------------------------------------------
    void paintGL();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief pure virtual initialize method we override in our base class to do our drawing
    /// this is only called one time, just after we have a valid GL context use this to init any global GL elements
    //----------------------------------------------------------------------------------------------------------------------
    void initializeGL();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the struct defining the elements of the coefficient matrix for use in the projection function
    //----------------------------------------------------------------------------------------------------------------------
    struct SevenPointLagrangianMatrixElement{
      float diag; // A(i,j,k)(i,j,k)
      float iUp;  // A(i,j,k)(i+1,j,k)
      float jUp;  // A(i,j,k)(i,j+1,k)
      float kUp;  // A(i,j,k)(i,j,k+1)
    };
  private:
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief event called by the timer to allow use to re-draw / animate
    //----------------------------------------------------------------------------------------------------------------------
    void timerEvent(QTimerEvent *);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief process key events
    //----------------------------------------------------------------------------------------------------------------------
    void keyPressEvent(QKeyEvent *);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief resize event
    //----------------------------------------------------------------------------------------------------------------------
    void resizeGL(int _w, int _h);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief a simple draw grid function
    /// @param[in] _size the size of the grid (width and height)
    /// @returns a pointer to the allocated VBO
    void  makeCubes(GLfloat _size);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief advection function
    /// @param[in] args vector of pointers to the data grids to update (velocity component, pressure, temperature, etc.)
    /// @param[in] isCentered vector of bools, true if corresponding grid data is stored at cell centers (pressure-like) or
    /// at faces (velocity-like); must have length equal to args
    /// @param[in] outsideValues the values of corresponding args quantities outside the simulation volume
    /// @param[in] oldIndex true if we use args[i][1] to compute args[i][0], false if it's the other way around
    /// @param[out] args updated with new quantity values
    /// @returns true on success, false on failure (if args and isCentered don't match)
    //----------------------------------------------------------------------------------------------------------------------
    bool advect(std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> args,
                              std::vector<bool> isCentered, std::vector<float> outsideValues);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief advection function (alternate constructor); takes args to contain u,v,w and p
    /// @returns true
    //----------------------------------------------------------------------------------------------------------------------
    bool advect();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief projection function
    //----------------------------------------------------------------------------------------------------------------------
    bool project(std::vector<std::vector<std::vector<SevenPointLagrangianMatrixElement>>> A, std::vector<std::vector<std::vector<float>>> z,
                 std::vector<std::vector<std::vector<float>>> d, std::vector<std::vector<std::vector<float>>> r,
                 std::vector<std::vector<std::vector<float>>> s, std::vector<std::vector<std::vector<float>>> precon,
                 std::vector<std::vector<std::vector<float>>> q);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function to add new animation frame of grid pressure data
    //----------------------------------------------------------------------------------------------------------------------
    void addPressureFrameData();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function to add new animation frame of grid velocity data
    //----------------------------------------------------------------------------------------------------------------------
    void addVelocityFrameData();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function to bake fluid simulation
    //----------------------------------------------------------------------------------------------------------------------
    void bake();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function to put the initial values into the pressure grid
    //----------------------------------------------------------------------------------------------------------------------
    void initializePressure(std::vector<std::vector<std::vector<float>>> pInitial);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function to load the initial values into the velocity grids (u,v,w)
    //----------------------------------------------------------------------------------------------------------------------
    void initializeVelocity(std::vector<std::vector<std::vector<float>>> uInitial, std::vector<std::vector<std::vector<float>>> vInitial, std::vector<std::vector<std::vector<float>>> wInitial);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function to load the initial values into the velocity grids (u,v,w)
    //----------------------------------------------------------------------------------------------------------------------
    float dotProduct(std::vector<std::vector<std::vector<float>>> aMatrix, std::vector<std::vector<std::vector<float>>> bMatrix);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function to load the initial values into the velocity grids (u,v,w)
    //----------------------------------------------------------------------------------------------------------------------
    bool applyA(std::vector<std::vector<std::vector<float>>> aMatrix, std::vector<std::vector<std::vector<float>>> targetMatrix,
                std::vector<std::vector<std::vector<SevenPointLagrangianMatrixElement>>> A);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function to apply preconditioner  to residual "vector" (3d matrix) and set sigma to dotproduct of z and r
    //----------------------------------------------------------------------------------------------------------------------
    bool applyPreconditioner(float& sigma, std::vector<std::vector<std::vector<SevenPointLagrangianMatrixElement>>> A,
                             std::vector<std::vector<std::vector<float>>> z, std::vector<std::vector<std::vector<float>>> d,
                             std::vector<std::vector<std::vector<float>>> r, std::vector<std::vector<std::vector<float>>> s,
                             std::vector<std::vector<std::vector<float>>> precon, std::vector<std::vector<std::vector<float>>> q);


    /// @brief a simple draw grid function
    /// @brief a pointer to our VBO data
    GLuint m_vboPointer=0;
    /// @brief store the size of the vbo data
    GLint m_vboSize=0;
    /// @brief size of the subVBO corresponding to a single cube
    GLint m_cubeSubVBOSize=0;
    /// @brief size of the simulation in grid cells in x, y and z, respectively
    GLint xSimSize=10, ySimSize=10, zSimSize=10;
    /// @brief ID of shader program
    GLint shaderProgramID=0;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief window width from resize event
    //----------------------------------------------------------------------------------------------------------------------
    int m_width;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief window height from resize event
    //----------------------------------------------------------------------------------------------------------------------
    int m_height;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief x offset to center OpenGL canvas in window and keep it square
    //----------------------------------------------------------------------------------------------------------------------
    int m_xOffset;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief y offset to center OpenGL canvas in window and keep it square
    //----------------------------------------------------------------------------------------------------------------------
    int m_yOffset;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief timer to move the stuff
    //----------------------------------------------------------------------------------------------------------------------
    QElapsedTimer timer;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief spin toggle
    //----------------------------------------------------------------------------------------------------------------------
    bool m_spin = false;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the 3d pressure grid values
    //----------------------------------------------------------------------------------------------------------------------
    std::vector<std::vector<std::vector<std::vector<float>>>> p;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the u, v and w velocities (along x, y and z axes, respectively; +y is up)
    //----------------------------------------------------------------------------------------------------------------------
    std::vector<std::vector<std::vector<std::vector<float>>>> u, v, w;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the grid temperatures (for smoke simulation)
    //----------------------------------------------------------------------------------------------------------------------
    std::vector<std::vector<std::vector<std::vector<float>>>> tm;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the timestep
    //----------------------------------------------------------------------------------------------------------------------
    float dt = 1E-2;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the frame duration (1/framerate)
    //----------------------------------------------------------------------------------------------------------------------
    float frameDuration;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the current simulation time
    //----------------------------------------------------------------------------------------------------------------------
    float simTime=0.0;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the total time to simulate
    //----------------------------------------------------------------------------------------------------------------------
    float totalSimTime=0.0;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the default gravitational acceleration
    //----------------------------------------------------------------------------------------------------------------------
    ngl::Vec3 g = ngl::Vec3(0.0,-9.81,0);

    // now project() stuff
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the tolerance for the maximum residual vector element in the projection function
    //----------------------------------------------------------------------------------------------------------------------
    float tol=1E-6;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the maximum number of iterations in the projection function
    //----------------------------------------------------------------------------------------------------------------------
    size_t maxIterations=100;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the tuning constant for Modified Incomplete Cholesky (in project())
    //----------------------------------------------------------------------------------------------------------------------
    const float tau = 0.97;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the density of the fluid
    //----------------------------------------------------------------------------------------------------------------------
    float rho = 1000;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief grid cell size
    //----------------------------------------------------------------------------------------------------------------------
    float dx;

  };

  #endif
