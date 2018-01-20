#ifndef OPENGLWINDOW_H_
#define OPENGLWINDOW_H_
#include <GL/glew.h>
#include <QOpenGLWindow>
#include <QElapsedTimer>
#include <ngl/Vec3.h>
#include <cstdlib>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "Bake.h"

enum FluidType {SMOKE=0, WATER=1, MIX=2};

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

    GLuint LoadShaders(const char * vertex_file_path,const char * fragment_file_path);

    GLuint programID;

    GLuint matrixID;
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
    /// @brief process key events
    //----------------------------------------------------------------------------------------------------------------------
    void mouseMoveEvent(QKeyEvent *);
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
    /// @brief a simple draw grid function
    /// @param[in] _size the size of the grid (width and height)
    /// @returns a pointer to the allocated VBO
    void  makePoints(GLfloat _size);
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
    void bake(float _size);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function to put the initial values into the pressure grid
    //----------------------------------------------------------------------------------------------------------------------
    void initializePressure(GridsHolder* gridsHolder);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function to put the initial values into the temperature grid
    //----------------------------------------------------------------------------------------------------------------------
    void initializeTemperature(GridsHolder* gridsHolder);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function to put the initial values into the temperature grid
    //----------------------------------------------------------------------------------------------------------------------
    void initializeSmokeConcentration(GridsHolder* gridsHolder);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function to load the initial values into the velocity grids (u,v,w)
    //----------------------------------------------------------------------------------------------------------------------
    void initializeVelocity(GridsHolder* gridsHolder);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function to load the solid cell positions into the GridsHolder object
    //----------------------------------------------------------------------------------------------------------------------
    void initializeSolid(GridsHolder* gridsHolder);

    // REWRITE THIS

    /// @brief a simple draw grid function
    /// @brief a pointer to our VBO data
    GLuint m_vboPointer=0;
    /// @brief store the size of the vbo data
    GLint m_vboSize=0;
    /// @brief size of the subVBO corresponding to a single cube
    GLint m_cubeSubVBOSize=0;
    /// @brief size of the subVBO corresponding to a single point (when doing the GL_POINTS render)
    GLint m_pointSubVBOSize=0;
    /// @brief size of the VBO portion corresponding to the tris of the solid cells
    GLint m_solidVboSize;
    /// @brief size of the VBO portion corresponding to the fluid points (when doing the GL_POINTS render)
    GLint m_fluidVboSize;
    /// @brief size of the simulation in grid cells in x, y and z, respectively
    GLint xSimSize=20, ySimSize=20, zSimSize=20;
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
    float frameDuration = 1.0 / 25.0;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the current simulation time
    //----------------------------------------------------------------------------------------------------------------------
    float simTime=0.0;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the total time to simulate
    //----------------------------------------------------------------------------------------------------------------------
    float totalSimTime=1.0;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the default gravitational acceleration
    //----------------------------------------------------------------------------------------------------------------------
    ngl::Vec3 g = ngl::Vec3(0.0,-9.81,0);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the bouyancy force coefficients
    //----------------------------------------------------------------------------------------------------------------------
    float at = 1, bt = 2;

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
    float rho = 0.0001;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief grid cell size
    //----------------------------------------------------------------------------------------------------------------------
    float dx;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief size of points for display
    //----------------------------------------------------------------------------------------------------------------------
    float pointSize;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief type of fluid (water or smoke)
    //----------------------------------------------------------------------------------------------------------------------
    FluidType fluidType = SMOKE;

    size_t m_normalOffset;

    size_t m_colorOffset;

    bool isPlaying = true;
    float timerOffset = 0.0;
    float lastPaused = 0.0;

    //typedef std::vector<std::vector<std::vector<std::vector<float>>>> Matrix4D;

    std::unique_ptr<GLfloat[]> tempColorData;
    std::unique_ptr<GLfloat[]> scColorData;
    std::unique_ptr<GLfloat[]> uColorData;
    std::unique_ptr<GLfloat[]> vColorData;
    std::unique_ptr<GLfloat[]> wColorData;
    std::unique_ptr<GLfloat[]> mainColorData;

    std::vector<GLfloat> solidFacesData;

    // Initial position : on +Z
    glm::vec3 position = glm::vec3( 0, 0, 5 );
    // Initial horizontal angle : toward -Z
    float horizontalAngle = 3.14f;
    // Initial vertical angle : none
    float verticalAngle = 0.0f;
    // Initial Field of View
    float initialFoV = 45.0f;

    float speed = 3.0f; // 3 units / second
    float mouseSpeed = 0.005f;

    glm::mat4 projectionMatrix;

    // something that can be fed into the shader
    class FrameData{

    public:
      enum GLfloatTransformationMethod {WHOLE_CUBE = 0, CUBE_WALL_X = 1, CUBE_WALL_Y = 2, CUBE_WALL_Z = 3, CENTER_POINTS = 4};

      FrameData(size_t xSize, size_t ySize, size_t zSize);
      //~FrameData();

      size_t xSize();
      size_t ySize();
      size_t zSize();
      size_t numFrames();

      size_t frameSize();

      // frame number start from zero
      bool set(size_t x, size_t y, size_t z,size_t frame,float value);
      float get(size_t x, size_t y, size_t z,size_t frame);

      bool setFrame();
      //getFrame(); -> shared pointer?

      bool addFrame(GridsHolder* gridsHolder, std::string gridName);
      bool addFrame(GridsHolder* gridsHolder, std::string gridName, size_t index);

      std::unique_ptr<GLfloat[]> dataToGLfloat(GLfloatTransformationMethod method);

    private:
      std::vector<float> data;
      std::string _name; // optional; default is "unnamed"
      size_t x_Size;
      size_t y_Size;
      size_t z_Size;
      size_t num_Frames;

    };

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief data holding pressure of each cell at each frame
    //----------------------------------------------------------------------------------------------------------------------
    std::unique_ptr<FrameData> pressureFrameData;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief data holding temperature of each cell at each frame
    //----------------------------------------------------------------------------------------------------------------------
    std::unique_ptr<FrameData> tempFrameData;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief data holding smoke concentration of each cell at each frame
    //----------------------------------------------------------------------------------------------------------------------
    std::unique_ptr<FrameData> scFrameData;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief data holding data to be displayed (pressure for water, concentration for smoke)
    //----------------------------------------------------------------------------------------------------------------------
    std::unique_ptr<FrameData> mainFrameData;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief data holding velocity of each cell wall at each frame
    //----------------------------------------------------------------------------------------------------------------------
    std::unique_ptr<FrameData> uFrameData;
    std::unique_ptr<FrameData> vFrameData;
    std::unique_ptr<FrameData> wFrameData;

    //test iterator, remove later
    size_t iter=0;

    //total number of frames, calculate later
    size_t totalFrames = 0;

    // This will identify our vertex buffer
    GLuint vertexbuffer;

    //temporary stuff you should delete later
    bool useTriangle = false;
};

#endif
