#ifndef OPENGLWINDOW_H_
#define OPENGLWINDOW_H_
#include <GL/glew.h>
#include <QOpenGLWindow>
#include <QElapsedTimer>

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
    class point3D{
    public:
      point3D(GLfloat x, GLfloat y, GLfloat z);
      point3D();
      GLfloat m_x, m_y, m_z;
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
    /// @param[in] _step sxstep the spacing between grid points
    /// @param[out] o_dataSize the size of the buffer allocated
    /// @returns a pointer to the allocated VBO
    void  makeCubes(GLfloat _size);
    /// @brief a simple second draw grid function
    /// @param[in] _size the size of the grid (width and height)
    /// @param[in] _step sxstep the spacing between grid points
    /// @param[out] o_dataSize the size of the buffer allocated
    /// @returns a pointer to the allocated VBO
    void  makeGrid1(GLfloat _size, size_t _steps);
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
    std::vector<std::vector<std::vector<float>>> p;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the u, v and w velocities (along x, y and z axes, respectively
    //----------------------------------------------------------------------------------------------------------------------
    std::vector<std::vector<std::vector<float>>> u, v, w;
  };

  #endif
