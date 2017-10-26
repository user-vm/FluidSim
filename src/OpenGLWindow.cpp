/*
 * Basic GL Window modified from the example here
 * http://qt-project.org/doc/qt-5.0/qtgui/openglwindow.html
 * adapted to use NGL
 */
#include "OpenGLWindow.h"
#include <QKeyEvent>
#include <QApplication>
#include <memory>
#include <iostream>

constexpr float gridSize=1.5;
constexpr int steps=24;

OpenGLWindow::point3D::point3D(GLfloat x, GLfloat y, GLfloat z){
  m_x = x;
  m_y = y;
  m_z = z;
}

OpenGLWindow::point3D::point3D(){
  point3D(0.0,0.0,0.0);
}

OpenGLWindow::OpenGLWindow()
{
  setTitle("Qt5 compat profile OpenGL 3.2");

}

OpenGLWindow::~OpenGLWindow()
{
  // now we have finished clear the device
  std::cout<<"deleting buffer\n";
  glDeleteBuffers(1,&m_vboPointer);
}



void OpenGLWindow::initializeGL()
{
  glewInit();

  glClearColor(0.0f, 0.0f, 0.0f, 1.0f);			   // Grey Background
  makeGrid(gridSize,steps);
}



void  OpenGLWindow::makeGrid( GLfloat _size, size_t _steps )
{
	// allocate enough space for our verts
	// as we are doing lines it will be 2 verts per line
	// and we need to add 1 to each of them for the <= loop
	// and finally muliply by 12 as we have 12 values per line pair
  int num_cubes=1;
  m_gridSubVBOSize=(_steps+1)*12;
  m_cubeSubVBOSize=3*3*2*6*num_cubes;
  m_vboSize= m_gridSubVBOSize + m_cubeSubVBOSize;
  std::unique_ptr<GLfloat []>vertexData( new GLfloat[m_vboSize]);
        // k is the index into our data set
  int k=-1;
        // claculate the step size for each grid value
  float step=_size/static_cast<float>(_steps);
        // pre-calc the offset for speed
        float s2=_size/2.0f;
        // assign v as our value to change each vertex pair
        float v=-s2;
        // loop for our grid values
  for(size_t i=0; i<=_steps; ++i)
        {
                // vertex 1 x,y,z
                vertexData[++k]=-s2; // x
                vertexData[++k]=v; // y
                vertexData[++k]=0.0; // z

		// vertex 2 x,y,z
		vertexData[++k]=s2; // x
		vertexData[++k]=v; // y
		vertexData[++k]=0.0; // z

		// vertex 3 x,y,z
		vertexData[++k]=v;
		vertexData[++k]=s2;
		vertexData[++k]=0.0;

		// vertex 4 x,y,z
		vertexData[++k]=v;
		vertexData[++k]=-s2;
		vertexData[++k]=0.0;
		// now change our step value
		v+=step;
	}

        std::vector<OpenGLWindow::point3D> verts=
          {
            //12 triangles, two for each face
            //face z=-0.5
            OpenGLWindow::point3D(0.5,0.5,-0.5),
            OpenGLWindow::point3D(-0.5,0.5,-0.5),
            OpenGLWindow::point3D(-0.5,-0.5,-0.5),
            OpenGLWindow::point3D(0.5,0.5,-0.5),
            OpenGLWindow::point3D(0.5,-0.5,-0.5),
            OpenGLWindow::point3D(-0.5,-0.5,-0.5),
            //face z=0.5
            OpenGLWindow::point3D(0.5,0.5,0.5),
            OpenGLWindow::point3D(-0.5,0.5,0.5),
            OpenGLWindow::point3D(-0.5,-0.5,0.5),
            OpenGLWindow::point3D(0.5,0.5,0.5),
            OpenGLWindow::point3D(0.5,-0.5,0.5),
            OpenGLWindow::point3D(-0.5,-0.5,0.5),
            //face x=-0.5
            OpenGLWindow::point3D(-0.5,0.5,0.5),
            OpenGLWindow::point3D(-0.5,-0.5,0.5),
            OpenGLWindow::point3D(-0.5,-0.5,-0.5),
            OpenGLWindow::point3D(-0.5,0.5,0.5),
            OpenGLWindow::point3D(-0.5,0.5,-0.5),
            OpenGLWindow::point3D(-0.5,-0.5,-0.5),
            //face x=0.5
            OpenGLWindow::point3D(0.5,0.5,0.5),
            OpenGLWindow::point3D(0.5,-0.5,0.5),
            OpenGLWindow::point3D(0.5,-0.5,-0.5),
            OpenGLWindow::point3D(0.5,0.5,0.5),
            OpenGLWindow::point3D(0.5,0.5,-0.5),
            OpenGLWindow::point3D(0.5,-0.5,-0.5),
            //face y=-0.5
            OpenGLWindow::point3D(0.5,-0.5,0.5),
            OpenGLWindow::point3D(-0.5,-0.5,0.5),
            OpenGLWindow::point3D(-0.5,-0.5,-0.5),
            OpenGLWindow::point3D(0.5,-0.5,0.5),
            OpenGLWindow::point3D(0.5,-0.5,-0.5),
            OpenGLWindow::point3D(-0.5,-0.5,-0.5),
            //face y=0.5
            OpenGLWindow::point3D(0.5,0.5,0.5),
            OpenGLWindow::point3D(-0.5,0.5,0.5),
            OpenGLWindow::point3D(-0.5,0.5,-0.5),
            OpenGLWindow::point3D(0.5,0.5,0.5),
            OpenGLWindow::point3D(0.5,0.5,-0.5),
            OpenGLWindow::point3D(-0.5,0.5,-0.5),
        };
        for(size_t i=0;i<verts.size();i++){
            vertexData[++k]=verts[i].m_x;
            vertexData[++k]=verts[i].m_y;
            vertexData[++k]=verts[i].m_z;
          }

        // now we will create our VBO first we need to ask GL for an Object ID
  glGenBuffers(1, &m_vboPointer);
        // now we bind this ID to an Array buffer
  glBindBuffer(GL_ARRAY_BUFFER, m_vboPointer);
        // finally we stuff our data into the array object
        // First we tell GL it's an array buffer
        // then the number of bytes we are storing (need to tell it's a sizeof(FLOAT)
        // then the pointer to the actual data
        // Then how we are going to draw it (in this case Statically as the data will not change)
  glBufferData(GL_ARRAY_BUFFER, m_vboSize*sizeof(GL_FLOAT) , vertexData.get(), GL_STATIC_DRAW);

}

void OpenGLWindow::paintGL()
{
  // set the viewport
  glViewport(0,0,m_width,m_height);
  // clear the colour and depth buffers ready to draw.
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // enable  vertex array drawing
  glEnableClientState(GL_VERTEX_ARRAY);
  // bind our VBO data to be the currently active one
  //glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,0);
  glVertexPointer(3,GL_FLOAT,0,0);



  glDrawArrays( GL_LINES, 0, m_gridSubVBOSize/3);


  // tell GL how this data is formated in this case 3 floats tightly packed starting at the begining
  // of the data (0 = stride, 0 = offset)
  // draw the VBO as a series of GL_LINES starting at 0 in the buffer and _vboSize*GLfloat
  //glDrawArrays( GL_LINES, 0, m_vboSize);
  //glEnableClientState(GL_VERTEX_ARRAY);
  glBindBuffer(GL_ARRAY_BUFFER, m_vboPointer);
  //glVertexPointer(3,GL_FLOAT,0,0);
  glPushMatrix();
  //glTranslatef(1.0,0,3.0);
  glRotatef(30.0,0.0,1.0,1.0);
  //glScalef(0.5,0.5,0.5);
  glDrawArrays(GL_TRIANGLES, m_gridSubVBOSize/3, m_cubeSubVBOSize/3);
  glPopMatrix();
  // now turn off the VBO client state as we have finished with it
  glDisableClientState(GL_VERTEX_ARRAY);
}

void OpenGLWindow::timerEvent(QTimerEvent *)
{
  update();
}

void OpenGLWindow::keyPressEvent(QKeyEvent *_event)
{
  switch (_event->key())
  {
   case Qt::Key_Escape : QApplication::exit(EXIT_SUCCESS); break;
  }
}

void OpenGLWindow::resizeGL(int _w, int _h)
{

  m_width  = static_cast<int>( _w * devicePixelRatio() );
  m_height = static_cast<int>( _h * devicePixelRatio() );
}
