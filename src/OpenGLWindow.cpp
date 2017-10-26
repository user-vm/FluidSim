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
#include <ngl/Vec3.h>

constexpr float cubeSize=0.2;

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

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_MULTISAMPLE);

  //ngl::Vec3 from(0,1,2);
  //ngl::Vec3 to(0,0,0);
  //ngl::Vec3 up(0,1,0);

  //m_cam.set(from,to,up);
  // set the shape using FOV 45 Aspect Ratio based on Width and Height
  // The final two are near and far clipping planes of 0.5 and 10
  //m_cam.setShape(45,720.0f/576.0f,0.001f,150);

  //get vertex shader

  FILE* vertexShaderFile = fopen("shaders/PhongVertex.glsl","r");

  //Getting File Size
  fseek( vertexShaderFile, 0, SEEK_END );
  long fileSize = ftell( vertexShaderFile );
  rewind( vertexShaderFile );

  char* vertexShaderString = (char*)malloc( sizeof( char) * (fileSize+1) );
  fread( vertexShaderString, sizeof( char ), fileSize, vertexShaderFile );
  vertexShaderString[fileSize] = '\0';
  fclose( vertexShaderFile );

  GLuint vertexShaderID = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource( vertexShaderID, 1, (const GLchar**)&vertexShaderFile, NULL);
  glCompileShader(vertexShaderID);

  //now fragment shader

  FILE* fragmentShaderFile = fopen("shaders/PhongFragment.glsl","r");

  //Getting File Size
  fseek( fragmentShaderFile, 0, SEEK_END );
  fileSize = ftell( fragmentShaderFile );
  rewind( fragmentShaderFile );

  char* fragmentShaderString = (char*)malloc( sizeof( char) * (fileSize+1) );
  fread( fragmentShaderString, sizeof( char ), fileSize, fragmentShaderFile );
  fragmentShaderString[fileSize] = '\0';
  fclose( fragmentShaderFile );

  GLuint fragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource( fragmentShaderID, 1, (const GLchar**)&fragmentShaderFile, NULL);
  glCompileShader(fragmentShaderID);

  //make the shader program and link

  shaderProgramID = glCreateProgram();
  glAttachShader(shaderProgramID, vertexShaderID);
  glAttachShader(shaderProgramID, fragmentShaderID);
  glLinkProgram(shaderProgramID);

  glUseProgram(shaderProgramID);
  glUniform4f(glGetUniformLocation(shaderProgramID,"material.ambient"),0.5,0.5,0.5,1.0);
  glUniform4f(glGetUniformLocation(shaderProgramID,"material.diffuse"),0.5,0.5,0.5,1.0);
  glUniform4f(glGetUniformLocation(shaderProgramID,"material.specular"),0.5,0.5,0.5,1.0);
  glUniform1f(glGetUniformLocation(shaderProgramID,"material.shininess"),0.5);

  glUniform4f(glGetUniformLocation(shaderProgramID,"light.position"),-4.0,-4.0,-4.0,1.0);
  glUniform4f(glGetUniformLocation(shaderProgramID,"light.position"),0.5,0.5,0.5,1.0);
  glUniform4f(glGetUniformLocation(shaderProgramID,"light.position"),0.5,0.5,0.5,1.0);
  glUniform4f(glGetUniformLocation(shaderProgramID,"light.position"),0.5,0.5,0.5,1.0);

  makeCubes(cubeSize);
}



void  OpenGLWindow::makeCubes( GLfloat _size)
{
  // allocate enough space for our verts
  // as we are doing lines it will be 2 verts per line
  // and we need to add 1 to each of them for the <= loop
  // and finally muliply by 12 as we have 12 values per line pair
  qint64 k=-1;
  m_cubeSubVBOSize=3*3*2*6*3;//3 coordinates per vertex, 3 vertices per tri, 2 tris per quad, 6 quads per cube, 3 data types (vertex, normal, colour)
  m_vboSize= m_cubeSubVBOSize*simSize*simSize*simSize*3;

  std::unique_ptr<GLfloat []>vertexData( new GLfloat[m_vboSize]);

  std::vector<ngl::Vec3> verts=
  {
    //12 triangles, two for each face
    //face z=-0.5
    ngl::Vec3(0.5,0.5,-0.5),
    ngl::Vec3(-0.5,0.5,-0.5),
    ngl::Vec3(-0.5,-0.5,-0.5), //+1
    ngl::Vec3(0.5,0.5,-0.5),
    ngl::Vec3(0.5,-0.5,-0.5),
    ngl::Vec3(-0.5,-0.5,-0.5), //-1
    //face z=0.5
    ngl::Vec3(0.5,0.5,0.5),
    ngl::Vec3(-0.5,0.5,0.5),
    ngl::Vec3(-0.5,-0.5,0.5), //
    ngl::Vec3(0.5,0.5,0.5),
    ngl::Vec3(0.5,-0.5,0.5),
    ngl::Vec3(-0.5,-0.5,0.5),
    //face x=-0.5
    ngl::Vec3(-0.5,0.5,0.5),
    ngl::Vec3(-0.5,-0.5,0.5),
    ngl::Vec3(-0.5,-0.5,-0.5),
    ngl::Vec3(-0.5,0.5,0.5),
    ngl::Vec3(-0.5,0.5,-0.5),
    ngl::Vec3(-0.5,-0.5,-0.5),
    //face x=0.5
    ngl::Vec3(0.5,0.5,0.5),
    ngl::Vec3(0.5,-0.5,0.5),
    ngl::Vec3(0.5,-0.5,-0.5),
    ngl::Vec3(0.5,0.5,0.5),
    ngl::Vec3(0.5,0.5,-0.5),
    ngl::Vec3(0.5,-0.5,-0.5),
    //face y=-0.5
    ngl::Vec3(0.5,-0.5,0.5),
    ngl::Vec3(-0.5,-0.5,0.5),
    ngl::Vec3(-0.5,-0.5,-0.5),
    ngl::Vec3(0.5,-0.5,0.5),
    ngl::Vec3(0.5,-0.5,-0.5),
    ngl::Vec3(-0.5,-0.5,-0.5),
    //face y=0.5
    ngl::Vec3(0.5,0.5,0.5),
    ngl::Vec3(-0.5,0.5,0.5),
    ngl::Vec3(-0.5,0.5,-0.5),
    ngl::Vec3(0.5,0.5,0.5),
    ngl::Vec3(0.5,0.5,-0.5),
    ngl::Vec3(-0.5,0.5,-0.5),
  };

  std::vector<ngl::Vec3> normals=
  {
    +1, -1,
  }

  ngl::Vec3 normalHolder;
  bool xCoeff=false, yCoeff=false, zCoeff=false;

  //we will interleave the vertex, normal and colour data

  for(size_t d1=0;d1<simSize;d1++)
    for(size_t d2=0;d2<simSize;d2++)
      for(size_t d3=0;d3<simSize;d3++)
        for(size_t i=0;i<verts.size();i+=3){
          //figure out what the normals will be

          xCoeff=false;
          yCoeff=false;
          zCoeff=false;

          if(verts[i].m_x==verts[i+1].m_x&&verts[i+1].m_x==verts[i+2].m_x)
            xCoeff = true;
          else if(verts[i].m_y==verts[i+1].m_y&&verts[i+1].m_y==verts[i+2].m_y)
            yCoeff = true;
          else if(verts[i].m_z==verts[i+1].m_z&&verts[i+1].m_z==verts[i+2].m_z)
            zCoeff = true;

          //add vertices in three by three, because each vertex in a tri has the same normal vector

          for(t=0;t<2;t++){
            //vertex position
            vertexData[++k]=(verts[i+t].m_x+d1*1.1)*_size;
            vertexData[++k]=(verts[i+t].m_y+d2*1.1)*_size;
            vertexData[++k]=(verts[i+t].m_z+d3*1.1)*_size;

            //vertex normal; factor of two because each quad center is 0.5 away from cube origin
            vertexData[++k]=verts[i+t].m_x * 2 * xCoeff;
            vertexData[++k]=verts[i+t].m_y * 2 * yCoeff;
            vertexData[++k]=verts[i+t].m_z * 2 * zCoeff;

            //vertex color
            vertexData[++k]=d1*1.0/simSize;
            vertexData[++k]=d2*1.0/simSize;
            vertexData[++k]=d3*1.0/simSize;

              }
          //vertex
          vertexData[++k]=(verts[i].m_x+d1*1.1)*_size;
          vertexData[++k]=(verts[i].m_y+d2*1.1)*_size;
          vertexData[++k]=(verts[i].m_z+d3*1.1)*_size;
          //normals; we used the knowledge that they must point outwards
          if(verts[i].m_x==verts[i+1].m_x&&verts[i+1].m_x==verts[i+2].m_x){
            vertexData[++k]=Vec3(verts[i].m_x*2,0,0);
            vertexData[++k]=Vec3(verts[i].m_x*2,0,0);
            vertexData[++k]=Vec3(verts[i].m_x*2,0,0);
            }
          else
            if(verts[i].m_y==verts[i+1].m_y&&verts[i+1].m_y==verts[i+2].m_y){
              vertexData[++k]=Vec3(0,verts[i].m_y*2,0);
              vertexData[++k]=Vec3(0,verts[i].m_y*2,0);
              vertexData[++k]=Vec3(0,verts[i].m_y*2,0);
              }
            else
              if(verts[i].m_z==verts[i+1].m_z&&verts[i+1].m_z==verts[i+2].m_z){
                vertexData[++k]=Vec3(0,0,verts[i].m_z*2);
                vertexData[++k]=Vec3(0,0,verts[i].m_z*2);
                vertexData[++k]=Vec3(0,0,verts[i].m_z*2);
              }
          //colours
          vertexData[++k]=Vec3()
    }

  // now we will create our VBO first we need to ask GL for an Object ID
  glGenBuffers(1, &m_vboPointer);
  // now we bind this ID to an Array buffer
  //glUseProgram(shaderProgramID);
  glBindBuffer(GL_ARRAY_BUFFER, m_vboPointer);
  // finally we stuff our data into the array object
  // First we tell GL it's an array buffer
  // then the number of bytes we are storing (need to tell it's a sizeof(FLOAT)
  // then the pointer to the actual data
  // Then how we are going to draw it (in this case Statically as the data will not change)
  glBufferData(GL_ARRAY_BUFFER, m_vboSize*sizeof(GL_FLOAT) , vertexData.get(), GL_DYNAMIC_DRAW);

}

void OpenGLWindow::paintGL()
{
  //TODO: update this shit to account for added normal and color data
  //it's also possible that the shaders might work with that added

  // set the viewport
  glViewport(0,0,m_width,m_height);
  // clear the colour and depth buffers ready to draw.
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // enable  vertex array drawing
  glEnableClientState(GL_VERTEX_ARRAY);
  // bind our VBO data to be the currently active one
  //glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,0);
  glUseProgram(shaderProgramID);
  glVertexPointer(3,GL_FLOAT,0,0);

  //glDrawArrays( GL_LINES, 0, m_gridSubVBOSize/3);


  // tell GL how this data is formated in this case 3 floats tightly packed starting at the begining
  // of the data (0 = stride, 0 = offset)
  // draw the VBO as a series of GL_LINES starting at 0 in the buffer and _vboSize*GLfloat
  //glDrawArrays( GL_LINES, 0, m_vboSize);
  //glEnableClientState(GL_VERTEX_ARRAY);
  //glBindBuffer(GL_ARRAY_BUFFER, m_vboPointer);
  //glVertexPointer(3,GL_FLOAT,0,0);
  glPushMatrix();
  //glTranslatef(1.0,0,3.0);
  glRotatef(30.0,0.0,1.0,1.0);
  //glScalef(0.5,0.5,0.5);
  glDrawArrays(GL_TRIANGLES, 0, m_vboSize/3);
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
