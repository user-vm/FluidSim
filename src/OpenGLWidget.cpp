#include "OpenGLWidget.h"

/*
 * Basic GL Widget modified from the example here
 * http://qt-project.org/doc/qt-5.0/qtgui/openglwindow.html
 * adapted to use NGL
 */
#include "OpenGLWidget.h"
#include <QKeyEvent>
#include <QApplication>
#include <memory>
#include <iostream>
#include <ngl/Vec3.h>
#include <QElapsedTimer>
#include <cmath>
#include <fstream>
#include <sstream>
#include <glm/gtc/type_ptr.hpp>

OpenGLWidget::OpenGLWidget(QWidget *parent) : QOpenGLWidget(parent)
{
  //setTitle("Qt5 compat profile OpenGL 3.2");

}

OpenGLWidget::~OpenGLWidget()
{
  // now we have finished clear the device
  std::cout<<"deleting buffer\n";
  glDeleteBuffers(1,&vertexbuffer);
}

void OpenGLWidget::initializeGL()
{

  // Direction : Spherical coordinates to Cartesian coordinates conversion
  glm::vec3 direction(
        cos(verticalAngle) * sin(horizontalAngle),
        sin(verticalAngle),
        cos(verticalAngle) * cos(horizontalAngle)
        );

  // Right vector
  glm::vec3 right = glm::vec3(
        sin(horizontalAngle - 3.14f/2.0f),
        0,
        cos(horizontalAngle - 3.14f/2.0f)
        );

  // Up vector
  glm::vec3 up = glm::cross( right, direction );

  position = initialPosition;
  verticalAngle = initialVerticalAngle;
  horizontalAngle = initialHorizontalAngle;
  focusPoint = initialFocusPoint;

  glewInit();

  //glPushAttrib(GL_ALL_ATTRIB_BITS);
  glClearColor(0.5f, 0.5f, 0.5f, 1.0f);			   // Grey Background

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_MULTISAMPLE);

  GLint m_viewport[4];

  glGetIntegerv(GL_VIEWPORT,m_viewport);
  pointSize = GLfloat(m_viewport[3])/2.0*cubeSize*10;

  std::cout<<m_width<<" "<<m_height<<"\n";

  viewMatrix = glm::lookAt(position,           // Camera is here
                           focusPoint,         // and looks here : at the same position, plus "direction"
                           up);                // Head is up (set to 0,-1,0 to look upside-down)

  modelMatrix = glm::mat4(1.0);

  //std::cout<<m_viewport[0]<<" "<<m_viewport[1]<<" "<<m_viewport[2]<<" "<<m_viewport[3]<<"\n";

  std::cout<<"\nviewMatrix=\n";
  for(int i=0;i<4;i++){
    std::cout<<"[";
    for(int j=0;j<4;j++)
      std::cout<<viewMatrix[i][j]<<" ";
    std::cout<<"]\n";}

  bake(cubeSize);
  makePoints();

  std::cout<<"m_width = "<<m_width<<"; m_height = "<<m_height<<"\n";
}

void OpenGLWidget::reset(){

  //glewInit();

  cubeSize /=2;

  //glPopAttrib();
  glDeleteBuffers(1,&vertexbuffer); //this actually doesn't work for some reason
  //glDeleteProgram(programID);
  //glDeleteProgram(pointProgramID);
  //glFlush();

  // Direction : Spherical coordinates to Cartesian coordinates conversion
  glm::vec3 direction(
        cos(verticalAngle) * sin(horizontalAngle),
        sin(verticalAngle),
        cos(verticalAngle) * cos(horizontalAngle)
        );

  // Right vector
  glm::vec3 right = glm::vec3(
        sin(horizontalAngle - 3.14f/2.0f),
        0,
        cos(horizontalAngle - 3.14f/2.0f)
        );

  // Up vector
  glm::vec3 up = glm::cross( right, direction );

  position = initialPosition;
  verticalAngle = initialVerticalAngle;
  horizontalAngle = initialHorizontalAngle;
  focusPoint = initialFocusPoint;

  //glewInit();

  //glClearColor(0.5f, 0.5f, 0.5f, 1.0f);			   // Grey Background

  //glEnable(GL_DEPTH_TEST);
  //glEnable(GL_MULTISAMPLE);

  //GLint m_viewport[4];

  //glGetIntegerv(GL_VIEWPORT,m_viewport);

  std::cout<<m_width<<" "<<m_height<<"\n";

  viewMatrix = glm::lookAt(position,           // Camera is here
                           focusPoint,         // and looks here : at the same position, plus "direction"
                           up);                // Head is up (set to 0,-1,0 to look upside-down)

  modelMatrix = glm::mat4(1.0);

  //std::cout<<m_viewport[0]<<" "<<m_viewport[1]<<" "<<m_viewport[2]<<" "<<m_viewport[3]<<"\n";

  std::cout<<"\nviewMatrix=\n";
  for(int i=0;i<4;i++){
    std::cout<<"[";
    for(int j=0;j<4;j++)
      std::cout<<viewMatrix[i][j]<<" ";
    std::cout<<"]\n";}

  //pointSize = GLfloat(m_viewport[3])/2.0*cubeSize*10;

  bake(cubeSize);
  makePoints();

  std::cout<<"m_width = "<<m_width<<"; m_height = "<<m_height<<"\n";
}

void OpenGLWidget::makePoints()
  {

    // Dark blue background
    glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

    // Create and compile our GLSL program from the shaders
    programID = LoadShaders( "shaders/SimpleVertexShader.glsl", "shaders/SimpleFragmentShader.glsl" );
    pointProgramID = LoadShaders("shaders/PointVertexShader.glsl", "shaders/PointFragmentShader.glsl");//, "shaders/PointGeometryShader.glsl");

    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    // Get a handle for our "MVP" uniform
    matrixID = glGetUniformLocation(programID, "MVP");
    matrixPointID = glGetUniformLocation(pointProgramID, "MVP");
    colorSolidXID = glGetUniformLocation(programID, "colorSolidX");
    colorSolidYID = glGetUniformLocation(programID, "colorSolidY");
    colorSolidZID = glGetUniformLocation(programID, "colorSolidZ");
    maximumSizeID = glGetUniformLocation(pointProgramID, "maxSize");
    maximumSizeCutoffID = glGetUniformLocation(pointProgramID, "maxSizeCutoff");
    pointColorID = glGetUniformLocation(pointProgramID, "pointColor");

    std::cout<<glGetError()<<"\n";

    static const GLfloat g_vertex_buffer_data[] = {
            -1.0f, -1.0f, 0.0f,
             1.0f, -1.0f, 0.0f,
             0.0f,  1.0f, 0.0f,
    };
    /*
    for(size_t d1=0;d1<solidFacesData.size();d1++)
      solidFacesData[d1]+= 0.25 * (std::rand()*2.0/RAND_MAX-0.5); //adds a random float value between -0.25 and 0.25
    */
    // Generate 1 buffer, put the resulting identifier in vertexbuffer
    glGenBuffers(1, &vertexbuffer);
    // The following commands will talk about our 'vertexbuffer' buffer
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    // Give our vertices to OpenGL.
    if(useTriangle)
      glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_DYNAMIC_DRAW);
    else{
      glBufferData(GL_ARRAY_BUFFER, (solidFacesData.size()+mainColorData.size())*sizeof(GLfloat), NULL, GL_STATIC_DRAW);
      glBufferSubData(GL_ARRAY_BUFFER, 0, solidFacesData.size()*sizeof(GLfloat),solidFacesData.data());
      glBufferSubData(GL_ARRAY_BUFFER, solidFacesData.size()*sizeof(GLfloat), mainColorData.size()*sizeof(GLfloat),mainColorData.data());}

    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
}

// from opengl-tutorial.org (distributed under a DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE)
GLuint OpenGLWidget::LoadShaders(const char * vertex_file_path,const char * fragment_file_path){

	// Create the shaders
	GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

	// Read the Vertex Shader code from the file
	std::string VertexShaderCode;
	std::ifstream VertexShaderStream(vertex_file_path, std::ios::in);
	if(VertexShaderStream.is_open()){
		std::stringstream sstr;
		sstr << VertexShaderStream.rdbuf();
		VertexShaderCode = sstr.str();
		VertexShaderStream.close();
	}else{
		printf("Impossible to open %s. Are you in the right directory ? Don't forget to read the FAQ !\n", vertex_file_path);
		getchar();
		return 0;
	}

	// Read the Fragment Shader code from the file
	std::string FragmentShaderCode;
	std::ifstream FragmentShaderStream(fragment_file_path, std::ios::in);
	if(FragmentShaderStream.is_open()){
		std::stringstream sstr;
		sstr << FragmentShaderStream.rdbuf();
		FragmentShaderCode = sstr.str();
		FragmentShaderStream.close();
	}

	GLint Result = GL_FALSE;
	int InfoLogLength;


	// Compile Vertex Shader
	printf("Compiling shader : %s\n", vertex_file_path);
	char const * VertexSourcePointer = VertexShaderCode.c_str();
	glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL);
	glCompileShader(VertexShaderID);

	// Check Vertex Shader
	glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 0 ){
		std::vector<char> VertexShaderErrorMessage(InfoLogLength+1);
		glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
		printf("%s\n", &VertexShaderErrorMessage[0]);
	}



	// Compile Fragment Shader
	printf("Compiling shader : %s\n", fragment_file_path);
	char const * FragmentSourcePointer = FragmentShaderCode.c_str();
	glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL);
	glCompileShader(FragmentShaderID);

	// Check Fragment Shader
	glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 0 ){
		std::vector<char> FragmentShaderErrorMessage(InfoLogLength+1);
		glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
		printf("%s\n", &FragmentShaderErrorMessage[0]);
	}



	// Link the program
	printf("Linking program\n");
	GLuint ProgramID = glCreateProgram();
	glAttachShader(ProgramID, VertexShaderID);
	glAttachShader(ProgramID, FragmentShaderID);
	glLinkProgram(ProgramID);

	// Check the program
	glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
	glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 0 ){
		std::vector<char> ProgramErrorMessage(InfoLogLength+1);
		glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
		printf("%s\n", &ProgramErrorMessage[0]);
	}


	glDetachShader(ProgramID, VertexShaderID);
	glDetachShader(ProgramID, FragmentShaderID);

	glDeleteShader(VertexShaderID);
	glDeleteShader(FragmentShaderID);

	return ProgramID;
}

// modified from opengl-tutorial.org (original distributed under a DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE)
GLuint OpenGLWidget::LoadPointShaders(const char * vertex_file_path,const char * fragment_file_path, const char * geometry_file_path){

	// Create the shaders
	GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
	GLuint GeometryShaderID = glCreateShader(GL_GEOMETRY_SHADER);

	// Read the Vertex Shader code from the file
	std::string VertexShaderCode;
	std::ifstream VertexShaderStream(vertex_file_path, std::ios::in);
	if(VertexShaderStream.is_open()){
		std::stringstream sstr;
		sstr << VertexShaderStream.rdbuf();
		VertexShaderCode = sstr.str();
		VertexShaderStream.close();
	}else{
		printf("Impossible to open %s. Are you in the right directory ? Don't forget to read the FAQ !\n", vertex_file_path);
		getchar();
		return 0;
	}

	// Read the Fragment Shader code from the file
	std::string FragmentShaderCode;
	std::ifstream FragmentShaderStream(fragment_file_path, std::ios::in);
	if(FragmentShaderStream.is_open()){
		std::stringstream sstr;
		sstr << FragmentShaderStream.rdbuf();
		FragmentShaderCode = sstr.str();
		FragmentShaderStream.close();
	}

	// Read the Fragment Shader code from the file
	std::string GeometryShaderCode;
	std::ifstream GeometryShaderStream(fragment_file_path, std::ios::in);
	if(FragmentShaderStream.is_open()){
		std::stringstream sstr;
		sstr << GeometryShaderStream.rdbuf();
		GeometryShaderCode = sstr.str();
		GeometryShaderStream.close();
	}

	GLint Result = GL_FALSE;
	int InfoLogLength;


	// Compile Vertex Shader
	printf("Compiling shader : %s\n", vertex_file_path);
	char const * VertexSourcePointer = VertexShaderCode.c_str();
	glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL);
	glCompileShader(VertexShaderID);

	// Check Vertex Shader
	glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 0 ){
		std::vector<char> VertexShaderErrorMessage(InfoLogLength+1);
		glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
		printf("%s\n", &VertexShaderErrorMessage[0]);
	}



	// Compile Fragment Shader
	printf("Compiling shader : %s\n", fragment_file_path);
	char const * FragmentSourcePointer = FragmentShaderCode.c_str();
	glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL);
	glCompileShader(FragmentShaderID);

	// Check Fragment Shader
	glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 0 ){
		std::vector<char> FragmentShaderErrorMessage(InfoLogLength+1);
		glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
		printf("%s\n", &FragmentShaderErrorMessage[0]);
	}

	// Compile Geometry Shader
	printf("Compiling shader : %s\n", geometry_file_path);
	char const * GeometrySourcePointer = GeometryShaderCode.c_str();
	glShaderSource(GeometryShaderID, 1, &GeometrySourcePointer , NULL);
	glCompileShader(GeometryShaderID);

	// Check Geometry Shader
	glGetShaderiv(GeometryShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(GeometryShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 0 ){
		std::vector<char> GeometryShaderErrorMessage(InfoLogLength+1);
		glGetShaderInfoLog(GeometryShaderID, InfoLogLength, NULL, &GeometryShaderErrorMessage[0]);
		printf("%s\n", &GeometryShaderErrorMessage[0]);
	}

	// Link the program
	printf("Linking program\n");
	GLuint ProgramID = glCreateProgram();
	glAttachShader(ProgramID, VertexShaderID);
	glAttachShader(ProgramID, FragmentShaderID);
	glAttachShader(ProgramID, GeometryShaderID);
	glLinkProgram(ProgramID);

	// Check the program
	glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
	glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 0 ){
		std::vector<char> ProgramErrorMessage(InfoLogLength+1);
		glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
		printf("%s\n", &ProgramErrorMessage[0]);
	}


	glDetachShader(ProgramID, VertexShaderID);
	glDetachShader(ProgramID, FragmentShaderID);
	glDetachShader(ProgramID, GeometryShaderID);

	glDeleteShader(VertexShaderID);
	glDeleteShader(FragmentShaderID);
	glDeleteShader(GeometryShaderID);


	return ProgramID;
}

//from opengl-tutorial.org

int a_value=0;

void OpenGLWidget::paintGL()
{

  glViewport(m_xOffset,m_yOffset,m_width,m_height);
  // Clear the screen
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_ACCUM_BUFFER_BIT);

  // Use our shader
  glUseProgram(programID);

  //sleep(1);

  projectionMatrix = glm::perspective(glm::radians(initialFoV), m_width*1.0f / m_height, 0.01f, 100.0f);
/*
  std::cout<<"\nprojectionMatrix=\n";
  for(int i=0;i<4;i++){
    std::cout<<"[";
    for(int j=0;j<4;j++)
      std::cout<<projectionMatrix[i][j]<<" ";
    std::cout<<"]\n";}

  std::cout<<"\nviewMatrix=\n";
  for(int i=0;i<4;i++){
    std::cout<<"[";
    for(int j=0;j<4;j++)
      std::cout<<viewMatrix[i][j]<<" ";
    std::cout<<"]\n";}

  std::cout<<"\nmodelMatrix=\n";
  for(int i=0;i<4;i++){
    std::cout<<"[";
    for(int j=0;j<4;j++)
      std::cout<<modelMatrix[i][j]<<" ";
    std::cout<<"]\n";}*/

  glm::mat4 MVP = glm::mat4(1.0);
  MVP = projectionMatrix * viewMatrix * modelMatrix;
/*
  std::cout<<"\nMVP=\n";
  for(int i=0;i<4;i++){
    std::cout<<"[";
    for(int j=0;j<4;j++)
      std::cout<<MVP[i][j]<<" ";
    std::cout<<"]\n";}*/

  //MVP = glm::rotate(MVP, glm::radians(timer.elapsed()/1000.0f*10.0f), glm::vec3(1.0f, 0.0f, 0.0f));

  glUniformMatrix4fv(matrixID, 1, GL_FALSE, &MVP[0][0]);

  std::cout<<glGetError()<<"x1a\n";

  glUniform3fv(colorSolidXID,1, glm::value_ptr(colorSolidX));
  glUniform3fv(colorSolidYID,1, glm::value_ptr(colorSolidY));
  glUniform3fv(colorSolidZID,1, glm::value_ptr(colorSolidZ));

  std::cout<<glGetError()<<"x1\n";

  //GLfloat transMatVals[16];
  //glGetUniformfv(programID,matrixID,transMatVals);
  //trans = glm::rotate(trans, glm::radians(180.0f), glm::vec3(0.0f, 0.0f, 1.0f));

  // 1rst attribute buffer : vertices
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);

  std::cout<<glGetError()<<"x\n";

  //glEnableClientState(GL_COLOR_ARRAY);
  glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer); //this needs to be here because otherwise the deleted buffer comes back from the dead and gets used

  //first the vertex positions
  glVertexAttribPointer(
     0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
     3,                  // size
     GL_FLOAT,           // type
     GL_FALSE,           // normalized?
     4*sizeof(GLfloat),  // stride
     (void*)0            // array buffer offset
  );

  //now the color
  glVertexAttribPointer(
        1,
        1,
        GL_INT,
        GL_FALSE,
        4*sizeof(GLfloat),
        (void*)(3*sizeof(GLfloat)));

  if(useTriangle)
    glDrawArrays(GL_TRIANGLES, 0, 3);
  else
    glDrawArrays(GL_TRIANGLES, 0, solidFacesData.size()/4); // Starting from vertex 0; 3 vertices total -> 1 triangle

  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);

  glUseProgram(pointProgramID);

  glUniformMatrix4fv(matrixPointID, 1, GL_FALSE, &MVP[0][0]);
  glUniform1f(maximumSizeID, maximumSize);
  glUniform1f(maximumSizeCutoffID, maximumSizeCutoff);
  glUniform3fv(pointColorID, 1, glm::value_ptr(pointColor));

  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);

  //now the vertex positions of the points
  glVertexAttribPointer(
        0,
        3,
        GL_FLOAT,
        GL_FALSE,
        4*sizeof(GLfloat),
        (void*)(solidFacesData.size()*sizeof(GLfloat)));
        /*
        7*sizeof(GLfloat),
        (void*)(solidFacesData.size()*sizeof(GLfloat)));*/

  //now the color of the points
  glVertexAttribPointer(
        1,
        1,
        GL_FLOAT,
        GL_FALSE,
        4*sizeof(GLfloat),
        (void*)((3+solidFacesData.size())*sizeof(GLfloat)));
        /*
        7*sizeof(GLfloat),
        (void*)((3+solidFacesData.size())*sizeof(GLfloat)));*/

  //glEnable(GL_BLEND);
  //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_PROGRAM_POINT_SIZE);

  uint currentFrame = int(timer.elapsed()/(frameDuration*1000.0)-timerOffset)%totalFrames;

  std::cout<<pointFrameOffset[currentFrame];
  //std::cout<<

  if(useTriangle)
    glDrawArrays(GL_TRIANGLES, 0, 3);
  else
    glDrawArrays(GL_POINTS,pointFrameOffset[currentFrame],pointFrameOffset[currentFrame+1]-pointFrameOffset[currentFrame]);

  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);

  std::cout<<glGetError()<<"\n";

}

void OpenGLWidget::timerEvent(QTimerEvent *)
{
  //std::cout<<timer.elapsed()<<"\n";
  update();
  //return;
  /*
  int frameOffset = 0;//int(timer.elapsed()/(frameDuration*1000)) % pressureFrameData->numFrames();
  glBufferSubData(GL_ARRAY_BUFFER,(GLintptr)(m_colorOffset*sizeof(GLfloat)), //+pressureFrameData->frameSize()*frameOffset*sizeof(float)),
                  pressureFrameData->frameSize()*4*3*2*6*sizeof(GLfloat)*10,&((mainColorData.get())[frameOffset*4*3*2*6]));
  update();*/
}

void OpenGLWidget::mousePressEvent(QMouseEvent *_event){

  std::cout<<"mousePressEvent at "<<timer.elapsed()<<"\n";

  this->setFocus();

  switch(_event->button()){
    case Qt::LeftButton:{
        mousePosOnLeftClick = ngl::Vec2(_event->x(),_event->y());
        posOnLeftClick = position;
        horizontalAngleOnLeftClick = horizontalAngle;
        verticalAngleOnLeftClick = verticalAngle;
        std::cout<<"LB pressed at "<<timer.elapsed()<<"\n";
        std::cout<<"pos= "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";
        break;
      }
    case Qt::RightButton:{
        mousePosOnRightClick = ngl::Vec2(_event->x(),_event->y());
        posOnRightClick = position;
        horizontalAngleOnRightClick = horizontalAngle;
        verticalAngleOnRightClick = verticalAngle;
        break;
      }
    case Qt::MiddleButton:{
        mousePosOnMiddleClick = ngl::Vec2(_event->x(),_event->y());
        posOnMiddleClick = position;
        fpOnMiddleClick = focusPoint;
        horizontalAngleOnMiddleClick = horizontalAngle;
        verticalAngleOnMiddleClick = verticalAngle;
        break;
      }
    default: {}
    }
}

void OpenGLWidget::wheelEvent(QWheelEvent *_event){

  QPoint numDegrees = _event->angleDelta();

  if(!numDegrees.isNull()){
      QPoint numSteps = numDegrees/15;

      // Direction : Spherical coordinates to Cartesian coordinates conversion
      glm::vec3 direction(
            cos(verticalAngle) * sin(horizontalAngle),
            sin(verticalAngle),
            cos(verticalAngle) * cos(horizontalAngle)
            );

      // Right vector
      glm::vec3 right = glm::vec3(
            sin(horizontalAngle - 3.14f/2.0f),
            0,
            cos(horizontalAngle - 3.14f/2.0f)
            );

      // Up vector
      glm::vec3 up = glm::cross( right, direction );

      position += direction * float(numSteps.y()*scrollSpeed);
      //std::cout<<numSteps.x()<<" "<<numSteps.y()<<"\n";

      viewMatrix = glm::lookAt(position,           // Camera is here
                               focusPoint,         // and looks here : at the same position, plus "direction"
                               up);                // Head is up (set to 0,-1,0 to look upside-down)

      std::cout<<"\nviewMatrix=\n";
      for(int i=0;i<4;i++){
        std::cout<<"[";
        for(int j=0;j<4;j++)
          std::cout<<viewMatrix[i][j]<<" ";
        std::cout<<"]\n";}

    }
}

//based on controls.cpp from
void OpenGLWidget::mouseMoveEvent(QMouseEvent *_event)
{
  switch(_event->buttons()){

    //left mouse button rotates camera around focus point
    case Qt::LeftButton :{
        // timer is called only once, the first time this function is called
        //static double lastTime = timer.elapsed()/1000.0;

        // Compute time difference between current and last frame
        //double currentTime = timer.elapsed()/1000.0;
        //float deltaTime = float(currentTime - lastTime);

        // Get mouse position
        double xPos, yPos;
        xPos = _event->x();
        yPos = _event->y();

        //glfwGetCursorPos(window, &xpos, &ypos);

        // Reset mouse position for next frame
        //glfwSetCursorPos(window, 1024/2, 768/2);

        // Compute new orientation
        horizontalAngle = horizontalAngleOnLeftClick + mouseSpeed * float(mousePosOnLeftClick.m_x-xPos);
        verticalAngle   = verticalAngleOnLeftClick - mouseSpeed * float(mousePosOnLeftClick.m_y-yPos);

        // Direction : Spherical coordinates to Cartesian coordinates conversion
        glm::vec3 direction(
              cos(verticalAngle) * sin(horizontalAngle),
              sin(verticalAngle),
              cos(verticalAngle) * cos(horizontalAngle)
              );

        // Right vector
        glm::vec3 right = glm::vec3(
              sin(horizontalAngle - 3.14f/2.0f),
              0,
              cos(horizontalAngle - 3.14f/2.0f)
              );

        // Up vector
        glm::vec3 up = glm::cross( right, direction );
        position = glm::distance(position,focusPoint)*direction + focusPoint;

        viewMatrix = glm::lookAt(position,           // Camera is here
                                 focusPoint,         // and looks here : at the same position, plus "direction"
                                 up);                // Head is up (set to 0,-1,0 to look upside-down)
        break;

      }
    case Qt::MiddleButton :{

        double xPos, yPos;
        xPos = _event->x();
        yPos = _event->y();

        glm::vec3 direction(
              cos(verticalAngle) * sin(horizontalAngle),
              sin(verticalAngle),
              cos(verticalAngle) * cos(horizontalAngle)
              );

        glm::vec3 right = glm::vec3(
              sin(horizontalAngle - 3.14f/2.0f),
              0,
              cos(horizontalAngle - 3.14f/2.0f)
              );

        glm::vec3 up = glm::cross( right, direction );


        position = posOnMiddleClick + mouseSpeed*(right*float(xPos-mousePosOnMiddleClick.m_x) + up*float(yPos-mousePosOnMiddleClick.m_y));
        focusPoint = fpOnMiddleClick + position - posOnMiddleClick; //fpOnMiddleClick - mouseSpeed*(right*float(xPos-mousePosOnMiddleClick.m_x) - up*float(yPos-mousePosOnMiddleClick.m_y));

        std::cout<<"position "<<position[0]<<" "<<position[1]<<" "<<position[2]<<" "<<", focusPoint "<<focusPoint[0]<<" "<<focusPoint[1]<<" "<<focusPoint[2]<<"\n";

        viewMatrix = glm::lookAt(position,           // Camera is here
                                 focusPoint,         // and looks here : at the same position, plus "direction"
                                 up);                // Head is up (set to 0,-1,0 to look upside-down)
        break;
      }
    case Qt::RightButton:{

        double yPos;
        yPos = _event->y();

        // Direction : Spherical coordinates to Cartesian coordinates conversion
        glm::vec3 direction(
              cos(verticalAngle) * sin(horizontalAngle),
              sin(verticalAngle),
              cos(verticalAngle) * cos(horizontalAngle)
              );

        // Right vector
        glm::vec3 right = glm::vec3(
              sin(horizontalAngle - 3.14f/2.0f),
              0,
              cos(horizontalAngle - 3.14f/2.0f)
              );

        // Up vector
        glm::vec3 up = glm::cross( right, direction );

        position = posOnRightClick - mouseSpeed * direction * float(yPos-mousePosOnRightClick.m_y);
        //std::cout<<numSteps.x()<<" "<<numSteps.y()<<"\n";

        viewMatrix = glm::lookAt(position,           // Camera is here
                                 focusPoint,         // and looks here : at the same position, plus "direction"
                                 up);                // Head is up (set to 0,-1,0 to look upside-down)

      }
    default:{}
    }

}

void OpenGLWidget::keyPressEvent(QKeyEvent *_event)
{
  switch (_event->key())
    {
    case Qt::Key_Escape : QApplication::exit(EXIT_SUCCESS); break;
    case ' ' : {
        //pause or play
        if(isPlaying){
            isPlaying = false;
            lastPaused = timer.elapsed();
          }
        else{
            isPlaying = true;
            timerOffset += timer.elapsed() - lastPaused;
          }
        reset();
        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
        break;
      }
    case ',':{
        reset();
        break;
      }
    case Qt::Key_Left : {
        if(isPlaying){
            isPlaying = false;
            lastPaused = timer.elapsed();
          }
        timerOffset += frameDuration*1000;
        break;
      }
    case Qt::Key_Right : {
        if(isPlaying){
            isPlaying = false;
            lastPaused = timer.elapsed();
          }
        timerOffset -= frameDuration*1000;
        break;
      }
    /*case '.':{
        position = glm::vec3(0.0,0.0,0.0);
        glm::vec3 direction(
              cos(verticalAngle) * sin(horizontalAngle),
              sin(verticalAngle),
              cos(verticalAngle) * cos(horizontalAngle)
              );
        glm::vec3 right = glm::vec3(
              sin(horizontalAngle - 3.14f/2.0f),
              0,
              cos(horizontalAngle - 3.14f/2.0f)
              );

        glm::vec3 up = glm::cross( right, direction );
        viewMatrix = glm::lookAt(position,           // Camera is here
                                 position+direction, // and looks here : at the same position, plus "direction"
                                 up);                // Head is up (set to 0,-1,0 to look upside-down)
        break;
      }*/
    default:{}
    }
}

void OpenGLWidget::resizeGL(int _w, int _h)
{

  m_width  = static_cast<int>(_w * devicePixelRatio() );//((_w<_h)?_w:_h) * devicePixelRatio() );
  m_height = static_cast<int>(_h * devicePixelRatio() );//((_w<_h)?_w:_h) * devicePixelRatio() );
  m_xOffset = (_w - m_width)/2;
  m_yOffset = (_h - m_height) / 2;
  //std::cout<<"Resized to"
}
/*
void OpenGLWidget::addPressureFrameData(GridsHolder grids){

  pressureFrameData.resize(pressureFrameData.size() + x_Size*y_Size*z_Size);

}

void OpenGLWidget::addVelocityFrameData(){

}
*/
void OpenGLWidget::resetCamera(){

  position = initialPosition;
  verticalAngle = initialVerticalAngle;
  horizontalAngle = initialHorizontalAngle;
  focusPoint = initialFocusPoint;

  // Direction : Spherical coordinates to Cartesian coordinates conversion
  glm::vec3 direction(
        cos(verticalAngle) * sin(horizontalAngle),
        sin(verticalAngle),
        cos(verticalAngle) * cos(horizontalAngle)
        );

  // Right vector
  glm::vec3 right = glm::vec3(
        sin(horizontalAngle - 3.14f/2.0f),
        0,
        cos(horizontalAngle - 3.14f/2.0f)
        );

  // Up vector
  glm::vec3 up = glm::cross( right, direction );

  viewMatrix = glm::lookAt(position,           // Camera is here
                           focusPoint,         // and looks here : at the same position, plus "direction"
                           up);                // Head is up (set to 0,-1,0 to look upside-down)

  std::cout<<"\nviewMatrix=\n";
  for(int i=0;i<4;i++){
    std::cout<<"[";
    for(int j=0;j<4;j++)
      std::cout<<viewMatrix[i][j]<<" ";
    std::cout<<"]\n";}
}

void OpenGLWidget::initializePressure(GridsHolder *gridsHolder){

  // here we will populate the pressure grid with its initial values

  TwoStepMatrix3D* p = gridsHolder->getTwoStepMatrix3DByName("p");

  size_t x_Size = p->xSize();
  size_t y_Size = p->ySize();
  size_t z_Size = p->zSize();

  p->setOld(x_Size/2, y_Size/2 ,z_Size/2, 10.0);
  p->setOld(x_Size/2+1, y_Size/2 ,z_Size/2, 10.0);
  p->setOld(x_Size/2, y_Size/2+1 ,z_Size/2, 10.0);
  p->setOld(x_Size/2, y_Size/2 ,z_Size/2+1, 10.0);

  /*
  for(size_t i=0;i<x_Size;i++)
    for(size_t j=0;j<y_Size;j++)
      for(size_t k=0;k<z_Size;k++)

        p->setOld(i,j,k,std::fabs(((x_Size*1.0)/2.0-i)/x_Size));*/
}

void OpenGLWidget::initializeTemperature(GridsHolder *gridsHolder){

  // here we will populate the pressure grid with its initial values

  TwoStepMatrix3D* T = gridsHolder->getTwoStepMatrix3DByName("T");

  size_t x_Size = T->xSize();
  size_t y_Size = T->ySize();
  size_t z_Size = T->zSize();

  T->setOld(x_Size/2, y_Size/2 ,z_Size/2, 10.0);
  T->setOld(x_Size/2+1, y_Size/2 ,z_Size/2, 10.0);
  T->setOld(x_Size/2, y_Size/2+1 ,z_Size/2, 10.0);
  T->setOld(x_Size/2+1, y_Size/2+1 ,z_Size/2, 10.0);
  T->setOld(x_Size/2, y_Size/2 ,z_Size/2+1, 10.0);
  T->setOld(x_Size/2+1, y_Size/2 ,z_Size/2+1, 10.0);
  T->setOld(x_Size/2, y_Size/2+1 ,z_Size/2+1, 10.0);
  T->setOld(x_Size/2+1, y_Size/2+1 ,z_Size/2+1, 10.0);

  T->setOld(0,0,z_Size-1,5.0);
  T->setOld(x_Size-1,0,z_Size-1,5.0);
  T->setOld(x_Size-1,0,0,5.0);
  T->setOld(0,0,0,5.0);

}

void OpenGLWidget::initializeSmokeConcentration(GridsHolder *gridsHolder){

  // here we will populate the pressure grid with its initial values

  TwoStepMatrix3D* sc = gridsHolder->getTwoStepMatrix3DByName("sc");

  size_t x_Size = sc->xSize();
  size_t y_Size = sc->ySize();
  size_t z_Size = sc->zSize();

  sc->setOld(x_Size/2, y_Size/2 ,z_Size/2, 10.0);
  sc->setOld(x_Size/2+1, y_Size/2 ,z_Size/2, 10.0);
  sc->setOld(x_Size/2, y_Size/2+1 ,z_Size/2, 10.0);
  sc->setOld(x_Size/2+1, y_Size/2+1 ,z_Size/2, 10.0);
  sc->setOld(x_Size/2, y_Size/2 ,z_Size/2+1, 10.0);
  sc->setOld(x_Size/2+1, y_Size/2 ,z_Size/2+1, 10.0);
  sc->setOld(x_Size/2, y_Size/2+1 ,z_Size/2+1, 10.0);
  sc->setOld(x_Size/2+1, y_Size/2+1 ,z_Size/2+1, 10.0);

  sc->setOld(0,0,z_Size-1,5.0);
  sc->setOld(x_Size-1,0,z_Size-1,5.0);
  sc->setOld(x_Size-1,0,0,5.0);
  sc->setOld(0,0,0,5.0);

}

void OpenGLWidget::initializeVelocity(GridsHolder *gridsHolder){

  //does nothing, all initial velocities will be zero


  // here we will populate the velocity grids with their initial values
  TwoStepMatrix3D* u = gridsHolder->getTwoStepMatrix3DByName("u");
  //TwoStepMatrix3D* v = gridsHolder->getTwoStepMatrix3DByName("v");
  //TwoStepMatrix3D* w = gridsHolder->getTwoStepMatrix3DByName("w");

  size_t x_Size = u->xSize();
  size_t y_Size = u->ySize();
  size_t z_Size = u->zSize();

  for(size_t i=0;i<x_Size;i++)
    for(size_t j=0;j<y_Size;j++)
      for(size_t k=0;k<z_Size;k++)

        u->setOld(i,j,k,std::fabs(((x_Size*1.0)/2.0-i)/x_Size));

}

void OpenGLWidget::initializeSolid(GridsHolder *gridsHolder){

  ngl::Vec3 solidSize = gridsHolder->getSolidDims();

  for(size_t i=0;i<solidSize.m_x;i++)
    for(size_t j=0;j<solidSize.m_z;j++)
      gridsHolder->setSolid(i,0,j,true);

  for(size_t i=2;i<6;i++)
    for(size_t j=3;j<8;j++)
      gridsHolder->setSolid(i,2,j,true);

  gridsHolder->setSolid(1,1,1,true);
  gridsHolder->setSolid(solidSize.m_x-1,1,solidSize.m_z-1,true);
  gridsHolder->setSolid(1,1,solidSize.m_z-1,true);
  gridsHolder->setSolid(solidSize.m_x-1,1,1,true);
}

/*  x_Size = v->xSize();
  y_Size = v->ySize();
  z_Size = v->zSize();

  for(size_t i=0;i<x_Size;i++)
    for(size_t j=0;j<y_Size;j++)
      for(size_t k=0;k<z_Size;k++)

        v->setOld(i,j,k,std::fabs(((y_Size*1.0)/2.0-j)/y_Size));

  x_Size = u->xSize();
  y_Size = u->ySize();
  z_Size = u->zSize();

  for(size_t i=0;i<x_Size;i++)
    for(size_t j=0;j<y_Size;j++)
      for(size_t k=0;k<z_Size;k++)

        w->setOld(i,j,k,std::fabs(((z_Size*1.0)/2.0-k)/z_Size));

}*/

void OpenGLWidget::bake(float _size){

  float dx = cubeSize;

  ngl::Vec3 gDep = g;

  switch(fluidType){
    case WATER:{
        at = 0.0;
        bt = 0.0;
        gDep = g;
        break;
      }
    case SMOKE:{
        gDep = 0.0;
        break;
      }
    case MIX:{}
    }

  //^^check if it works (break statement position)^^

  float simTime = 0.0;

  // need to initialize matrices
/*
  std::vector<std::unique_ptr<GridTuple>> gridsToMake = {std::unique_ptr<GridTuple>(new GridTuple("u",GRID_3D_TWOSTEP,xSimSize+1,ySimSize,zSimSize)),
                                         std::unique_ptr<GridTuple>(new GridTuple("v",GRID_3D_TWOSTEP,xSimSize,ySimSize+1,zSimSize)),
                                         std::unique_ptr<GridTuple>(new GridTuple("w",GRID_3D_TWOSTEP,xSimSize,ySimSize,zSimSize+1)),
                                         std::unique_ptr<GridTuple>(new GridTuple("p",GRID_3D_TWOSTEP,xSimSize,ySimSize,zSimSize)),
                                         std::unique_ptr<GridTuple>(new GridTuple("r",GRID_3D,xSimSize,ySimSize,zSimSize)),
                                         std::unique_ptr<GridTuple>(new GridTuple("z",GRID_3D,xSimSize,ySimSize,zSimSize)),
                                         std::unique_ptr<GridTuple>(new GridTuple("s",GRID_3D,xSimSize,ySimSize,zSimSize)),
                                         std::unique_ptr<GridTuple>(new GridTuple("d",GRID_3D,xSimSize,ySimSize,zSimSize)),
                                         std::unique_ptr<GridTuple>(new GridTuple("precon",GRID_3D,xSimSize,ySimSize,zSimSize)),
                                         std::unique_ptr<GridTuple>(new GridTuple("q",GRID_3D,xSimSize,ySimSize,zSimSize))};
*/

  std::vector<std::unique_ptr<GridTuple>> gridsToMake;

  gridsToMake.push_back(std::unique_ptr<GridTuple>(new GridTuple("u",GRID_3D_TWOSTEP,xSimSize+1,ySimSize,zSimSize)));
  gridsToMake.push_back(std::unique_ptr<GridTuple>(new GridTuple("v",GRID_3D_TWOSTEP,xSimSize,ySimSize+1,zSimSize)));
  gridsToMake.push_back(std::unique_ptr<GridTuple>(new GridTuple("w",GRID_3D_TWOSTEP,xSimSize,ySimSize,zSimSize+1)));
  gridsToMake.push_back(std::unique_ptr<GridTuple>(new GridTuple("p",GRID_3D_TWOSTEP,xSimSize,ySimSize,zSimSize)));
  gridsToMake.push_back(std::unique_ptr<GridTuple>(new GridTuple("T",GRID_3D_TWOSTEP,xSimSize,ySimSize,zSimSize)));
  gridsToMake.push_back(std::unique_ptr<GridTuple>(new GridTuple("sc",GRID_3D_TWOSTEP,xSimSize,ySimSize,zSimSize)));
  gridsToMake.push_back(std::unique_ptr<GridTuple>(new GridTuple("r",GRID_3D,xSimSize,ySimSize,zSimSize)));
  gridsToMake.push_back(std::unique_ptr<GridTuple>(new GridTuple("z",GRID_3D,xSimSize,ySimSize,zSimSize)));
  gridsToMake.push_back(std::unique_ptr<GridTuple>(new GridTuple("s",GRID_3D,xSimSize,ySimSize,zSimSize)));
  gridsToMake.push_back(std::unique_ptr<GridTuple>(new GridTuple("d",GRID_3D,xSimSize,ySimSize,zSimSize)));
  gridsToMake.push_back(std::unique_ptr<GridTuple>(new GridTuple("precon",GRID_3D,xSimSize,ySimSize,zSimSize)));
  gridsToMake.push_back(std::unique_ptr<GridTuple>(new GridTuple("q",GRID_3D,xSimSize,ySimSize,zSimSize)));
  gridsToMake.push_back(std::unique_ptr<GridTuple>(new GridTuple("A",GRID_3D_7PL,xSimSize,ySimSize,zSimSize)));

  ngl::Vec3 solidGridSize = ngl::Vec3(xSimSize,ySimSize,zSimSize);

  std::unique_ptr<GridsHolder> grids = std::unique_ptr<GridsHolder>(new GridsHolder(std::move(gridsToMake), dx, dt, tol, maxIterations, rho, gDep, at, bt, solidGridSize));

  //addPressureFrameData();
  //addVelocityFrameData();
  size_t currentFrame = 0;

  //initializePressure(grids.get());
  initializeVelocity(grids.get());
  initializeTemperature(grids.get());
  initializeSmokeConcentration(grids.get());
  initializeSolid(grids.get());

  //pressureFrameData = std::unique_ptr<FrameData>(new FrameData(xSimSize,ySimSize,zSimSize));
  //tempFrameData = std::unique_ptr<FrameData>(new FrameData(xSimSize,ySimSize,zSimSize));
  mainFrameData = std::unique_ptr<FrameData>(new FrameData(xSimSize,ySimSize,zSimSize));
  //uFrameData = std::unique_ptr<FrameData>(new FrameData(xSimSize+1,ySimSize,zSimSize));
  //vFrameData = std::unique_ptr<FrameData>(new FrameData(xSimSize,ySimSize+1,zSimSize));
  //wFrameData = std::unique_ptr<FrameData>(new FrameData(xSimSize,ySimSize,zSimSize+1));

  //pressureFrameData->addFrame(grids.get(),"p");
  //tempFrameData->addFrame(grids.get(),"T");
  if(fluidType==WATER){
    mainFrameData->addFrame(grids.get(),"p");
  }
  else{
    mainFrameData->addFrame(grids.get(),"sc");
  }
  //uFrameData->addFrame(grids.get(),"u");
  //vFrameData->addFrame(grids.get(),"v");
  //wFrameData->addFrame(grids.get(),"w");

  float tempDt = dt;

  // prepare A matrix; since walls are not implemented yet, all diag values will be 6, and all other values will be -1, except at the (xSimSize-1, ySimSize-1, zSimSize-1) corner

  while(simTime<=totalSimTime-dt){

      TwoStepMatrix3D* p;
      TwoStepMatrix3D* T;
      TwoStepMatrix3D* sc;

      bool makeFrame = false;

      if(size_t((simTime + dt)/frameDuration) > currentFrame){
          makeFrame = true;
          tempDt = (currentFrame + 1) * frameDuration - simTime;
        }

      // only need information for two consectutive time steps
      // ADVECT RUNS WITH [0] AS OLD, BODY REWRITES [1], PROJECT RUNS WITH [1] AS OLD, REPEAT

      // body function is just updating the velocities to account for gravity
      if(fluidType==WATER){
          grids.get()->advect({"u","v","w","p"}, tempDt);
          grids.get()->body(tempDt);
          p = grids.get()->getTwoStepMatrix3DByName("p");
          grids.get()->project(tempDt);
          simTime += tempDt;
          if(makeFrame){
              mainFrameData->addFrame(grids.get(),"p");
              currentFrame++;
              std::cout<<currentFrame<<"/"<<totalSimTime/frameDuration<<"\n";
            }
      }
      else{
          grids.get()->advect({"u","v","w","p","T","sc"}, tempDt);
          grids.get()->bodyBuoy(tempDt);
          T = grids.get()->getTwoStepMatrix3DByName("T");
          sc = grids.get()->getTwoStepMatrix3DByName("sc");
          T->swap();
          sc->swap();
          grids.get()->project(tempDt);
          simTime += tempDt;
          if(makeFrame){
              mainFrameData->addFrame(grids.get(),"sc");
              currentFrame++;
              std::cout<<currentFrame<<"/"<<totalSimTime/frameDuration<<"\n";
            }
      }

      //TwoStepMatrix3D* u = grids.get()->getTwoStepMatrix3DByName("u");
      //TwoStepMatrix3D* v = grids.get()->getTwoStepMatrix3DByName("v");
      //TwoStepMatrix3D* w = grids.get()->getTwoStepMatrix3DByName("w");
      //TwoStepMatrix3D* p = grids.get()->getTwoStepMatrix3DByName("p");
      //TwoStepMatrix3D* T = grids.get()->getTwoStepMatrix3DByName("T");
      //TwoStepMatrix3D* sc = grids.get()->getTwoStepMatrix3DByName("sc");

      tempDt = dt;
    }

  solidFacesData = grids->solidToFaces(-xSimSize/2.0*_size,-ySimSize/2.0*_size,-zSimSize/2.0*_size,xSimSize/2.0*_size,ySimSize/2.0*_size,zSimSize/2.0*_size);

  grids.reset();

  //tempColorData = tempFrameData->dataToGLfloat(FrameData::WHOLE_CUBE);
  mainColorData = mainFrameData->dataToGLfloat(FrameData::CENTER_POINTS,
                                               pointFrameOffset,
                                               ngl::Vec3(-xSimSize/2.0*_size,-ySimSize/2.0*_size,-zSimSize/2.0*_size),
                                               ngl::Vec3(xSimSize/2.0*_size,ySimSize/2.0*_size,zSimSize/2.0*_size));
  //uColorData = uFrameData->dataToGLfloat(FrameData::CUBE_WALL_X);
  //vColorData = vFrameData->dataToGLfloat(FrameData::CUBE_WALL_Y);
  //wColorData = wFrameData->dataToGLfloat(FrameData::CUBE_WALL_Z);

  totalFrames = currentFrame+1;

  std::cout<<"Finished bake.\n";

  startTimer(frameDuration*1000);
  timer.start();

  lastPaused = 0.0;
  isPlaying = true;
}

size_t OpenGLWidget::FrameData::xSize(){

  return x_Size;
}

size_t OpenGLWidget::FrameData::ySize(){

  return y_Size;
}

size_t OpenGLWidget::FrameData::zSize(){

  return z_Size;
}

bool OpenGLWidget::FrameData::addFrame(GridsHolder* gridsHolder, std::string gridName){

  return addFrame(gridsHolder, gridName, num_Frames);
}

bool OpenGLWidget::FrameData::addFrame(GridsHolder* gridsHolder, std::string gridName, size_t atFrame){

  if(atFrame > num_Frames){
      std::cout<<"Index of frame to add is not in scope of FrameData object \""<<_name<<"\"";
    return false;}

  // grid data saves only two timesteps; that's why you need to copy them
  // consider file caching in the future? (whatever that means)

  bool doSizesMatch = true;

  GridType gridType = gridsHolder->getTypeByName(gridName);

  //gridsHolder.getAnyByName();

  if(gridType == GRID_3D){
      Matrix3D* grid = gridsHolder->getMatrix3DByName(gridName);

      if(grid->xSize() != x_Size){
          std::cout<<"X-dimensions of FrameData object \""<<_name<<"\" and grid object \""<<gridName<<"\" do not match.\n";
          doSizesMatch = false;
        }

      if(grid->ySize() != y_Size){
          std::cout<<"Y-dimensions of FrameData object \""<<_name<<"\" and grid object \""<<gridName<<"\" do not match.\n";
          doSizesMatch = false;
        }

      if(grid->zSize() != z_Size){
          std::cout<<"Z-dimensions of FrameData object \""<<_name<<"\" and grid object \""<<gridName<<"\" do not match.\n";
          doSizesMatch = false;
        }

      if(!doSizesMatch)
        return false;

      data.resize(data.size() + x_Size*y_Size*z_Size);

      size_t pos = atFrame * x_Size * y_Size * z_Size;

      for(size_t i=0;i<x_Size;i++)
        for(size_t j=0;j<y_Size;j++)
          for(size_t k=0;k<z_Size;k++,pos++)

            data[pos] = grid->get(i,j,k);

      num_Frames++;
      return true;
    }

  if(gridType == GRID_3D_TWOSTEP){
      TwoStepMatrix3D* grid = gridsHolder->getTwoStepMatrix3DByName(gridName);

      if(grid->xSize() != x_Size){
          std::cout<<"X-dimensions of FrameData object \""<<_name<<"\" and grid object \""<<gridName<<"\" do not match.\n";
          doSizesMatch = false;
        }

      if(grid->ySize() != y_Size){
          std::cout<<"Y-dimensions of FrameData object \""<<_name<<"\" and grid object \""<<gridName<<"\" do not match.\n";
          doSizesMatch = false;
        }

      if(grid->zSize() != z_Size){
          std::cout<<"Z-dimensions of FrameData object \""<<_name<<"\" and grid object \""<<gridName<<"\" do not match.\n";
          doSizesMatch = false;
        }

      if(!doSizesMatch)
        return false;

      data.resize(data.size() + x_Size*y_Size*z_Size);

      size_t pos = atFrame * x_Size * y_Size * z_Size;

      for(size_t i=0;i<x_Size;i++)
        for(size_t j=0;j<y_Size;j++)
          for(size_t k=0;k<z_Size;k++,pos++)

            data[pos] = grid->getOld(i,j,k);

      num_Frames++;

      return true;
    }

  std::cout<<"Grid \""<<gridName<<"\" is neither of type Matrix3D nor TwoStepMatrix3D.\n";
  return false;
}

OpenGLWidget::FrameData::FrameData(size_t xSize, size_t ySize, size_t zSize){

  x_Size = xSize;
  y_Size = ySize;
  z_Size = zSize;
  num_Frames = 0;
}

size_t OpenGLWidget::FrameData::frameSize(){

  return x_Size * y_Size * z_Size;
}

std::vector<GLfloat> OpenGLWidget::FrameData::dataToGLfloat(GLfloatTransformationMethod method, std::vector<GLint> &pointFrameOffset, ngl::Vec3 minCoords, ngl::Vec3 maxCoords){

  //return std::unique_ptr<GLfloat[]>(data.data());
  const int numCubeVars = 4 * 3 * 2 * 6; //4 color values/vertex * 3 vertices/tri * 2 tris/quad * 6 quads/cube

  std::vector<GLfloat> floatData = std::vector<GLfloat>();
  pointFrameOffset.resize(num_Frames+1);

  float maxVal = (*(std::max_element(data.begin(),data.end())));

  std::cout<<"MAX="<<maxVal<<"\n";

  size_t k2 = 0;

  std::cout<<"METHOD="<<method<<"\n";

  if(maxVal==0)
    return floatData;

  pointFrameOffset[0] = 0;

  if(method == WHOLE_CUBE) //this is the easiest to implement, start with this
    for(size_t i_fr = 0;i_fr<num_Frames;i_fr++)
      for(size_t i_x=0;i_x<x_Size;i_x++)
        for(size_t i_y=0;i_y<y_Size;i_y++)
          for(size_t i_z=0;i_z<z_Size;i_z++,k2++)
            for(size_t i_c=0;i_c<numCubeVars/4;i_c++){

                floatData.push_back(1);
                floatData.push_back(1);
                floatData.push_back(0);
                floatData.push_back((maxVal==0)?0:(data[k2]/maxVal));}

  if(method == CENTER_POINTS)
    for(size_t i_fr = 0;i_fr<num_Frames;i_fr++){
        for(size_t i_x=0;i_x<x_Size;i_x++)
          for(size_t i_y=0;i_y<y_Size;i_y++)
            for(size_t i_z=0;i_z<z_Size;i_z++,k2++){

              //if(data[k2]){
                  floatData.push_back((minCoords.m_x*(x_Size-i_x-0.5)+maxCoords.m_x*(i_x+0.5))/x_Size);
                  floatData.push_back((minCoords.m_y*(y_Size-i_y-0.5)+maxCoords.m_y*(i_y+0.5))/y_Size);
                  floatData.push_back((minCoords.m_z*(z_Size-i_z-0.5)+maxCoords.m_z*(i_z+0.5))/z_Size);
                  floatData.push_back((maxVal==0)?0:(data[k2]/maxVal));}

        pointFrameOffset[i_fr+1] = floatData.size()/4;}

  return floatData;

}

size_t OpenGLWidget::FrameData::numFrames(){
  return num_Frames;
}
