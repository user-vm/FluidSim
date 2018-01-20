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
#include <QElapsedTimer>
#include <cmath>
#include <fstream>
#include <sstream>

constexpr float cubeSize=0.05;

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

  glClearColor(0.5f, 0.5f, 0.5f, 1.0f);			   // Grey Background

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_MULTISAMPLE);

  GLint m_viewport[4];

  glGetIntegerv(GL_VIEWPORT,m_viewport);

  projectionMatrix = glm::perspective(glm::radians(initialFoV), 4.0f / 3.0f, 0.1f, 100.0f);

  std::cout<<m_viewport[0]<<" "<<m_viewport[1]<<" "<<m_viewport[2]<<" "<<m_viewport[3]<<"\n";

  pointSize = GLfloat(m_viewport[3])/2.0*cubeSize*10;

  bake(cubeSize);
  makePoints(cubeSize);
}

void  OpenGLWindow::makeCubes( GLfloat _size)
{
  // allocate enough space for our verts
  // as we are doing lines it will be 2 verts per line
  // and we need to add 1 to each of them for the <= loop
  // and finally muliply by 12 as we have 12 values per line pair
  m_cubeSubVBOSize=3*2*6*(6+totalFrames*4);// 6 + totalFrames coordinates per vertex for now, 3 vertices/tri, 2 tris/quad, 6 quads/cube //3 coordinates per vertex, 3 vertices per tri, 2 tris per quad, 6 quads per cube, 3 data types (vertex, normal, colour)
  m_vboSize= m_cubeSubVBOSize*xSimSize*ySimSize*zSimSize;

  // all vertex position data, all normal data, all color data
  m_normalOffset = (3*2*6)*3*xSimSize*ySimSize*zSimSize;
  m_colorOffset = (3*2*6)*6*xSimSize*ySimSize*zSimSize;

  qint64 kVertex=-1, kNormal=m_normalOffset-1, kColor = m_colorOffset-1, kPres = -1;

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
      ngl::Vec3(-0.5,0.5,-0.5)
    };

  ngl::Vec3 normalHolder;
  bool xCoeff=false, yCoeff=false, zCoeff=false;

  //we will interleave the vertex, normal and colour data

  for(size_t d1=0;d1<xSimSize;d1++)
    for(size_t d2=0;d2<ySimSize;d2++)
      for(size_t d3=0;d3<zSimSize;d3++)
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

            for(int t=0;t<3;t++){

                //vertex color
                for(int gc=0;gc<totalFrames;gc++){
                    vertexData[m_vboSize-4-(++kColor)]=mainColorData[++kPres];//0.5 + d1*0.5/xSimSize;// + verts[i+t].m_x;
                    vertexData[m_vboSize-2-(++kColor)]=mainColorData[++kPres];//0.5 + d2*0.5/ySimSize;// + verts[i+t].m_y;
                    vertexData[m_vboSize-(++kColor)]=mainColorData[++kPres];//0.5 + d3*0.5/zSimSize;// + verts[i+t].m_z;
                    vertexData[m_vboSize+2-(++kColor)]=mainColorData[++kPres];}//0.5;}

                //std::cout<<vertexData[k]<<"\n";

                //vertex position
                vertexData[m_normalOffset-3-(++kVertex)]=(verts[i+t].m_x+d1*1 - xSimSize/2.0)*_size;
                vertexData[m_normalOffset-1-(++kVertex)]=(verts[i+t].m_y+d2*1 - ySimSize/2.0)*_size;
                vertexData[m_normalOffset+1-(++kVertex)]=(verts[i+t].m_z+d3*1 - zSimSize/2.0)*_size;

                //vertex normal; factor of two because each quad center is 0.5 away from cube origin
                vertexData[m_colorOffset-3-(++kNormal)]=verts[i+t].m_x * 2 * xCoeff;
                vertexData[m_colorOffset-1-(++kNormal)]=verts[i+t].m_y * 2 * yCoeff;
                vertexData[m_colorOffset+1-(++kNormal)]=verts[i+t].m_z * 2 * zCoeff;

              }

          }

  for(size_t i=0;i<verts.size();i++){
      std::cout<<verts[i].m_x<<" "<<verts[i].m_y<<" "<<verts[i].m_z<<"\n";
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
  //glBufferData(GL_ARRAY_BUFFER, kColor*sizeof(GL_FLOAT) , vertexData.get(), GL_DYNAMIC_DRAW);
  glBufferData(GL_ARRAY_BUFFER, m_vboSize*sizeof(GL_FLOAT) , vertexData.get(), GL_DYNAMIC_DRAW);

  //m_vboSize = kColor;

}

void OpenGLWindow::makePoints(GLfloat _size)
  {
    // allocate enough space for our verts
    // as we are doing lines it will be 2 verts per line
    // and we need to add 1 to each of them for the <= loop
    // and finally muliply by 12 as we have 12 values per line pair
    m_pointSubVBOSize=3+totalFrames*4;//the 7 is from 4 color values// 6 + totalFrames coordinates per vertex for now, 3 vertices/tri, 2 tris/quad, 6 quads/cube //3 coordinates per vertex, 3 vertices per tri, 2 tris per quad, 6 quads per cube, 3 data types (vertex, normal, colour)
    m_fluidVboSize= m_pointSubVBOSize*xSimSize*ySimSize*zSimSize;
    m_solidVboSize = solidFacesData.size();

    m_vboSize = m_fluidVboSize + m_solidVboSize;

    // all vertex position data, all normal data, all color data
    m_colorOffset = 3*xSimSize*ySimSize*zSimSize;

    qint64 kVertex=-1, kColor = m_colorOffset-1, kPres = -1;

    std::unique_ptr<GLfloat []>vertexData( new GLfloat[m_vboSize]);

    for(size_t d1=0;d1<xSimSize;d1++){
      for(size_t d2=0;d2<ySimSize;d2++){
        for(size_t d3=0;d3<zSimSize;d3++){

            //add vertices in three by three, because each vertex in a tri has the same normal vector
            //if(gc==0){
            vertexData[++kVertex]=(d1*1 - xSimSize/2.0)*_size;
            vertexData[++kVertex]=(d2*1 - ySimSize/2.0)*_size;
            vertexData[++kVertex]=(d3*1 - zSimSize/2.0)*_size;//}

            //vertex color
            //std::cout<<"("<<mainColorData[kPres-4]<<" "<<mainColorData[kPres-3]<<" "<<mainColorData[kPres-2]<<" "<<mainColorData[kPres-1]<<") ";
            for(int gc=0;gc<totalFrames;gc++){
                vertexData[++kColor]=mainColorData[++kPres];//0.5 + d1*0.5/xSimSize;// + verts[i+t].m_x;
                vertexData[++kColor]=mainColorData[++kPres];//0.5 + d2*0.5/ySimSize;// + verts[i+t].m_y;
                vertexData[++kColor]=mainColorData[++kPres];//0.5 + d3*0.5/zSimSize;// + verts[i+t].m_z;
                vertexData[++kColor]=mainColorData[++kPres];}//0.5;}

        //std::cout<<"\n";
        }
      //std::cout<<"---new row---\n";
      }
    //std::cout<<"\n+++NEW FRAME+++\n\n";
    }

    std::cout<<"kColor = "<<kColor<<"; kVertex = "<<kVertex<<"\n";

    qint64 kSolid=m_fluidVboSize-1;
    /*
    for(size_t d1=0;d1<m_solidVboSize;d1++)
      //vertexData[++kSolid] = solidFacesData[d1];
      std::cout<<
    */
    // now we will create our VBO first we need to ask GL for an Object ID

    /*for(size_t d1=0;d1<kVertex;d1++)
      std::cout<<vertexData[d1]<<" ";*/

    // Dark blue background
    glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

    // Create and compile our GLSL program from the shaders
    programID = LoadShaders( "shaders/SimpleVertexShader.glsl", "shaders/SimpleFragmentShader.glsl" );

    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    // Get a handle for our "MVP" uniform
    matrixID = glGetUniformLocation(programID, "MVP");

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
    else
      glBufferData(GL_ARRAY_BUFFER, solidFacesData.size()*sizeof(GLfloat), solidFacesData.data(), GL_DYNAMIC_DRAW);
}

// from opengl-tutorial.org
GLuint OpenGLWindow::LoadShaders(const char * vertex_file_path,const char * fragment_file_path){

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

//from opengl-tutorial.org

int a_value=0;

void OpenGLWindow::paintGL()
{

  glViewport(m_xOffset,m_yOffset,m_width,m_height);
  // Clear the screen
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Use our shader
  glUseProgram(programID);

  glm::mat4 ViewMatrix = getViewMatrix();
  glm::mat4 ModelMatrix = glm::mat4(1.0);
  glm::mat4 MVP = ProjectionMatrix * ViewMatrix * ModelMatrix;
  MVP = glm::rotate(MVP, glm::radians(timer.elapsed()/1000.0f*10.0f), glm::vec3(1.0f, 0.0f, 0.0f));

  glUniformMatrix4fv(matrixID, 1, GL_FALSE, &MVP[0][0]);

  //GLfloat transMatVals[16];
  //glGetUniformfv(programID,matrixID,transMatVals);
  //trans = glm::rotate(trans, glm::radians(180.0f), glm::vec3(0.0f, 0.0f, 1.0f));

  // 1rst attribute buffer : vertices
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
  glVertexAttribPointer(
     0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
     3,                  // size
     GL_FLOAT,           // type
     GL_FALSE,           // normalized?
     0,                  // stride
     (void*)0            // array buffer offset
  );

  // Draw the triangle !
  //glRotatef(36,0,0,1); //doesn't work, probably need shader-side rotation
  if(useTriangle)
    glDrawArrays(GL_TRIANGLES, 0, 3);
  else
    glDrawArrays(GL_TRIANGLES, 0, solidFacesData.size()/3); // Starting from vertex 0; 3 vertices total -> 1 triangle
  glDisableVertexAttribArray(0);

  a_value++;

  std::cout<<glGetError()<<"\n";

}

void OpenGLWindow::timerEvent(QTimerEvent *)
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

//based on controls.cpp from
void OpenGLWindow::mouseMoveEvent(QMouseEvent *_event)
{
  // timer is called only once, the first time this function is called
  static double lastTime = timer.elapsed()/1000.0;

  // Compute time difference between current and last frame
  double currentTime = timer.elapsed()/1000.0;
  float deltaTime = float(currentTime - lastTime);

  double xpos, ypos;
  xPos = _event->x();
  yPos = _event->y();
  // Get mouse position

  glfwGetCursorPos(window, &xpos, &ypos);

  // Reset mouse position for next frame
  glfwSetCursorPos(window, 1024/2, 768/2);

  // Compute new orientation
  horizontalAngle += mouseSpeed * float(1024/2 - xpos );
  verticalAngle   += mouseSpeed * float( 768/2 - ypos );

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

  ViewMatrix = glm::lookAt(position,           // Camera is here
                           position+direction, // and looks here : at the same position, plus "direction"
                           up);                // Head is up (set to 0,-1,0 to look upside-down)
}

void OpenGLWindow::keyPressEvent(QKeyEvent *_event)
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
    }
}

void OpenGLWindow::resizeGL(int _w, int _h)
{

  m_width  = static_cast<int>( ((_w<_h)?_w:_h) * devicePixelRatio() );
  m_height = static_cast<int>( ((_w<_h)?_w:_h) * devicePixelRatio() );
  m_xOffset = (_w - m_width)/2;
  m_yOffset = (_h - m_height) / 2;
}
/*
void OpenGLWindow::addPressureFrameData(GridsHolder grids){

  pressureFrameData.resize(pressureFrameData.size() + x_Size*y_Size*z_Size);

}

void OpenGLWindow::addVelocityFrameData(){

}
*/
void OpenGLWindow::initializePressure(GridsHolder *gridsHolder){

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

void OpenGLWindow::initializeTemperature(GridsHolder *gridsHolder){

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

void OpenGLWindow::initializeSmokeConcentration(GridsHolder *gridsHolder){

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

void OpenGLWindow::initializeVelocity(GridsHolder *gridsHolder){

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

void OpenGLWindow::initializeSolid(GridsHolder *gridsHolder){

  ngl::Vec3 solidSize = gridsHolder->getSolidDims();

  for(size_t i=0;i<solidSize.m_x;i++)
    for(size_t j=0;j<solidSize.m_z;j++)
      gridsHolder->setSolid(i,0,j,true);
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

void OpenGLWindow::bake(float _size){

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
          //grids.get()->project(tempDt);
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
  mainColorData = mainFrameData->dataToGLfloat(FrameData::CENTER_POINTS);
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

size_t OpenGLWindow::FrameData::xSize(){

  return x_Size;
}

size_t OpenGLWindow::FrameData::ySize(){

  return y_Size;
}

size_t OpenGLWindow::FrameData::zSize(){

  return z_Size;
}

bool OpenGLWindow::FrameData::addFrame(GridsHolder* gridsHolder, std::string gridName){

  return addFrame(gridsHolder, gridName, num_Frames);
}

bool OpenGLWindow::FrameData::addFrame(GridsHolder* gridsHolder, std::string gridName, size_t atFrame){

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

OpenGLWindow::FrameData::FrameData(size_t xSize, size_t ySize, size_t zSize){

  x_Size = xSize;
  y_Size = ySize;
  z_Size = zSize;
  num_Frames = 0;
}

size_t OpenGLWindow::FrameData::frameSize(){

  return x_Size * y_Size * z_Size;
}

std::unique_ptr<GLfloat[]> OpenGLWindow::FrameData::dataToGLfloat(GLfloatTransformationMethod method){

  //return std::unique_ptr<GLfloat[]>(data.data());
  const int numCubeVars = 4 * 3 * 2 * 6; //4 color values/vertex * 3 vertices/tri * 2 tris/quad * 6 quads/cube

  std::unique_ptr<GLfloat[]> floatData(new GLfloat[frameSize()*numFrames()*((method==WHOLE_CUBE)?numCubeVars:4)]);

  float maxVal = (*(std::max_element(data.begin(),data.end())));

  std::cout<<"MAX="<<maxVal<<"\n";

  size_t k = 0, k2 = 0;

  std::cout<<"METHOD="<<method<<"\n";

  if(method == WHOLE_CUBE) //this is the easiest to implement, start with this
    for(size_t i_fr = 0;i_fr<num_Frames;i_fr++)
      for(size_t i_x=0;i_x<x_Size;i_x++)
        for(size_t i_y=0;i_y<y_Size;i_y++)
          for(size_t i_z=0;i_z<z_Size;i_z++,k2++)
            for(size_t i_c=0;i_c<numCubeVars/4;i_c++){

              floatData[k++] = 1;
              floatData[k++] = 1;
              floatData[k++] = 0;
              floatData[k++] = 1;}//(maxVal==0)?0:(data[k2]/maxVal);}

  if(method == CENTER_POINTS)
    for(size_t i_fr = 0;i_fr<num_Frames;i_fr++)
      for(size_t i_x=0;i_x<x_Size;i_x++)
        for(size_t i_y=0;i_y<y_Size;i_y++)
          for(size_t i_z=0;i_z<z_Size;i_z++,k2++){

            floatData[k++] = 1;
            floatData[k++] = 1;
            floatData[k++] = 0;
            floatData[k++] = 1;}//(maxVal==0)?0:(data[k2]/maxVal);}

  return floatData;

}

size_t OpenGLWindow::FrameData::numFrames(){
  return num_Frames;
}
