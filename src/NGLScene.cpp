#include <QMouseEvent>
#include <QGuiApplication>

#include "NGLScene.h"
#include <ngl/Camera.h>
#include <ngl/Light.h>
#include <ngl/Transformation.h>
#include <ngl/Material.h>
#include <ngl/NGLInit.h>
#include <ngl/SimpleVAO.h>
#include <ngl/VAOPrimitives.h>
#include <ngl/Random.h>
#include <ngl/ShaderLib.h>



NGLScene::NGLScene()
{
  setTitle("Qt5 Simple NGL Demo");
  m_data.resize(simSize);
  std::cout<<"OIHJVRIBWSHJOVERHUION";
}


NGLScene::~NGLScene()
{
  std::cout<<"Shutting down NGL, removing VAO's and Shaders\n";
  m_vao->removeVAO();
}

void NGLScene::resizeGL( int _w, int _h )
{
  m_cam.setShape( 45.0f, static_cast<float>( _w ) / _h, 0.05f, 350.0f );
  m_win.width  = static_cast<int>( _w * devicePixelRatio() );
  m_win.height = static_cast<int>( _h * devicePixelRatio() );
}

void NGLScene::initializeGL()
{
  // we need to initialise the NGL lib which will load all of the OpenGL functions, this must
  // be done once we have a valid GL context but before we call any GL commands. If we dont do
  // this everything will crash
  ngl::NGLInit::instance();

  glClearColor(0.4f, 0.4f, 0.4f, 1.0f);			   // Grey Background
  // enable depth testing for drawing

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_MULTISAMPLE);
  // Now we will create a basic Camera from the graphics library
  // This is a static camera so it only needs to be set once
  // First create Values for the camera position
  ngl::Vec3 from(0,1,22);
  ngl::Vec3 to(0,0,0);
  ngl::Vec3 up(0,1,0);

  m_cam.set(from,to,up);
  // set the shape using FOV 45 Aspect Ratio based on Width and Height
  // The final two are near and far clipping planes of 0.5 and 10
  m_cam.setShape(45,720.0f/576.0f,0.001f,150);

  // now to load the shader and set the values
  // grab an instance of shader manager
  ngl::ShaderLib *shader=ngl::ShaderLib::instance();
  // we are creating a shader called Phong to save typos
  // in the code create some constexpr
  constexpr auto shaderProgram = "Phong";
  constexpr auto vertexShader  = "PhongVertex";
  constexpr auto fragShader    = "PhongFragment";
  // create the shader program
  shader->createShaderProgram( shaderProgram );
  // now we are going to create empty shaders for Frag and Vert
  shader->attachShader( vertexShader, ngl::ShaderType::VERTEX );
  shader->attachShader( fragShader, ngl::ShaderType::FRAGMENT );
  // attach the source
  shader->loadShaderSource( vertexShader, "shaders/PhongVertex.glsl" );
  shader->loadShaderSource( fragShader, "shaders/PhongFragment.glsl" );
  // compile the shaders
  shader->compileShader( vertexShader );
  shader->compileShader( fragShader );
  // add them to the program
  shader->attachShaderToProgram( shaderProgram, vertexShader );
  shader->attachShaderToProgram( shaderProgram, fragShader );


  // now we have associated that data we can link the shader
  shader->linkProgramObject( shaderProgram );
  // and make it active ready to load values
  ( *shader )[ shaderProgram ]->use();
  // the shader will use the currently active material and light0 so set them
  ngl::Material m( ngl::STDMAT::GOLD );
  // load our material values to the shader into the structure material (see Vertex shader)
  m.loadToShader( "material" );
  shader->setUniform("viewerPos",m_cam.getEye().toVec3());
  // now create our light this is done after the camera so we can pass the
  // transpose of the projection matrix to the light to do correct eye space
  // transformations
  ngl::Mat4 iv=m_cam.getViewMatrix();
  iv.transpose();
  iv=iv.inverse();
  ngl::Light l(ngl::Vec3(0,1,0),ngl::Colour(1,1,1,1),ngl::Colour(1,1,1,1),ngl::LightModes::POINTLIGHT);
  l.setTransform(iv);
  // load these values to the shader as well
  l.loadToShader("light");
  /*shader->use("nglColourShader");
  shader->setUniform("Colour",1.0f,1.0f,1.0f,1.0f);
  glViewport(0,0,width(),height());
  startTimer(250);
  elapsedTimer.start();
  glPointSize(10);
  m_text.reset(new  ngl::Text(QFont("Arial",18)));
  m_text->setScreenSize(width(),height());*/
  // create the VAO but don't populate
  m_vao.reset( ngl::VAOFactory::createVAO(ngl::simpleVAO,GL_TRIANGLES));

  m_vao->bind();
  buildVAO();
  /*
  for(auto h : m_data){
      std::cout<<h.m_x<<" "<<h.m_y<<" "<<h.m_z<<"\n";
    }*/
  //std::cout<<"\n";
  m_vao->setData( ngl::SimpleVAO::VertexData(m_data.size()*sizeof(ngl::Vec3),m_data[0].m_x));

  // We must do this each time as we change the data.
  m_vao->setVertexAttributePointer(0,3,GL_FLOAT,0,0);
  m_vao->setNumIndices(m_data.size());

  std::vector <ngl::Vec3> normals;
  makeCubeNormals(normals);

  m_vao->setData(ngl::SimpleVAO::VertexData(normals.size()*sizeof(ngl::Vec3),normals[0].m_x));
  m_vao->setVertexAttributePointer(2,3,GL_FLOAT,0,0);
  m_vao->setNumIndices(m_data.size());

}

/*
    std::vector<ngl::Vec3> verts=
    {
        //12 triangles, two for each face
        //face z=-0.5
        ngl::Vec3(0.5,0.5,-0.5),
        ngl::Vec3(-0.5,0.5,-0.5),
        ngl::Vec3(-0.5,-0.5,-0.5),
        ngl::Vec3(0.5,0.5,-0.5),
        ngl::Vec3(0.5,-0.5,-0.5),
        ngl::Vec3(-0.5,-0.5,-0.5),
        //face z=0.5
        ngl::Vec3(0.5,0.5,0.5),
        ngl::Vec3(-0.5,0.5,0.5),
        ngl::Vec3(-0.5,-0.5,0.5),
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
    };*/

std::vector<ngl::Vec3> NGLScene::makeCube(){//int x,int y,int z){

  std::vector<ngl::Vec3> verts=
  {
    //12 triangles, two for each face
    //face z=-0.5
    ngl::Vec3(0.5,0.5,-0.5),
    ngl::Vec3(-0.5,0.5,-0.5),
    ngl::Vec3(-0.5,-0.5,-0.5),
    ngl::Vec3(0.5,0.5,-0.5),
    ngl::Vec3(0.5,-0.5,-0.5),
    ngl::Vec3(-0.5,-0.5,-0.5),
    //face z=0.5
    ngl::Vec3(0.5,0.5,0.5),
    ngl::Vec3(-0.5,0.5,0.5),
    ngl::Vec3(-0.5,-0.5,0.5),
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

  return verts;
}

void NGLScene::makeCubeNormals(std::vector<ngl::Vec3> &normals){

  qint64 i,j;
  ngl::Vec3 n;
  for(i=0;i<m_data.size();i++){
      if(i%2)
        ngl::Vec3 n=ngl::calcNormal(m_data[i+2],m_data[i+1],m_data[i]);
      else
        ngl::Vec3 n=ngl::calcNormal(m_data[i],m_data[i+1],m_data[i+2]);
      for(j=0;j<3;j++)
        normals.push_back(n);
    }
}

void NGLScene::buildVAO(){

  size_t x,y,z;
  std::vector<ngl::Vec3> cubeVertices;
  for(x=0;x<simSize;x++){
      for(y=0;y<simSize;y++){
          for(z=0;z<simSize;z++){
              cubeVertices = makeCube();
              for(auto it : cubeVertices){
                  it+=ngl::Vec3(x,y,z);
                  it*=simScale;
                  m_data.push_back(it);}
            }
        }
    }

}

void NGLScene::paintGL()
{
  // clear the screen and depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0,0,m_win.width,m_win.height);
  // Rotation based on the mouse position for our global
  // transform
  ngl::Mat4 rotX;
  ngl::Mat4 rotY;
  // create the rotation matrices
  rotX.rotateX(m_win.spinXFace);
  rotY.rotateY(m_win.spinYFace);
  // multiply the rotations
  m_mouseGlobalTX=rotY*rotX;
  // add the translations
  m_mouseGlobalTX.m_m[3][0] = m_modelPos.m_x;
  m_mouseGlobalTX.m_m[3][1] = m_modelPos.m_y;
  m_mouseGlobalTX.m_m[3][2] = m_modelPos.m_z;
  ngl::ShaderLib *shader=ngl::ShaderLib::instance();
  (*shader)["Phong"]->use();

  ngl::Mat4 MVP;
  MVP=m_cam.getVPMatrix()*m_mouseGlobalTX;

  shader->setUniform("MVP",MVP);

  m_vao->bind();
  m_vao->draw();

  m_vao->unbind();

}


void NGLScene::timerEvent(QTimerEvent *_event)
{
  NGL_UNUSED(_event);
  qint64 i = elapsedTimer.elapsed()/1000;

  update();
}

//----------------------------------------------------------------------------------------------------------------------
void NGLScene::mouseMoveEvent( QMouseEvent* _event )
{
  // note the method buttons() is the button state when event was called
  // that is different from button() which is used to check which button was
  // pressed when the mousePress/Release event is generated
  if ( m_win.rotate && _event->buttons() == Qt::LeftButton )
    {
      int diffx = _event->x() - m_win.origX;
      int diffy = _event->y() - m_win.origY;
      m_win.spinXFace += static_cast<int>( 0.5f * diffy );
      m_win.spinYFace += static_cast<int>( 0.5f * diffx );
      m_win.origX = _event->x();
      m_win.origY = _event->y();
      update();
    }
  // right mouse translate code
  else if ( m_win.translate && _event->buttons() == Qt::RightButton )
    {
      int diffX      = static_cast<int>( _event->x() - m_win.origXPos );
      int diffY      = static_cast<int>( _event->y() - m_win.origYPos );
      m_win.origXPos = _event->x();
      m_win.origYPos = _event->y();
      m_modelPos.m_x += INCREMENT * diffX;
      m_modelPos.m_y -= INCREMENT * diffY;
      update();
    }
}


//----------------------------------------------------------------------------------------------------------------------
void NGLScene::mousePressEvent( QMouseEvent* _event )
{
  // that method is called when the mouse button is pressed in this case we
  // store the value where the maouse was clicked (x,y) and set the Rotate flag to true
  if ( _event->button() == Qt::LeftButton )
    {
      m_win.origX  = _event->x();
      m_win.origY  = _event->y();
      m_win.rotate = true;
    }
  // right mouse translate mode
  else if ( _event->button() == Qt::RightButton )
    {
      m_win.origXPos  = _event->x();
      m_win.origYPos  = _event->y();
      m_win.translate = true;
    }
}

//----------------------------------------------------------------------------------------------------------------------
void NGLScene::mouseReleaseEvent( QMouseEvent* _event )
{
  // that event is called when the mouse button is released
  // we then set Rotate to false
  if ( _event->button() == Qt::LeftButton )
    {
      m_win.rotate = false;
    }
  // right mouse translate mode
  if ( _event->button() == Qt::RightButton )
    {
      m_win.translate = false;
    }
}

//----------------------------------------------------------------------------------------------------------------------
void NGLScene::wheelEvent( QWheelEvent* _event )
{

  // check the diff of the wheel position (0 means no change)
  if ( _event->delta() > 0 )
    {
      m_modelPos.m_z += ZOOM;
    }
  else if ( _event->delta() < 0 )
    {
      m_modelPos.m_z -= ZOOM;
    }
  update();
}
//----------------------------------------------------------------------------------------------------------------------

void NGLScene::keyPressEvent(QKeyEvent *_event)
{
  // this method is called every time the main window recives a key event.
  // we then switch on the key value and set the camera in the GLWindow
  switch (_event->key())
    {
    // escape key to quite
    case Qt::Key_Escape : QGuiApplication::exit(EXIT_SUCCESS); break;
      // turn on wirframe rendering
    case Qt::Key_W : glPolygonMode(GL_FRONT_AND_BACK,GL_LINE); break;
      // turn off wire frame
    case Qt::Key_S : glPolygonMode(GL_FRONT_AND_BACK,GL_FILL); break;
      // show full screen
    case Qt::Key_F : showFullScreen(); break;
      // show windowed
    case Qt::Key_N : showNormal(); break;
    default : break;
    }
  // finally update the GLWindow and re-draw
  //if (isExposed())
  update();
}
