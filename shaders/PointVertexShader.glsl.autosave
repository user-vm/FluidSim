#version 330 core

// Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 vertexPosition_modelspace;
layout(location = 1) in vec4 vertexColorIn;

out vec4 vertexColorOut;

uniform mat4 MVP;

void main(){

    //gl_PointSize = 5.0;
    gl_Position =  MVP * vec4(vertexPosition_modelspace,1);

    vertexColorOut = vertexColorIn;
}

