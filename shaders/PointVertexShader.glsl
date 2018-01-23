#version 330 core

// Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 vertexPosition_modelspace;
layout(location = 1) in float pointSize;

out vec3 vertexColorOut;

uniform mat4 MVP;
uniform float maxSize;
uniform float maxSizeCutoff;
uniform vec3 pointColor;

void main(){

    //gl_PointSize = 5.0;
    gl_Position =  MVP * vec4(vertexPosition_modelspace,1.0);

    if(pointSize*maxSize <= maxSizeCutoff)
        gl_PointSize = pointSize*maxSize;
    else
        gl_PointSize = maxSize;

    vertexColorOut = pointColor;
}

