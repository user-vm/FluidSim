#version 330 core

// Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 vertexPosition_modelspace;

uniform vec3 lineColor;
uniform mat4 MVP;

out vec3 out_color;

void main(void)
{
    gl_Position =  MVP * vec4(vertexPosition_modelspace,1.0);

    out_color = lineColor;
}
