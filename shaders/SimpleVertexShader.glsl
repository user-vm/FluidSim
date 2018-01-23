#version 330 core

// Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 vertexPosition_modelspace;
layout(location = 1) in int vertexAxisIn;

out vec3 vertexColorOut;

uniform mat4 MVP;
uniform vec3 colorSolidX;
uniform vec3 colorSolidY;
uniform vec3 colorSolidZ;

void main(){

    gl_Position =  MVP * vec4(vertexPosition_modelspace,1);

    //vertexAxisIn < 0.0 ->x axis; vertexAxisIn == 0.0 -> y axis; vertexAxisIn > 0.0 -> z axis
    if(vertexAxisIn==0.0)
        vertexColorOut = colorSolidY;
    else
        if(vertexAxisIn>0.0)
            vertexColorOut = colorSolidZ;
        else
            vertexColorOut = colorSolidX;
}

