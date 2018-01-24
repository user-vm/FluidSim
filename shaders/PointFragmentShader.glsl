#version 330 core

in vec3 vertexColorOut;

// Output data
out vec4 gl_FragColor;

void main()
{

    gl_FragColor = vec4(vertexColorOut,1.0);
}
