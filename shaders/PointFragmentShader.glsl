#version 330 core

in vec3 pointColor;

// Output data
out vec4 gl_FragColor;

void main()
{

    gl_FragColor = vec4(pointColor,1.0);
}
