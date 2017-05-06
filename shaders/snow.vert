#version 330

uniform mat4 modelViewProjection;

in vec4 in_position;
in vec4 in_normal;

out vec4 vertex;
out vec4 normal;

void main() {
  gl_PointSize = 2.0;
  gl_Position = modelViewProjection * in_position;

  vertex = in_position;
  normal = in_normal;
}
