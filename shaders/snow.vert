#version 330

uniform mat4 modelViewProjection;

in vec4 in_position;
in float in_mass;

out vec4 vertex;
out float mass;

void main() {
  gl_PointSize = 2.0;
  gl_Position = modelViewProjection * in_position;

  vertex = in_position;
  mass = in_mass;
}
