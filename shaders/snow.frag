#version 330

uniform vec4 in_color;

in vec4 vertex;
in float mass;

out vec4 color;

void main() {
  color = in_color;
  color.w *= (mass - 0.5) * 2;
}
