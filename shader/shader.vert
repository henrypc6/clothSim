#version 400

layout(location=0) in vec3 in_position;
layout(location=1) in vec2 in_texcoord;
layout(location=2) in vec3 in_color;
layout(location=3) in float in_alpha;
layout(location=4) in vec3 in_normal;

out vec3 vposition;
out vec2 vtexcoord;
out float valpha;
out vec3 vnormal;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main()
{
	gl_PointSize = 1.0;

	vposition = in_position;
	vtexcoord = in_texcoord;
	valpha = in_alpha;
	vnormal = in_normal;


	vec4 pos = projection * model * view * vec4(in_position,1.0);
	gl_Position = pos;
}
