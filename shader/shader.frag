#version 400

layout(location=0) out vec4 out_fragcolor;

in vec3 vposition;
in vec2 vtexcoord;
in float valpha;
in vec3 vnormal;

uniform sampler2D Tex1;
uniform vec3 eye;


//
//	Lights
//
struct DirLight {
	vec3 direction;
	vec3 intensity;
};


//
//	Diffuse color
//
vec3 diffuse(DirLight light, vec3 color, vec3 N){

	vec3 L = -1.f*normalize(light.direction);

	// Material values
	vec3 mamb = color;
	vec3 mdiff = color;

	// Color coeffs
	vec3 ambient = 0.3f * mamb * light.intensity;
	float diff = max( dot(N, L), 0.f );
	vec3 diffuse = diff * mdiff * light.intensity;

	vec3 r = reflect(L, N);
	vec3 v = normalize(eye - vposition);
	vec3 specular = vec3(0.0);
	if (dot(N, L) > 0.0) {
		specular = 0.1 * pow( max(dot(r, v), 0.0), 20)* color * light.intensity;
	}

	return ( ambient + diffuse + specular);
}

void main(){

	DirLight light0, light1;
	light0.direction = vposition - vec3(0, 18, 0);
	light0.intensity = vec3(1, 1, 1);
	light1.direction = vposition - eye;
	light0.intensity = vec3(1, 1, 1);


	vec4 texColor = texture( Tex1, vtexcoord );
	vec3 color = diffuse(light0, vec3(texColor), vnormal);

	color += diffuse(light1, vec3(texColor), vnormal);

	out_fragcolor = vec4(color, valpha);
	// out_fragcolor = vec4(1, 1, 0, 1);
} 
