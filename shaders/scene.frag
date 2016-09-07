#version 400
layout( location = 0 ) out vec4 FragColor;

in vec3 Position;
in vec3 Normal;

uniform vec3 LP;

const float shininess = 16.0;

void main () {

	vec3 normal  = normalize( Normal );							
	vec3 lightDir = normalize(LP - Position);
	vec3 viewDir  = normalize(-Position);
		
	float lightIntensity = 0.6f/length(lightDir);
	lightDir = normalize(lightDir);

	vec3 white = vec3(1.0f, 1.0f, 1.0f);

	//Diffuse part-----------
	float diff = max(dot(lightDir, normal), 0.0);
	vec3 diffuse = diff * white * lightIntensity;

	//specular part-------------
	vec3 H = normalize(lightDir + viewDir);
	float NdH = max(dot(H, normal), 0.0);
	float spec = pow(NdH, shininess);
	vec3 specular = spec * white;

	// Ambient-------------
	vec3 ambient = 0.3*lightIntensity * white;// * texcolor * lightIntensity;
	
	vec3 resultLight = ambient + diffuse;// + specular;
	FragColor = vec4(resultLight, 1.0);
}

