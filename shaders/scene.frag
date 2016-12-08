#version 400
layout( location = 0 ) out vec4 FragColor;

in vec3 Position;
in vec3 Normal;

uniform vec4 Color;

const float shininess = 50.0;

void main () {

	vec3 normal   = normalize(Normal);							
	vec3 lightDir = normalize(vec3(0.0, 2.0, 0.0) - Position);
	vec3 viewDir  = normalize(-Position);
		
	float lightIntensity = 0.6/length(lightDir);
	lightDir = normalize(lightDir);

	vec3 white = vec3(1.0, 1.0, 1.0);

	//Diffuse part-----------
	float diff = max(dot(lightDir, normal), 0.0);
	vec3 diffuse = diff * vec3(Color) * lightIntensity;

	//specular part-------------
	//vec3 H = normalize(lightDir + viewDir);
	//float NdH = max(dot(H, normal), 0.0);
	//float spec = pow(NdH, shininess);
	//vec3 specular = spec * white;

	// Ambient-------------
	vec3 ambient = 0.3*lightIntensity * vec3(Color);
	
	vec3 resultLight = ambient + diffuse;// + specular;
	FragColor = vec4(resultLight, Color.a);
}

