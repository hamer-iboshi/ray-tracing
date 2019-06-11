#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>

#define MAX_RAY_DEPTH 5
#define X 0
#define Y 1
#define Z 2
#define max(a,b) (a > b) ? a : b

typedef struct {
	int direction[3];
} Ray;

typedef struct {
	float center[3]; //position of the sphere
	float radius, radius2; //sphere radius and radius^2
	float transparency, reflection; //surface transparency and reflectivity
	float surface_color[3], emission_color[3]; //surface color and emission (light)
	char *label;
} Object;

void subtractXYZ(float *result, float *p1, float *p2){
	result[0] = p1[0] - p2[0];
	result[1] = p1[1] - p2[1];
	result[2] = p1[2] - p2[2];
}

float dot(float *v1, float *v2){
	return v1[X] * v2[X] + v1[Y] * v2[Y] + v1[Z] * v2[Z];
}
//functions of vec3
void normalize(float *vec){
	float nor2 = dot(vec,vec);
	if (nor2 > 0) {
		float inv_nor = 1 / sqrt(nor2);
		vec[X] = vec[X] * inv_nor;
		vec[Y] = vec[Y] * inv_nor;
		vec[Z] = vec[Z] * inv_nor;
	}
}


//functions of Sphere
bool intersect(Object *obj, float *rayorig, float *raydir, float *t0, float *t1){
	float l[3], tca, d2, thc;
	subtractXYZ(l,obj->center,rayorig); //center - rayorig
	tca = dot(l,raydir); //l.dot(raydir)
	if (tca < 0) return false;
	d2 = dot(l,l)- (tca * tca) ; //l.dot(l)
	if(d2 > obj->radius2) return false;
	thc = sqrt(obj->radius2 - d2);
	*t0 = tca - thc;
	*t1 = tca + thc;
	return true;
}


float mix(float a, float b, float mx){
	return b * mx + a * (1 - mx);
}

/* main trace function, takes a ray as argument (defined by its origin and direction).Test if this ray intersects any of the geometry in the scene.
If the ray intersects an object, compute the intersection point, the normal at the intersection point, and shade this point using this information.
Shading depends on the surface property (is it transparent, reflective, diffuse). The function returns a color for the ray.
If the ray intersects an object that is the color of the object at the intersection point, otherwise it returns the background color. */
float* trace(float ray_orig[3], float ray_dir[3], int depth, Object *obj,int obj_size){
	float tnear = INFINITY;
	int i;
	Object *object = NULL;
	static float color[3];
	// find intersection of this ray with the sphere in the scene
	for (i=0; i < obj_size; i++) {
		float t0 = INFINITY, t1 = INFINITY;
		if (intersect(&obj[i],ray_orig,ray_dir,&t0,&t1)) {
			if (t0 < 0) t0 = t1;
			if (t0 < tnear) {
				tnear = t0;
				object = &obj[i];
			}
		}
	}

	// if there's no intersection return black or background color
	if (!object) {
		for (i=0; i<3; i++)
			color[i] = 2;
		// printf("out\n");
			return color;
	}
	float surfacecolor[3] = { 0.0, 0.0, 0.0 };
	float phit[3], nhit[3];
	for (i=0; i<3; i++)
		phit[i] = ray_orig[i] + ray_dir[i] * tnear;
	subtractXYZ(nhit,phit,object->center);
	normalize(nhit); //normalize normal direction
	// If the normal and the view direction are not opposite to each other reverse the normal direction.
	//That also means we are inside the sphere so set the inside bool to true. Finally reverse the sign of IdotN which we want positive.
	float bias = 1e-4; //add some vias to the poin to be traced
	bool inside = false;
	if (dot(ray_dir,nhit) > 0) {
		for (i=0; i<3; i++){
			nhit[i] = -nhit[i];
		}
		inside = true;
	}
	float ray_orig2[3];
	float refldir[3];

	if ((object->transparency > 0 || object->reflection > 0) && depth < MAX_RAY_DEPTH) {
		float facingratio = -dot(ray_dir,nhit);
		// change the mix value to tweak the effect
		float fresneleffect = mix(pow(1-facingratio,4),1,0.1);
		//compute reflection direction
		float dot_raydir_nhit = dot(ray_dir,nhit);
		for(i=0; i<3; i++){
			refldir[i]=ray_dir[i] - nhit[i]*2*dot_raydir_nhit;
			ray_orig2[i] = phit[i] + nhit[i] * bias;
		}
		normalize(refldir);
		float * reflection = trace(ray_orig2, refldir, depth+1, obj, obj_size);
		float * refraction = NULL;
		// if the sphere is also transparent compute refraction ray
		// if (object->transparency) {
		// 	float ior = 1.1, eta = (inside) ? ior : 1/ior; //inside the surface?
		// 	float cosi = -dot(nhit,ray_dir);
		// 	float k = 1 - eta*eta*(1 - cosi*cosi);
		// 	float ksqrt = sqrt(k);
		// 	for(i=0; i<3; i++){
		// 		refldir[i] = ray_dir[i]*eta + nhit[i] * (eta * cosi - ksqrt);
		// 		ray_orig2[i] = phit[i] + nhit[i] * bias;
		// 	}
		// 	normalize(refldir);
		// 	printf("refdir [%f %f %f]\n",refldir[0],refldir[1],refldir[2]);
		// 	printf("origin [%f %f %f]\n",ray_orig2[0],ray_orig2[1],ray_orig2[2]);
		// 	refraction = trace(ray_orig2, refldir, depth+1, obj, obj_size);
		// }
		float refraction_aux[3] = { 0.0, 0.0, 0.0 };
		if(refraction == NULL) refraction = refraction_aux;
		for(i=0; i<3; i++)
			surfacecolor[i] = (reflection[i]*fresneleffect + refraction[i] * (1 - fresneleffect) * object->transparency) * object->surface_color[i];
	} else {
		//diffuse object
		for (i=0; i < obj_size; i++){
			if (obj[i].emission_color[X] > 0) {

				//light
				float transmission[3] = { 1.0, 1.0, 1.0 };
				float light_direction[3];
				for(int j=0; j<3; j++)
					light_direction[j] = obj[i].center[j] - phit[j];
				normalize(light_direction);
				for(int j=0; j < obj_size; j++){
					if( i != j){
						float t0, t1;
						for(int k=0; k<3; k++)
							ray_orig2[k] = phit[k] + nhit[k] * bias;
						if(intersect(&obj[j],ray_orig2,light_direction, &t0, &t1)){
							for(int k=0; k<3; k++)
								transmission[k] = 0;
							break;
						}
					}
				}
				for(int j=0; j<3; j++){
					surfacecolor[j] += object->surface_color[j] * transmission[j] *( dot(nhit,light_direction)) * obj[i].emission_color[j];
				}
			}
		}

	}
	for(i=0; i<3; i++)
		color[i] = surfacecolor[i] + object->emission_color[i];
	return color;
}

void render(Object *obj, int obj_size){
	int width = 3840, height = 2160;
	// int width = 640, height = 480;
	float ***image = malloc(width*sizeof(float**));
	float inv_width = 1 / (0.0+width), inv_height = 1 / (0.0+height);
	float fov = 30, aspectratio = width / (0.0+height);
	float angle = tan(M_PI * 0.5 * (fov / (0.0 + 180)));
	int i, j, k;
	for (i = 0; i < width; i++) {
		image[i] = (float **) malloc(height*sizeof(float*));
		for (j = 0; j < height; j++) {
		  image[i][j] = (float*) malloc(3*sizeof(float));
		}
	}
	//trace rays
	float xx , yy;
	float raydir[3], *pixel, rayorig[3] = { 0, 0, 0 }; //cam origin
	for (j=0; j<height; j++) {
		for (i=0; i<width; i++) {
			xx = (2 * ((i + 0.5) * inv_width) - 1) * angle * aspectratio;
			yy = (1 - 2 * ((j + 0.5) * inv_height)) * angle;
			raydir[0] = xx; raydir[1] = yy; raydir[2] = -1;
			normalize(raydir);
			pixel = trace(rayorig, raydir, 0, obj, obj_size);
			for(k=0; k<3; k++)
				image[i][j][k]=pixel[k];
		}
	}

	//save result to a PPM image
	FILE *fp = fopen("first.ppm", "wb"); /* b - binary mode */
	(void) fprintf(fp, "P6\n%d %d\n255\n", width, height);
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			static unsigned char color[3];
			color[0] = image[i][j][0]*255;
			color[1] = image[i][j][1]*255;
			color[2] = image[i][j][2]*255;
			(void) fwrite(color, 1, 3, fp);
		}
	}
	(void) fclose(fp);
	for (i = 0; i < width; i++) {
		for (j = 0; j < height; j++) {
		  free(image[i][j]);
		}
		free(image[i]);
	}
	free(image);
}

void createObject(char *label, Object *obj,float *position, float radius, float *surface_color,float reflectivity, float transparency,float *emission_color){
	for(int i=0; i<3; i++){
		obj->center[i] = position[i];
		obj->surface_color[i] = surface_color[i];
		obj->emission_color[i] = emission_color[i];
	}
	obj->label = label;
	obj->radius = radius;
	obj->radius2 = radius*radius;
	obj->transparency = transparency;
	obj->reflection = reflectivity;
}


int main(){
	Object *objects = malloc(6*sizeof(Object));
	//circle_center, radius, surface_color,reflectivity, transparency, emission_color
	float emission_color[3] = {0.0, 0.0, 0.0},ligth_ecolor[3] = {3.0 , 3.0, 3.0};
	float position1[3] = {0.0, -10004, -20}, surfacecolor1[3] = {0.20, 0.20, 0.20};
	char label1[] = "ground";
	createObject(label1, &objects[0],position1, 10000, surfacecolor1, 0, 0.0, emission_color);
	char label2[] = "red";
	float position2[3] = {0.0, 0, -20}, surfacecolor2[3] = {1.00, 0.32, 0.36};
	createObject(label2, &objects[1],position2, 4, surfacecolor2, 1, 0.5, emission_color); //red
	char label3[] = "yellow";
	float position3[3] = {5.0, -1, -15}, surfacecolor3[3] = {0.90, 0.76, 0.46};
	createObject(label3, &objects[2],position3, 2, surfacecolor3, 1, 0.0, emission_color); //yellow
	char label4[] = "blue";
	float position4[3] = {5.0, 0, -25}, surfacecolor4[3] = {0.65, 0.77, 0.97};
	createObject(label4,&objects[3],position4, 3, surfacecolor4, 1, 0.0, emission_color); //blue
	char label5[] = "black";
	float position5[3] = {-5.5, 0, -15}, surfacecolor5[3] = {0.90, 0.90, 0.90};
	createObject(label5,&objects[4],position5, 3, surfacecolor5, 1, 0.0, emission_color); //black
	// light
	char label6[] = "light";
	float position6[3] = {0.0, 20, -30}, surfacecolor6[3] = {0.00, 0.00, 0.00};
	createObject(label6,&objects[5],position6, 3, surfacecolor6, 0, 0, ligth_ecolor);
	render(objects, 6);
	return 0;
}
