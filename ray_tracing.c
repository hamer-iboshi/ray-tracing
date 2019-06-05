#include <stdio.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>
#define MAX_RAY_DEPTH 3

typedef struct {
	int direction[3];
} Ray;

typedef struct {
	double center[3]; //position of the sphere
	double radius, radius2; //sphere radius and radius^2
	double transparency, reflection; //surface transparency and reflectivity
	double surface_color[3], emission_color[3]; //surface color and emission (light)
} Object;


void subtractXYZ(double *result, double *p1, double *p2){
	result[0] = p1[0] - p2[0];
	result[1] = p1[1] - p2[1];
	result[2] = p1[2] - p2[2];
}

//functions of vec3
void normalize(double *vec){
	double nor2 = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]; 
	if (nor2 > 0) {
		double inv_nor = 1 / sqrt(nor2);
		vec[0] *= inv_nor;
		vec[2] *= inv_nor;
		vec[1] *= inv_nor;
	}
}

//functions of Sphere
bool intersect(Object *obj, double *rayorig, double *raydir, double *t0, double *t1){
	double l[3], tca, d2, thc;
	subtractXYZ(l,obj->center,rayorig);
	tca = raydir[0] * l[0] + raydir[1] * l[1] + raydir[2] * l[2]; //l.dot(raydir)
	if (tca < 0) return false;
	d2 = (l[0] * l[0] + l[1] * l[1] + l[2] * l[2]) - (tca * tca) ; //l.dot(l)
	if(d2 > obj->radius2) return false;
	thc = sqrt(obj->radius2 - d2);
	*t0 = tca - thc;
	*t1 = tca + thc;
	return true;
}



int Trace(Ray *ray, int depth){

	return 0;
}


int main(){
	/*
	int eyePosition[3] = {0,0,0};
	int lightPosition[3] = {0,0,0};
	int light_brightness = 0;

	for(int j=0; j < img_height; j++){
		for(int i=0; i < img_width; i++){
			//compute primary ray direction
			Ray primRay;
			computePrimRay(i, j, &primRay);
			Point pHit;
			Normal nHit;
			double minDist = INFINITY;
			Object object = NULL;

			//for each object if we use more than one object
			for (int k=0; k < n_objetcts; k++) {
				if (intersect(objects[k],primRay, &pHit, &nHit)) {
					double distance = Distance(eyePosition, pHit);
					if (distance < minDist) {
						object = objects[k];
						minDist = distance;
					}
				}
			}
			
			bool isShadow = false;				
			if (object != NULL) {
				//compute illumnination
				Ray shadowRay;
				subtractXYZ(&shadowRay.direction,&lightPosition,pHit); //direction = lightPosition - pHit
				for ( k = 0; k < n_objetcts; k++) {
					if (insersect(objects[k], shadowRay, &pHit, &nHit)) {
						isShadow = true;
						break;
					}
				}
			}
			if (!isShadow) {
				pixels[i][j] = object->color * light_brightness;
			} else {
				pixels[i][j] = 0;
			}
		}
	}

*/
	return 0;
}
