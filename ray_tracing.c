#include <stdio.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>

#define MAX_RAY_DEPTH 3
#define X 0
#define Y 1
#define Z 2
#define max(a,b) (a > b) ? a : b

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

double dot(double *v1, double *v2){
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}
//functions of vec3
void normalize(double *vec){
	double nor2 = dot(vec,vec); 
	if (nor2 > 0) {
		double inv_nor = 1 / sqrt(nor2);
		vec[0] *= inv_nor;
		vec[1] *= inv_nor;
		vec[2] *= inv_nor;
	}
}


//functions of Sphere
bool intersect(Object *obj, double *rayorig, double *raydir, double *t0, double *t1){
	double l[3], tca, d2, thc;
	subtractXYZ(l,obj->center,rayorig);
	tca = dot(l,raydir); //l.dot(raydir)
	if (tca < 0) return false;
	d2 = dot(l,l)- (tca * tca) ; //l.dot(l)
	if(d2 > obj->radius2) return false;
	thc = sqrt(obj->radius2 - d2);
	*t0 = tca - thc;
	*t1 = tca + thc;
	return true;
}


double mix(double a, double b, double mx){
	return b * mx + a * (1 - mx);
} 

/* This is the main trace function. It takes a ray as argument (defined by its origin and direction).Test if this ray intersects any of the geometry in the scene. 
If the ray intersects an object, compute the intersection point, the normal at the intersection point, and shade this point using this information. 
Shading depends on the surface property (is it transparent, reflective, diffuse). The function returns a color for the ray. 
If the ray intersects an object that is the color of the object at the intersection point, otherwise it returns the background color. */
double* trace(double *ray_orig, double *ray_dir, int depth, Object *obj,int obj_size){
	double tnear = INFINITY;
	int i;
	Object *object = NULL;
	static double color[3];
	// find intersection of this ray with the sphere in the scene
	for (i=0; i < obj_size; i++) {
		double t0 = INFINITY, t1 = INFINITY;
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
			return color;
	}
	double surfacecolor[3] = { 0.0, 0.0, 0.0 };
	double phit[3], nhit[3];
	for (i=0; i<3; i++)
		phit[i] = ray_orig[i] + ray_dir[i] * tnear;
	for (i=0; i<3; i++)
		nhit[i] = phit[i] - object->center[i];
	normalize(nhit); //normalize normal direction
	// If the normal and the view direction are not opposite to each other reverse the normal direction. 
	//That also means we are inside the sphere so set the inside bool to true. Finally reverse the sign of IdotN which we want positive.
	double bias = 1e-4; //add some vias to the poin to be traced
	bool inside = false;
	if (dot(ray_dir,nhit) > 0) {
		for (i=0; i<3; i++)
			nhit[i] = -nhit[i];
		inside = true;
	}
	double ray_orig2[3];
	if ((object->transparency > 0 || object->reflection > 0) && depth < MAX_RAY_DEPTH) {
		double facingratio = -dot(ray_dir,nhit);
		//change the mix value to tweak the effect
		double fresneleffect = mix(pow(1-facingratio,3),1,0.1);
		//compute reflection direction
		double refldir[3];
		double dot_raydir_nhit = dot(ray_dir,nhit);
		for(i=0; i<3; i++)
			refldir[i]-=nhit[i]*2*dot_raydir_nhit;
			ray_orig2[i] = phit[i] + nhit[i] * bias;
		double * reflection = trace(ray_orig2, refldir, depth+1, obj, obj_size);
		double * refraction = NULL;
		// if the sphere is also transparent compute refraction ray
		if (object->transparency) {
			double ior = 1.1, eta = (inside) ? ior : 1/ior; //inside the surface?
			double cosi = -dot(nhit,ray_dir);
			double k = 1 - eta*eta*(1 - cosi*cosi);
			double refrdir[3], ksqrt = sqrt(k);
			for(i=0; i<3; i++)
				refrdir[i] = ray_dir[i]*eta + nhit[i] * (eta * cosi - ksqrt);
				ray_orig2[i] = phit[i] + nhit[i] * bias;
			normalize(refrdir);
			refraction = trace(ray_orig2, refldir, depth+1, obj, obj_size);
		} 
			double refraction_aux[3] = { 0.0, 0.0, 0.0 };
			if(refraction == NULL) refraction = refraction_aux;
			for(i=0; i<3; i++)
				surfacecolor[i] = (reflection[i]*fresneleffect + refraction[i] * (1 - fresneleffect) * object->transparency) * object->surface_color[i];
	} else {
		//diffuse object
		for (i=0; i < obj_size; i++){
			if (obj->emission_color[X] > 0) {
				//light
				double transmission[3] = { 1.0, 1.0, 1.0 };
				double light_direction[3];
				for(int j=0; j<3; j++)
					light_direction[j] = obj[i].center[j] - phit[j];
				normalize(light_direction);
				for(int j=0; j < obj_size; j++){
					if( i != j){
						double t0, t1;
						for(int k=0; k<3; k++)
							ray_orig2[k] = phit[k] + nhit[k] * bias;
						if(intersect(&obj[i],ray_orig2,light_direction, &t0, &t1)){
							for(int k=0; k<3; k++)
								transmission[i] = 0;
							break;
						}
					}
				}
				for(int j=0; j<3; j++)
					surfacecolor[i] += object->surface_color[i] * transmission[i] * max(0.0,dot(nhit,light_direction)) * object->emission_color[i];
			}
		}
		
	}
	for(i=0; i<3; i++)
		color[i] = surfacecolor[i] + object->emission_color[i];
	return color;
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
					if (intersect(objects[k], shadowRay, &pHit, &nHit)) {
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
