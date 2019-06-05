#include <stdio.h>
#include <limits.h>

int main(){
	
	int eyePosition[3] = {0,0,0};
	int lightPosition[3] = {0,0,0};


	for(int j=0; j < img_height; j++){
		for(int i=0; i < img_width; i++){
			//compute primary ray direction
			Ray primRay;
			computePrimRay(i, j, &primRay);
			Point pHit;
			Normal nHit;
			double minDist = DBL_MAX;
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
				subtractXYZ(shadowRay.direction,lightPosition,pHit); //direction = lightPosition - pHit
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


	return 0;
}
