# include <stdio.h>
# include <math.h>


double T = 0;
double dt = 0.01;

struct coord{
	double x;
	double y;
	double z;
};

struct coord vector_product(struct coord v1,struct coord v2){
	struct coord result;

	result.x = (v1.y * v2.z) - (v1.z * v2.y);
	result.y = (v1.z * v2.x) - (v1.x * v2.z);
	result.z = (v1.x * v2.y) - (v1.y * v2.x);

	return result;
}

int main(){
	
	FILE *fptr_s;
	fptr_s = fopen("20240509_demo.txt","w");

	double lamda = 0.1;
	double gamma = 1;
	double T = 0;
	double dt = 0.01;

	struct coord H = {0,0,1};
	struct coord S = {1,0,0};
	struct coord S_prime = {0,0,0};
	struct coord d_S = {0,0,0};
	struct coord d_S_prime = {0,0,0};

	struct coord S_anal = {0,0,0};


for (int i=0; i<5000; i++) {
    struct coord k1, k2, k3, k4;
    struct coord S_temp;

    struct coord SxH = vector_product(S, H);
    struct coord SxSxH = vector_product(S, vector_product(S, H));

    k1.x = dt * (-gamma / (1 + pow(lamda, 2.0)) * (SxH.x + lamda * SxSxH.x));
    k1.y = dt * (-gamma / (1 + pow(lamda, 2.0)) * (SxH.y + lamda * SxSxH.y));
    k1.z = dt * (-gamma / (1 + pow(lamda, 2.0)) * (SxH.z + lamda * SxSxH.z));

    S_temp.x = S.x + 0.5 * k1.x;
    S_temp.y = S.y + 0.5 * k1.y;
    S_temp.z = S.z + 0.5 * k1.z;

    SxH = vector_product(S_temp, H);
    SxSxH = vector_product(S_temp, vector_product(S_temp, H));

    k2.x = dt * (-gamma / (1 + pow(lamda, 2.0)) * (SxH.x + lamda * SxSxH.x));
    k2.y = dt * (-gamma / (1 + pow(lamda, 2.0)) * (SxH.y + lamda * SxSxH.y));
    k2.z = dt * (-gamma / (1 + pow(lamda, 2.0)) * (SxH.z + lamda * SxSxH.z));

    S_temp.x = S.x + 0.5 * k2.x;
    S_temp.y = S.y + 0.5 * k2.y;
    S_temp.z = S.z + 0.5 * k2.z;

    SxH = vector_product(S_temp, H);
    SxSxH = vector_product(S_temp, vector_product(S_temp, H));

    k3.x = dt * (-gamma / (1 + pow(lamda, 2.0)) * (SxH.x + lamda * SxSxH.x));
    k3.y = dt * (-gamma / (1 + pow(lamda, 2.0)) * (SxH.y + lamda * SxSxH.y));
    k3.z = dt * (-gamma / (1 + pow(lamda, 2.0)) * (SxH.z + lamda * SxSxH.z));

    S_temp.x = S.x + k3.x;
    S_temp.y = S.y + k3.y;
    S_temp.z = S.z + k3.z;

    SxH = vector_product(S_temp, H);
    SxSxH = vector_product(S_temp, vector_product(S_temp, H));

    k4.x = dt * (-gamma / (1 + pow(lamda, 2.0)) * (SxH.x + lamda * SxSxH.x));
    k4.y = dt * (-gamma / (1 + pow(lamda, 2.0)) * (SxH.y + lamda * SxSxH.y));
    k4.z = dt * (-gamma / (1 + pow(lamda, 2.0)) * (SxH.z + lamda * SxSxH.z));

    S.x += (k1.x + 2 * k2.x + 2 * k3.x + k4.x) / 6;
    S.y += (k1.y + 2 * k2.y + 2 * k3.y + k4.y) / 6;
    S.z += (k1.z + 2 * k2.z + 2 * k3.z + k4.z) / 6;

		S_anal.x = (1/cosh(lamda*gamma*H.z*T/(1+pow(lamda,2))))*cos(gamma*H.z*T/(1+pow(lamda,2)));
		S_anal.y = (1/cosh(lamda*gamma*H.z*T/(1+pow(lamda,2))))*sin(gamma*H.z*T/(1+pow(lamda,2)));
		S_anal.z = sinh(lamda*gamma*H.z*T/(1+pow(lamda,2)))/cosh(lamda*gamma*H.z*T/(1+pow(lamda,2)));

		
		T = T + dt;

		fprintf(fptr_s, "%f, %f ,%f, %f, %f, %f, %f\n",T,S.x,S.y,S.z,S_anal.x-S.x,S_anal.y-S.y,S_anal.z-S.z);
	}
	fclose(fptr_s);
	return 0;
}
