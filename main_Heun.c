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
	fptr_s = fopen("20240509_demo_Heun.txt","w");

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


	for (int i=0; i<5000; i++){

		struct coord SxH;
		struct coord SxSxH;
		SxH = vector_product(S,H);
		SxSxH = vector_product( S , vector_product( S, H) );

		d_S.x = -gamma/(1+pow(lamda,2.0)) * (SxH.x + lamda*SxSxH.x );
		d_S.y = -gamma/(1+pow(lamda,2.0)) * (SxH.y + lamda*SxSxH.y );
		d_S.z = -gamma/(1+pow(lamda,2.0)) * (SxH.z + lamda*SxSxH.z );

		S_prime.x = S.x + d_S.x *dt;
		S_prime.y = S.y + d_S.y *dt;
		S_prime.z = S.z + d_S.z *dt;

		struct coord S_primexH;
		struct coord S_primexS_primexH;
		S_primexH = vector_product(S_prime,H);
		S_primexS_primexH = vector_product(S_prime, vector_product(S_prime, H));

		
		d_S_prime.x = -gamma/(1+pow(lamda,2.0)) * (S_primexH.x + lamda*S_primexS_primexH.x );
		d_S_prime.y = -gamma/(1+pow(lamda,2.0)) * (S_primexH.y + lamda*S_primexS_primexH.y );
		d_S_prime.z = -gamma/(1+pow(lamda,2.0)) * (S_primexH.z + lamda*S_primexS_primexH.z );
		
		S.x = S.x + 0.5*(d_S.x + d_S_prime.x)*dt;
		S.y = S.y + 0.5*(d_S.y + d_S_prime.y)*dt;
		S.z = S.z + 0.5*(d_S.z + d_S_prime.z)*dt;

		S_anal.x = (1/cosh(lamda*gamma*H.z*T/(1+pow(lamda,2))))*cos(gamma*H.z*T/(1+pow(lamda,2)));
		S_anal.y = (1/cosh(lamda*gamma*H.z*T/(1+pow(lamda,2))))*sin(gamma*H.z*T/(1+pow(lamda,2)));
		S_anal.z = sinh(lamda*gamma*H.z*T/(1+pow(lamda,2)))/cosh(lamda*gamma*H.z*T/(1+pow(lamda,2)));

		
		T = T + dt;

		fprintf(fptr_s, "%f, %f ,%f, %f, %f, %f, %f\n",T,S.x,S.y,S.z,S_anal.x-S.x,S_anal.y-S.y,S_anal.z-S.z);
	}
	fclose(fptr_s);
	return 0;
}
