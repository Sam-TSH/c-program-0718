# include <stdio.h>
# include <math.h>


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

struct coord Get_k_LLG_equation(struct coord S_LLG,struct coord H_LLG, double gamma, double lamda){

  struct coord SxH = vector_product(S_LLG, H_LLG);
  struct coord SxSxH = vector_product(S_LLG, vector_product(S_LLG, H_LLG));

  struct coord k_LLG;
  k_LLG.x = -gamma / (1 + pow(lamda, 2.0)) * (SxH.x + lamda * SxSxH.x);
  k_LLG.y = -gamma / (1 + pow(lamda, 2.0)) * (SxH.y + lamda * SxSxH.y);
  k_LLG.z = -gamma / (1 + pow(lamda, 2.0)) * (SxH.z + lamda * SxSxH.z);

  return k_LLG;
}


int main(){
	
	FILE *fptr_s;
	fptr_s = fopen("20240718_demo_rk4.txt","w");

  double T = 0;
  double dt = 0.01;

  double lamda = 0.1;
  double gamma = 1;

  struct coord H = {0,0,1};
  struct coord S = {1,0,0};
  struct coord S_anal = {0,0,0};


  for (int i=0; i<5000; i++) {

    struct coord k1, k2, k3, k4;
    struct coord S_temp;

    k1 = Get_k_LLG_equation(S,H,gamma,lamda);

    S_temp.x = S.x + 0.5 * dt * k1.x;
    S_temp.y = S.y + 0.5 * dt * k1.y;
    S_temp.z = S.z + 0.5 * dt * k1.z;
    k2 = Get_k_LLG_equation(S_temp,H,gamma,lamda);

    S_temp.x = S.x + 0.5 * dt * k2.x;
    S_temp.y = S.y + 0.5 * dt * k2.y;
    S_temp.z = S.z + 0.5 * dt * k2.z;
    k3 = Get_k_LLG_equation(S_temp,H,gamma,lamda);

    S_temp.x = S.x + dt * k3.x;
    S_temp.y = S.y + dt * k3.y;
    S_temp.z = S.z + dt * k3.z;
    k4 = Get_k_LLG_equation(S_temp,H,gamma,lamda);

    S.x += dt * (k1.x + 2 * k2.x + 2 * k3.x + k4.x) / 6;
    S.y += dt * (k1.y + 2 * k2.y + 2 * k3.y + k4.y) / 6;
    S.z += dt * (k1.z + 2 * k2.z + 2 * k3.z + k4.z) / 6;

		S_anal.x = (1/cosh(lamda*gamma*H.z*T/(1+pow(lamda,2))))*cos(gamma*H.z*T/(1+pow(lamda,2)));
		S_anal.y = (1/cosh(lamda*gamma*H.z*T/(1+pow(lamda,2))))*sin(gamma*H.z*T/(1+pow(lamda,2)));
		S_anal.z = sinh(lamda*gamma*H.z*T/(1+pow(lamda,2)))/cosh(lamda*gamma*H.z*T/(1+pow(lamda,2)));

		
		T += dt;

		fprintf(fptr_s, "%f, %f ,%f, %f, %f, %f, %f\n",T,S.x,S.y,S.z,S_anal.x-S.x,S_anal.y-S.y,S_anal.z-S.z);
	}
	fclose(fptr_s);
	return 0;
}