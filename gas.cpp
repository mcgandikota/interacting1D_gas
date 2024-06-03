using namespace std;
#include <math.h>
#include <cstdio>
#include <random>
#include <iostream>

void initialize(int N);
void mdrun(void);
void forces(void);
float x[20000];
float force[20000];
float deltaT=0.0001;
int N;
float J=1.0;
float k;
int signK;
int skip=0;
int print=1000;
float D=1.000;
float Gamma=1; //viscosity 
float T=1; //temperature

random_device rd;
default_random_engine generator;
normal_distribution<double> distribution(0.0,2*Gamma*T);

int main(int argc, char *argv[]){
if (argc!=4) {printf("./a.out k N steps\n");exit(0);}
k = atof(argv[1]);
if (k>=0) {signK=1;}
else {signK=-1;}
N = atoi(argv[2]);
long int steps = atoi(argv[3]);

initialize(N);
generator.seed( rd() );

	for (int i=0; i<skip; i++){
	mdrun();
	}

	for(int i=0; i<steps; i++){
	mdrun();
		if(i%print==0){
		cout<<N<<endl<<endl;
			for(int j=0; j<N; j++){
			cout<<"1 "<<x[j]<<" 0.0 0.0"<<endl;
			}
		}
	}

			for(int j=0; j<N-1; j++){
				cout <<x[j+1]-x[j]<<endl;
			}

return 0;
}

void initialize(int N){

	for (int i=0; i<N; i++){
	x[i] = -N/2 + 0.5 + i;  //for even number N this gives a symmetric distribution about x=0
	}
}

void mdrun(void){
	forces();
	float random;

	for (int i=0; i<N; i++){
	random=distribution(generator);
	random=random/sqrt(deltaT);  
	random=random*sqrt(2*D);//coefficient??????????????
	x[i] = x[i] + (force[i]+random)*deltaT;
	}

}


void forces(void){
	float harmonic;
	float interactions; 
	int i,j;

	for (i=0; i<N; i++){
	harmonic = -x[i];
	interactions = 0;
	
		for (j=0;j<N;j++){
			if(j!=i){
			interactions += pow(abs(x[i]-x[j]),-k-2)*(x[i]-x[j]);
	//cout<<i<<" "<<j<<" "<<interactions<<endl;
			}
		}
	interactions *= -1./2.*J*signK*k; //need -1 for k<0 and 1 for k>0
	//cout<<interactions<<endl;
	force[i]=harmonic+interactions;
	//cout<<force[i]<<endl;
	//if(force[i]!=force[i]) {cout<<"Blew up\n"<<endl; exit(0);}
	}
}




