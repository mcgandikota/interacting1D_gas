using namespace std;
#include <math.h>
#include <cstdio>
#include <random>
#include <iostream>

void initialize(int N);
void mdrun(void);
void forces(void);
void histogram(void);

float x[20000];
float force[20000];
float deltaT=0.0001;
int N;
float J=1.0;
float k;
int signK;
int skip=10000;
int print=1000;
float D=1.000;
float Gamma=1; //viscosity 
float T=1; //temperature
float bin_size=1.;
int Nbins=10000;  //even number only
int bin[20000];

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

FILE *movie;
movie=fopen("out.dump","w"); //rewwrites file
fclose(movie);

initialize(N);
generator.seed( rd() );

	for (int i=0; i<skip; i++){
	mdrun();
	}

	for(int i=0; i<steps; i++){
	mdrun();
		if(i%print==0){
			movie=fopen("out.dump","a");
			fprintf(movie,"%d\n\n",N);
			fclose(movie);
			for(int j=0; j<N; j++){
			movie=fopen("out.dump","a");
			fprintf(movie,"1 %f 0.0 0.0\n",x[j]);
			fclose(movie);
			//cout<<"1 "<<x[j]<<" 0.0 0.0"<<endl;
			}
		histogram();
		}
	}

			//for(int j=0; j<N-1; j++){
		//		cout <<x[j+1]-x[j]<<endl;
		//	}

FILE *out;
out=fopen("histogram.dat","w"); //rewwrites file
fclose(out);

	for(int i=0; i<(int)Nbins/2; i++){
	out=fopen("histogram.dat","a");
	fprintf(out,"%f %f \n",i*bin_size,bin[i]*1.0/(steps/print));
	fclose(out);
	}
	for(int i=(int)Nbins/2; i<Nbins; i++){
	out=fopen("histogram.dat","a");
	fprintf(out,"%f %f \n",-(i-Nbins/2.+1)*bin_size,bin[i]*1.0/(steps/print));
	fclose(out);
	}
		
return 0;
}

void initialize(int N){
float space=40;

	for (int i=0; i<N; i++){
	x[i] = -N/2 + 0.5*space + i*space;  //for even number N this gives a symmetric distribution about x=0
	}


	for (int i=0; i<Nbins; i++){ //initialize bins for histogram for position probability density
	bin[i]=0;
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
	interactions *= 1./2.*J*signK*k; //need -1 for k<0 and 1 for k>0
	//cout<<interactions<<endl;
	force[i]=harmonic+interactions;
	//cout<<force[i]<<endl;
	//if(force[i]!=force[i]) {cout<<"Blew up\n"<<endl; exit(0);}
	}
}

void histogram(void){
int n;

        for(int i=0;i<N;i++){
        n = (int)floor(abs(x[i])/bin_size);
		if (x[i]<0.0) {
		n += (int) Nbins/2.; 
		//cout<<"negative "<<n<<" "<<bin[n]+1<<endl;
		}
        bin[n]++;
        }
}
