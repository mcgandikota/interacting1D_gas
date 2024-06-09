using namespace std;
#include <math.h>
#include <cstdio>
#include <random>
#include <iostream>

void initialize(int N);
void mdrun(void);
void forces(void);
void histogram(void);
void energy(void);

float x[20000];
float velocity[20000];
float force[20000];
float deltaT=0.00001;
int N;
float J=1.0;
float k;
int signK;
//int skip=3000000;
int skip=400000;
int print=1000;
float Gamma=1; //viscosity 
float T=0.50; //temperature
float bin_size=0.4;
int Nbins=100000;  //even number only
int bin[20000000];
float potential_energy;
float harmonic_energy;
float interaction_energy;
float kinetic_energy;

random_device rd;
default_random_engine generator;
normal_distribution<double> distribution(0.0,sqrt(Gamma*T)); //coefficient??????????????, you are assuming equilibrium here

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

FILE *energy_out;
energy_out=fopen("energy.dat","w"); //rewwrites file
fclose(energy_out);

initialize(N);
generator.seed( rd() );

	cout<<"Equilibration start"<<endl;
	for (int i=0; i<skip; i++){
	mdrun();
	}
	cout<<"Equilibration end"<<endl;

	for(int i=0; i<steps; i++){
	mdrun();
		if(i%print==0){
		cout<<i<<" ";
		//print movie
		movie=fopen("out.dump","a");
		fprintf(movie,"%d\n\n",N);
		fclose(movie);
			for(int j=0; j<N; j++){
			movie=fopen("out.dump","a");
			fprintf(movie,"1 %f 0.0 0.0\n",x[j]);
			fclose(movie);
			}

		//print energies
		energy();
		energy_out=fopen("energy.dat","a"); //rewwrites file
		fprintf(energy_out,"%f %f %f\n",potential_energy,kinetic_energy,potential_energy+kinetic_energy);
		fclose(energy_out);
		
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
	fprintf(out,"%f %f \n",i*bin_size/N,bin[i]*1.0/(steps/print));
	fclose(out);
	}
	for(int i=(int)Nbins/2; i<Nbins; i++){
	out=fopen("histogram.dat","a");
	fprintf(out,"%f %f \n",-(i-Nbins/2.+1)*bin_size/N,bin[i]*1.0/(steps/print));
	fclose(out);
	}
		
cout<<endl;
return 0;
}

void initialize(int N){
float space=0.3;

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
	velocity[i]=force[i]+random;
	random=random/sqrt(deltaT);  
	x[i] = x[i] + (force[i]+random)*deltaT;
	}

}


void forces(void){
	float harmonic;
	float interactions; 
	int i,j;

	for (i=0; i<N; i++){
	harmonic = -x[i];
	interactions = 0.;
	
		for (j=0;j<N;j++){
			if(j!=i){
			interactions += pow(abs(x[i]-x[j]),-k-2)*(x[i]-x[j]);
			}
		}
	interactions *= 1./2.*J*signK*k; //need -1 for k<0 and 1 for k>0
	force[i]=harmonic+interactions;
	if(force[i]!=force[i]) {cout<<"Blew up\n"<<endl; exit(0);}
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

void energy(void){
//calculate energy 
potential_energy=0.;
harmonic_energy=0.;
interaction_energy=0.;
kinetic_energy=0.;

        for(int i=0;i<N;i++){
	kinetic_energy += velocity[i]*velocity[i];
	}
	kinetic_energy *= 1/2.;

        for(int i=0;i<N;i++){
	harmonic_energy += x[i]*x[i];
        	for(int j=0;j<N;j++){
			if(j!=i){
			interaction_energy += pow(abs(x[i]-x[j]),-k);
			}
		}

	}
harmonic_energy *= 1/2.;
interaction_energy *= J*signK/2.;
potential_energy=harmonic_energy+interaction_energy;

kinetic_energy/= N;
potential_energy/= N;
}
