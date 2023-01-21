/* Numerical simulation of the 2D Inertial Spin Model
Davide Venturelli & Emanuele Loffredo
Last update: 25.04.2022
Compile as $ gcc -o ism ism.c -lm -Wall -Wextra -O3 -march=native
Execute as ./ism <temperature> <flock speed> <inertia>
*/

// SIMULATION PARAMETERS
#define L 20				// # of individuals per side (i.e. unit density)
#define N (L*L)
#define rc 1.5
#define NP 9				// # of parameters
#define cost 1/(2*sqrt(3))
#define N_ext 1				// # of individuals subject to perturbation
#define t0 0				// time at which the perturbation begins
#define tau 0.1				// duration of perturbation
#define A0 0				// amplitude of perturbation
#define n_c 7				// # of nearest-neighbours for topological interaction
#define T_term	100			// thermalization time

// LIBRARIES & TYPES
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_sf.h>
//#include <gsl/gsl_eigen.h>

typedef struct{
	double tmax;
	double dt;
	double eta;
	double chi;
	double T;
	double J;
	} input; 

typedef struct{
	double x;
	double y;
	int tile;				// Position on the mosaic
	double delta;			// Temporary variable for storing distances
	} vector;
	
typedef struct{
	int *list;				// Bucket list of the birds it contains
	int numData;			// Number of items it contains
	int count;				// How many steps before the next update?
	} cell;

// GLOBAL VARIABLES	
// Variables: r(t) position, o(t) polar angle, u(t) angular velocity.
// Variables that have to be updated by external functions at every step are local, but their pointers are global.
double *op, *Op, *up, *Up;
vector *pos;
input *parameter;
double gauss[2*N], An[N], fn[N];
//double *adjacencym;
double fn1;
double v0;					// Flock speed
int near[24*N];				// List of nearest neighbours
int * windx, * windy;		// Winding numbers
vector * initial;			// Stores initial positions
int p=(L+1)*L/2;			// Index of the bird in the center of the box

// FUNCTION PROTOTYPES
void boxmuller(double *address);
void neighborhood(int *neighbours);
int x(int k);  										// gives x coord of k-th particle	
int y(int k);  										// gives y coord of k-th particle
int indice(int x, int y);  							// gives i-th particle from coordinates
void updatemosaic(cell *mosaic);
void fill(cell *mosaic, int i);
void empty(cell *mosaic);
void startlattice();
void startrandom();
void init();
void initUP();
void inithalf();
void Swap(double **a, double **b);
void evolve(double *factor, cell *mosaic);
void evolve_field(double *factor, cell *mosaic, double t, double shift);
vector polarization();
void observable(double tfin, double factor[], cell *mosaic, cell *localbird);
void trajectory(double tfin, double factor[], cell *mosaic);
void turn_angle(double tfin, double factor[], cell *mosaic, double angle);
double meanvic(cell *mosaic);													// mean number of neighbours
void wandering(double tfin, double factor[], cell *mosaic);
double dist(int i, int j);
void evolve_field_topological(double *factor, cell *mosaic, cell * bird, double t, double shift);
void quicksort(int * data, int first, int last);
void update_sortedlist(double *factor, cell * mosaic, cell * bird);
void printCells(cell * mosaic);
void help();
void update_wind();
double reshift(double posit);
void find_p(cell * mosaic, double v);

// EXTERNAL PERTURBATION
double H_1 (double t) {return (A0 /2)*(1+tanh((t-t0)/tau));} 					// field coupled with velocity (angle) to induce a turn
//double H_1 (double t) {return 0;} 											// use this to study wandering time

// MAIN BODY
int main(int argc, char** argv){
	if (argc != 4) help();

	// Seed of random function
	int seed, i;
	seed = 200*time(NULL);
	srand48 (seed);
	
	// Local variables for trajectories
	double o[N];
	double O[N];
	double u[N] = {0};
	double U[N];
	vector r[N];
	input par[6];
	int w_x[N];														// Winding # along x
	int w_y[N];														// Winding # along y
	vector init_pos[N]; 											// Stores initial positions
	
	// Assignation to their global pointers
	op = o;
	Op = O;
	up = u;
	Up = U;
	pos = r;
	parameter = par;
	windx = w_x;
	windy = w_y;
	initial = init_pos;
	
	// Input acquisition and one time operations
	parameter->tmax=1000;
	parameter->dt=0.01;
	parameter->eta=1;
	parameter->chi=atof(argv[3]);
	parameter->T=atof(argv[1]);	
	parameter->J=8;
	
	double tfin=parameter->tmax;
	double factor[NP];
	v0=atof(argv[2]);																	// read flock speed from input					
	
	factor[0]=parameter->dt;															// dt 
	factor[1]=0.5*parameter->dt * parameter->dt;										// dt^2/2
	factor[2]=factor[5]*parameter->dt;													// sigma*dt*sqrt(dt)		
	factor[3]=0.5*parameter->dt;														// dt/2
	factor[4]=parameter->dt*parameter->eta/parameter->chi; 								// 2*gamma*dt
    factor[5]=sqrt(2*parameter->T*parameter->eta)/(parameter->chi)*sqrt(parameter->dt);	// sigma*sqrt(dt)
    factor[6]=parameter->eta/parameter->chi;											// 2*gamma
    factor[7]=parameter->J/parameter->chi;												// J/chi
    factor[8]=1/parameter->chi;															// beta
	
	// Initialization
	initUP();
	//startrandom();														
	startlattice();
	neighborhood(near);					
	
	for(i=0; i<N; i++){																	// Stores initial positions
		(initial+i)->x = (pos+i)->x;
		(initial+i)->y = (pos+i)->y;
	}						
	
	cell mosaic[N];																		
	for (i=0;i<N;i++){																	// Allocates mosaic
		mosaic[i].list=(int *)malloc(L*L*sizeof(int));
	}
	updatemosaic(mosaic);																// Initializes mosaic
	
	cell localbird[N];																	// Ordered list of neighbouring birds
	for (i=0;i<N;i++){																	// Allocates sorted list
		localbird[i].list = (int *)malloc(N*sizeof(int));
		localbird[i].count = 0;
	}
	update_sortedlist(factor, mosaic, localbird);										// Initializes ordered list
	
	printf("Prima %d\n",p);

	// Evolution		
	for (i=0;i<(int)(T_term/parameter->dt);i++){										//Initial transient to discard
		evolve(factor, mosaic);
		//evolve_field_topological(factor, mosaic, localbird, 0,0);
	}

	find_p(mosaic, v0);
	printf("Dopo %d\n",p);
		
	//observable(tfin, factor, mosaic, localbird);
	//trajectory(tfin, factor, mosaic);
	//wandering(tfin, factor, mosaic);
	//turn_angle(tfin, factor, mosaic, M_PI/2);
	
	// Free memory and exit
	for (i=0;i<N;i++){
		free(mosaic[i].list);
		free(localbird[i].list);
	}
	exit(EXIT_SUCCESS);
}


// DEFINITION OF FUNCTIONS
// Creates the nearest-neighbour list (up to the second)
void neighborhood(int *neighbours){
	int jxr, jyu, jxl, jyd, k;
	int jxrr, jyuu, jxll, jydd;
	
	for (k=0;k<N;k++){
		
		// fillig with nearest neighbours
		jxr=(x(k)+1)%L;
		neighbours[24*k]=indice(jxr,y(k));
		jxl=(x(k)+L-1)%L;
		neighbours[24*k+1]=indice(jxl,y(k));
		jyu=(y(k)+1)%L;
		neighbours[24*k+2]=indice(x(k),jyu);
		jyd=(y(k)+L-1)%L;
		neighbours[24*k+3]=indice(x(k),jyd);
		// filling with next nearest neighbours
		neighbours[24*k+4]=indice(jxr,jyu);
		neighbours[24*k+5]=indice(jxr,jyd);
		neighbours[24*k+6]=indice(jxl,jyd);
		neighbours[24*k+7]=indice(jxl,jyu);
		// filling rest of the 5x5 cube
		jxrr=(x(k)+1+1)%L;
		jxll=(x(k)+L-1-1)%L;
		jyuu=(y(k)+1+1)%L;
		jydd=(y(k)+L-1-1)%L;
		
		neighbours[24*k+8]=indice(x(k),jyuu);
		neighbours[24*k+9]=indice(jxr,jyuu);
		neighbours[24*k+10]=indice(jxrr,jyuu);
		neighbours[24*k+11]=indice(jxrr,jyu);
		neighbours[24*k+12]=indice(jxrr,y(k));
		neighbours[24*k+13]=indice(jxrr,jyd);
		neighbours[24*k+14]=indice(jxrr,jydd);
		neighbours[24*k+15]=indice(jxr,jydd);
		neighbours[24*k+16]=indice(x(k),jydd);
		neighbours[24*k+17]=indice(jxl,jydd);
		neighbours[24*k+18]=indice(jxll,jydd);
		neighbours[24*k+19]=indice(jxll,jyd);
		neighbours[24*k+20]=indice(jxll,y(k));
		neighbours[24*k+21]=indice(jxll,jyu);
		neighbours[24*k+22]=indice(jxll,jyuu);
		neighbours[24*k+23]=indice(jxl,jyuu);
	}
}	

// Bucket list handling
void fill(cell *mosaic, int i){
	mosaic->list[(mosaic->numData)++]=i;}
	
void empty(cell *mosaic){
	mosaic->numData=0;}

// Construction/update of the mosaic
void updatemosaic(cell *mosaic){
	int i,j;
	for (j=0;j<N;j++){
		empty(mosaic +j);
		}
	
	for (i=0;i<N;i++){
		//j=indice((int)(fmod((((pos+i)->x)+L),L)),(int)(fmod((((pos+i)->y)+L),L)));
		j=indice((int)((pos+i)->x),(int)((pos+i)->y));
		fill(mosaic+j,i);
		(pos+i)->tile=j;							// assignes particle i to j-th cell of mosaic
		}
		
	update_wind();									// Triggers update of winding numbers
}

// Updates the winding numbers
void update_wind(){
	int i;
	for(i=0; i<N; i++){
		windx[i] = (int) ( ((pos+i)->x - (initial+i)->x )/L );
		windy[i] = (int) ( ((pos+i)->y - (initial+i)->y )/L );
	}
}

// Adaptation of sorting algorithm (p.417 Barone et al.) to our case.
// Recursive function. In the first call, set first==0 and last==(n-1), where n is the length of the list to sort.
void quicksort(int * data, int first, int last){
	double cut, tmp;
	int i,j;
	
	if(first < last){								// Checks if the vector contains at least 2 elements (if not, do nothing)
		cut = (pos+data[last])->delta;				// Set delta of the last element as upper threshold
		i = first;
		j = last;
		
		while(i <= j){
			// The two indices move towards each other until they meet. 
			// The vector is sorted so that all the elements <= cut are on the left, the others on the right.
			while( (i <= last) && ( (pos+data[i])->delta <= cut ) ) i++;
			while( (pos+data[j])->delta > cut ) j--;
			
			if( i < j ){							// Exchange the two positions
				tmp = data[i];
				data[i] = data[j];
				data[j] = tmp;
				i++;
				j--;
			}
		}
		
		if( i<= last ){								// Split the vector into 2 subvectors and sort them
			quicksort(data, first, j);
			quicksort(data, i, last);
		} else {
			quicksort(data, first, last-1);
		}
	}
}

// Updates the ordered list of nearest birds (topological interaction only)
void update_sortedlist(double *factor, cell *mosaic, cell * bird){
	int i, j, k, index;
	
	for(i=0; i<N; i++){										// Cycle over all birds
		if( bird[i].count > 0 ){							// Skip the update if it's not time yet
			(bird[i].count)--;		
			continue;
		}													// Update the list otherwise
		
		// Filling the list. I want it as long as (n_c+2): one is always me, and we need the first outer element.
		empty( (bird+i) );									// Empty the list of bird i
			
		k=0;												// Cycle over same tile (self-interaction is null)
		while( k < mosaic[((pos+i)->tile)].numData ){
			index = mosaic[((pos+i)->tile)].list[k];		// Copy the index of the bird just found
			k++;
			
			if( windx[i] != windx[index] || windy[i] != windy[index]) continue; // Check for open BC
			
			fill( bird+i, index );							// Add it to the sorted list
			(pos+index)->delta = dist(i,index);				// Saves in bird index its distance from i-th 
		}
		
		if( bird[i].numData < (n_c+2) ){					// Cycle over the square
			for (j=0;j<8;j++){
				k=0;
				while( k < mosaic[near[24*((pos+i)->tile) + j]].numData ){
					index = mosaic[near[24*((pos+i)->tile) + j]].list[k];
					k++;
					
					if( windx[i] != windx[index] || windy[i] != windy[index]) continue; // Check for open BC
					
					fill( bird+i, index );				
					(pos+index)->delta = dist(i,index);	
				}
			}
		}
		
		if( bird[i].numData < (n_c+2) ){					// Cycle over the outer crown (hopefully unneeded)
			for (j=8;j<24;j++){
				k=0;
				while( k < mosaic[near[24*((pos+i)->tile) + j]].numData ){
					index = mosaic[near[24*((pos+i)->tile) + j]].list[k];
					k++;
					
					if( windx[i] != windx[index] || windy[i] != windy[index]) continue; // Check for open BC
					
					fill( bird+i, index );				
					(pos+index)->delta = dist(i,index);	
				}
			}
		}
		
		if( bird[i].numData < (n_c+2) ){	
			printf("Your flock got too diluted, I can't correctly account for the interactions anymore.\n");
			exit(EXIT_FAILURE);
		}
		
		// Now that we have a list with (at least) n_c+2 birds, we need to sort them by their distances from i
		quicksort( bird[i].list, 0, n_c+1 );				// It doesn't sort the whole list, it's enough up to n_c+2
		
		// Finally, we give an estimate of the time this will remain valid and save it in the variable count.
		// Worst case scenario, the two birds on the border cross after count = |dist(i,n_c+1)-dist(i,n_c+2)| / (2*v0*dt) timesteps.
		bird[i].count = fabs( (pos+ bird[i].list[n_c] )->delta - (pos+ bird[i].list[n_c+1] )->delta ) / (2*v0*factor[0]);
	}
}

// Time evolution in the presence of external perturbing fields, with topological interactions.
// "shift" is the phaseshift added to the angle. The applied field needs to know "t", i.e. the physical time (note evolve() would be otherwise autonomous).
// Can be used without a field, just set it to zero and call evolve_field_topological(factor, mosaic, bird, 0, 0)
void evolve_field_topological(double *factor, cell *mosaic, cell * bird, double t, double shift){	
	int i, j, k, vic;
	double xn, xn1;
	//double aux; // PBC position update
	
	// Update angles
	for (i=0;i<N;i++){
		// Evaluates fn[i]
		boxmuller(gauss +2*i);
		fn[i]=0;
		//printf("# With bird %d : \t",i);//DEBUG
		for (k=1; k<(n_c+1); k++){								// Skip k=0 which is the bird itself
			vic = bird[i].list[k];
			//printf("%d\t", vic);//DEBUG
			fn[i] += sin( *( op+ vic ) - *(op+i));	
		}
		fn[i] *= factor[7];
		//printf("\n");//DEBUG
		
		// Evaluates An[i]
		An[i] = factor[1]*(fn[i]-factor[6]*(*(up+i))) + factor[2]*(0.5*gauss[2*i]+cost*gauss[2*i+1]);
		// Evaluates xn+1[i]
		*(Op+i) = *(op+i) + An[i] + factor[0]*( *(up+i) );
	}
	
	// Apply perturbation to the angles
	xn = *(op+ p) - shift;
	*(Op+p) += -factor[3]*(factor[8]*factor[0]*(H_1)(t)*sin(xn));
	
	// Update positions (Euler-Cromer)
	int mos_flag=0;												// Time to update the mosaic
	for (i=0; i<N;i++){
		(pos+i)->x = (pos+i)->x + v0*factor[0]*cos(*(Op+i));
		//(pos+i)->x = fmod((aux+L),L);							// Applying PBC
		(pos+i)->y = (pos+i)->y + v0*factor[0]*sin(*(Op+i));
		//(pos+i)->y = fmod((aux+L),L);							// fmod() is the func % with double
		
		// This can be optimized, because the error in constructing the mosaic is null, while actually it is enough that sortedlist be correct.
		if(!mos_flag){											// Flags as soon as a particle changes tile
			j=indice((int)((pos+i)->x),(int)((pos+i)->y));
			if( ((pos+i)->tile) != j ) mos_flag=1;
		}
	}
	
	if(mos_flag) {
		updatemosaic(mosaic);	
	}
	update_sortedlist(factor, mosaic, bird);	
	
	// Update angular velocities
	for (i=0;i<N;i++){	
		//Evaluates fn+1 to compute un+1[i]
		fn1=0;
		for (k=1; k<(n_c+1); k++){
			vic = bird[i].list[k];
			fn1 += sin( *(Op+vic) - *(Op+i) );
		}	
		fn1 *= factor[7];
		
		//Evaluates u_n+1[i]
		*(Up+i) = (1-factor[4])*(*(up+i)) + factor[3]*(fn1 + fn[i]) + factor[5]*gauss[2*i] - factor[6]*An[i];
	}
	
	// Apply perturbation to angular velocities (present even if the field is coupled to the angles only)
	xn = *(op + p) - shift;
	xn1 = *(Op + p) - shift;
	*(Up+p) += factor[3]*((1 - factor[4])*(-factor[8]*(H_1)(t)*sin(xn)) - factor[8]*(H_1)(t + factor[0])*sin(xn1)); 
	
	// Sends x_(n+1) into x_n
	Swap(&op, &Op);
	Swap(&up, &Up);
}

// Evolution with metric interaction
void evolve(double *factor, cell *mosaic){	
	int i, j, k, vic;
	for (i=0;i<N;i++){
		//Evaluates fn[i]
		boxmuller(gauss +2*i);
		fn[i]=0;
		
		for (j=0;j<24;j++){
			k=0;
			// cylcling over all the 5x5 crown
			while (k<mosaic[near[24*((pos+i)->tile) + j]].numData) {
				vic = mosaic[near[24*((pos+i)->tile) +j]].list[k];
				if(dist(i,vic)<=rc && windx[i]==windx[vic] && windy[i]==windy[vic]) {
					fn[i] += sin(*(op+vic) - *(op+i));
				}	
				k++;
			}
		}
		
		k=0;
		while (k<mosaic[((pos+i)->tile)].numData) {
			vic = mosaic[((pos+i)->tile)].list[k];
			if(windx[i]==windx[vic] && windy[i]==windy[vic]){
				fn[i] +=sin(*(op+vic) - *(op+i));	
			}
			k++;
		}
		
		fn[i] *= factor[7];
		
		// Evaluates An[i]
		An[i] = factor[1]*(fn[i]-factor[6]*(*(up+i))) + factor[2]*(0.5*gauss[2*i]+cost*gauss[2*i+1]);
		
		//Evaluates xn+1[i]
		*(Op+i)=*(op+i) +An[i]+factor[0]*(*(up+i));
	}
	
	for (i=0; i<N;i++){
		// Update positions
		(pos+i)->x= (pos+i)->x + v0*factor[0]*cos(*(Op+i));
		//(pos+i)->x = fmod((aux+L),L);							// applying PBC
		(pos+i)->y= (pos+i)->y + v0*factor[0]*sin(*(Op+i));
		//(pos+i)->y = fmod((aux+L),L);							// fmod() is the func % with double
	}
	
	updatemosaic(mosaic);									
	
	for (i=0;i<N;i++){	
		// Evaluates fn+1 to compute un+1[i]
		fn1=0;
		for (j=0;j<24;j++){	
			k=0;
			// cylcling over all the 5x5 crown
			while (k<mosaic[near[24*((pos+i)->tile) + j]].numData) {	
				vic = mosaic[near[24*((pos+i)->tile) +j]].list[k];
				if(dist(i,vic)<=rc && windx[i]==windx[vic] && windy[i]==windy[vic]) {
					fn1 +=sin(*(Op+vic) - *(Op+i));}
				k++;			
			}
		}
		
		k=0;
		while (k<mosaic[((pos+i)->tile)].numData) {
			vic = mosaic[((pos+i)->tile)].list[k];
			if(windx[i]==windx[vic] && windy[i]==windy[vic]){
				fn1 += sin(*(Op+vic) - *(Op+i));}
			k++;
		}
		fn1 *= factor[7];
		
		// Evaluates un+1[i]
		*(Up+i) = (1-factor[4])*(*(up+i)) + factor[3]*(fn1 + fn[i]) + factor[5]*gauss[2*i] - factor[6]*An[i];
	}
	
	Swap (&op, &Op);
	Swap (&up, &Up);
}

// Evolution in the presence of an external perturbing field, with metric interactions
void evolve_field(double *factor, cell *mosaic, double t, double shift){	
	int i, j, k, vic;
	double xn, xn1;
	for (i=0;i<N;i++){
		// Evaluates fn[i]
		boxmuller(gauss +2*i);
		fn[i]=0;
		for (j=0;j<24;j++){
			k=0;
			// cylcling over all the 5x5 crown
			while (k<mosaic[near[24*((pos+i)->tile) + j]].numData) {
				vic=mosaic[near[24*((pos+i)->tile) +j]].list[k];
				if(dist(i,vic)<=rc && windx[i]==windx[vic] && windy[i]==windy[vic]) {
					fn[i] +=sin(*(op+vic) - *(op+i));}
				k++;
			}
		}
		
		k=0;
		while (k<mosaic[((pos+i)->tile)].numData) {
			vic=mosaic[((pos+i)->tile)].list[k];
			if(windx[i]==windx[vic] && windy[i]==windy[vic]){
				fn[i] +=sin(*(op+vic) - *(op+i));}
			k++;
		}
		
		fn[i] *= factor[7];
		
		// Evaluates An[i]
		An[i] = factor[1]*(fn[i]-factor[6]*(*(up+i))) + factor[2]*(0.5*gauss[2*i]+cost*gauss[2*i+1]);
		
		// Evaluates xn+1[i]
		*(Op+i)=*(op+i) +An[i]+factor[0]*(*(up+i));
	}
	
	//Apply perturbation
	xn = *(op+ p) - shift;
	*(Op+p) += -factor[3]*(factor[8]*factor[0]*(H_1)(t)*sin(xn));
		
	for (i=0; i<N;i++){
		//Update positions
		(pos+i)->x= (pos+i)->x + v0*factor[0]*cos(*(Op+i));
		//(pos+i)->x = fmod((aux+L),L);							// applying PBC
		(pos+i)->y= (pos+i)->y + v0*factor[0]*sin(*(Op+i));
		//(pos+i)->y = fmod((aux+L),L);							// fmod() is the func % with double
	}
	
	updatemosaic(mosaic);			
	
	for (i=0;i<N;i++){	
		// Evaluates fn+1 to compute un+1[i]
		fn1=0;
		for (j=0;j<24;j++){	
			k=0;
			// cylcling over all the 5x5 crown
			while (k<mosaic[near[24*((pos+i)->tile) + j]].numData) {
				vic=mosaic[near[24*((pos+i)->tile) +j]].list[k];
				if(dist(i,vic)<=rc && windx[i]==windx[vic] && windy[i]==windy[vic]) {
					fn1 +=sin(*(Op+vic) - *(Op+i));}
				k++;			
			}
		}
		
		k=0;
		while (k<mosaic[((pos+i)->tile)].numData) {
			vic=mosaic[((pos+i)->tile)].list[k];
			if(windx[i]==windx[vic] && windy[i]==windy[vic]){
			fn1 +=sin(*(Op+vic) - *(Op+i));}
			k++;
		}
		
		fn1 *= factor[7];
		
		// Evaluates un+1[i]
		*(Up+i) = (1-factor[4])*(*(up+i)) + factor[3]*(fn1 + fn[i]) + factor[5]*gauss[2*i] - factor[6]*An[i];
	}
	
	// Apply perturbation
	xn = *(op + p) - shift;
	xn1 = *(Op + p) - shift;
	*(Up+p) += factor[3]*((1 - factor[4])*(-factor[8]*(H_1)(t)*sin(xn)) - factor[8]*(H_1)(t + factor[0])*sin(xn1)); 
	
	Swap (&op, &Op);
	Swap (&up, &Up);
}

// Lattice initialization
void startlattice(){			// position on the vertices of a 2d regular lattice
	int i;
	for (i=0;i<N;i++){
		(pos+i)->x=x(i);
		(pos+i)->y=y(i);
	}
}

void startrandom(){				// position randomly scattered in the 2d box
	int i;
	for (i=0;i<N;i++){
		(pos+i)->x = L*drand48();
		(pos+i)->y = L*drand48();
	}
}

void init(){ 					// initializes all velocities at random
	int i;
	for (i=0;i<N;i++){
		*(op+i)=2*M_PI*drand48();
	}	
}

void initUP(){					// initializes all velocities up
	int i;
	for (i=0;i<N;i++){
		*(op+i)=1.5*M_PI;
	}	
}

void inithalf(){				// initializes half velocities up
	int i;
	for (i=0;i<N;i++){
		*(op+i)=0;
		if (i<=(int)(N/80)) *(op+i)=2*M_PI*drand48();
	}	
}

// Sends x_(n+1) into x_n. 
void Swap(double **a, double **b){
    double *tmp = *a;			// Recall (double ** a) is a pointer to an array of doubles.
    *a = *b;					// The content of 'a' (i.e. the array, i.e. the address of its first element) is assigned to the pointer tmp
    *b = tmp;					// The array to which 'b' points is assigned to the one to which 'a' points					
}

// Measures polarization
vector polarization(){
	double c=0,s=0;
	int i;
	vector pol;
	for (i=0;i<N;i++){
		c +=cos(*(op+i));
		s +=sin(*(op+i));
	}
	pol.x=c/N;
	pol.y=s/N;
	pol.delta = sqrt(c*c+s*s)/N;
	
	return pol;
}

// Statistics on spatial distribution
// Prints mean # of neighbours within a distance rc from each bird
double meanvic(cell *mosaic){  							// # of neighbours & eigevalues of adjacency matrix
	int i, j, k;
	double vic=0.0;
	//gsl_matrix_set_zero(adjacencym);
	
	for (i=0;i<N;i++){									// Cycle over birds
		for (j=0;j<24;j++){
			k=0;
			// cylcling over all the 5x5 crown
			while (k<mosaic[near[24*((pos+i)->tile) + j]].numData) {
				if(dist(i,mosaic[near[24*((pos+i)->tile) +j]].list[k])<=rc) {
					vic +=1;
					//gsl_matrix_set(adjacencym,i,mosaic[near[24*(int)((pos+i)->tile) +j]].list[k],1);}
					}	
				k++;
			}
		}
		
		/*k=0;
		while (k<mosaic[((pos+i)->tile)].numData) {
			vic +=1;
			//gsl_matrix_set(adjacencym,i,mosaic[((pos+i)->tile)].list[k],1);		//wrongly setting n_{ii}=1
			//adjacencym[i][mosaic[((pos+i)->tile)].list[k]]=1;  
			k++;
		}
		
		//gsl_matrix_set(adjacencym,i,i,0);				// resetting n_{ii}=0 which means no self interaction
		*/
		vic += mosaic[((pos+i)->tile)].numData;			// Replaces the last while cycle
	}
	
	return (vic-N)/N;									// Subtract myself first
}


// Evolution and real-time plotting of observables
void observable(double tfin, double factor[], cell *mosaic, cell *localbird){
	int i;
	int imax=(int)(tfin/factor[0]);
	vector pol;
	double mean_vic;
	
	/*
	// To print trajectory
	double T=factor[5]*factor[5]/(2*factor[0]*factor[8]*factor[8]);
	int k;
	char filename[50];
	sprintf(filename, "Trajectory_T%.3lf_v%.2lf_chi%.2lf_J8.txt", T, v0, 1/factor[8]);
	FILE *fs;
	fs=fopen (filename,"w");
	*/

	for (i=0;i<imax;i++){
		
		if (i%100 ==0){		// to save values at real time
		pol=polarization();
		//printf("%lf %lf %lf %lf \n", i*factor[0], pol.x,pol.y, pol.delta);
		printf("%lf %lf ", i*factor[0], pol.delta);
		
		// To print the trajectory
		//for (k=0;k<N;k++) fprintf(fs,"%lf %lf %lf %lf \n", (pos+k)->x, (pos+k)->y, cos(*(op+k)), sin(*(op+k)));

		// Clustering check
		//printf("%lf ", i*factor[0]);
		mean_vic = meanvic(mosaic);						// Print average # of neighbours
		printf("%lf\n", mean_vic);							
		}
		
		//evolve(factor, mosaic);
		evolve_field_topological(factor, mosaic, localbird, 0,0);
	}
	
	// If you were printing the trajectory
	// fclose(fs);
}

// Evolves and prints out the whole trajectory (positions and angles)
void trajectory(double tfin, double factor[], cell *mosaic){
	int i, k;
	int imax = (int)(tfin/factor[0]);
	
	double T=factor[5]*factor[5]/(2*factor[0]*factor[8]*factor[8]);
	char filename[50];
	sprintf(filename, "Windnum_T%.3lf_v%.2lf_chi%.2lf_J8.txt", T, v0, 1/factor[8]);
	FILE *fs;
	fs=fopen (filename,"w");
	
	for (i=0;i<imax;i++){
		if (i%100==0){
			for (k=0;k<N;k++){
				printf("%lf %lf %lf %lf \n", (pos+k)->x, (pos+k)->y, cos(*(op+k)), sin(*(op+k)));
				fprintf(fs, "%d %d %d \n",k, windx[k], windy[k]);
			}
		}
		evolve(factor, mosaic);
	}
}

// Computes the pol angle phi(t) when no field is applied
void wandering(double tfin, double factor[], cell *mosaic){	
	int i, k;
	int imax = (int)(tfin/factor[0]);
	double mean, sub;
	
	for (i=0;i<imax;i++){
		if (i%100==0){
			mean =0.0;
			for (k=0;k<N;k++) {
				mean += *(op + k);
			}
			
		mean /=N;
		if (i==0) sub =mean;
		printf("%lf %lf \n", i*factor[0], mean-sub);
		}
		
		evolve(factor, mosaic);
	}
}

// Computes euclidean distance between particle i and j
double dist(int i, int j){				
	double dx, dy, dist;
	dx = ((pos+i)->x )- ((pos+j)->x );
	dx -= L*round((double)(dx/L));
	dy = ((pos+i)->y )- ((pos+j)->y );
	dy -= L*round((double)(dy/L));
	
	dist = sqrt(dx*dx + dy*dy);	
	return dist;
}

// Simulates a turn induced by an external perturbing field
void turn_angle(double tfin, double factor[], cell *mosaic, double angle){
	int i, j, k, ifree, imax, kmax;
	int count;
	double window, free, shift, phi, phi_loc;
	double cosines, mean, time, sub=0.0, meanperp;
	vector pol;
	char filename[56];
	char filename2[56];
	char filename3[56];
	
	window = 5000;									//Duration of each turn (if window=tfin only one turn happens) (10000)
	free = 1000;									//Relaxation time after the turn (1000)
	kmax = (int)(tfin/window);						//# of turns
	imax = (int)((window - free)/ factor[0]);		//# of steps in each window
	ifree = (int)(free/factor[0]);					//# of relaxation steps
	pol = polarization();							//Initial polarization
	phi = atan2(pol.y, pol.x);
	//printf("Initial angle %.8lf \n", phi);
	
	sprintf(filename, "Pol_%d_T0.005_v%.2lf_chi%.2lf_J50.txt", A0, v0, 1/factor[8]);
	sprintf(filename2, "Trajectory_%d_T0.005_v%.2lf_chi%.2lf_J50.txt", A0, v0, 1/factor[8]);
	sprintf(filename3, "p130_%d_T0.005_v%.2lf_chi%.2lf_J50.txt", A0, v0, 1/factor[8]);
	
	FILE *ft;
	ft=fopen (filename,"w");
	
	// FILE *fd;
	// fd=fopen ("Turn_30_T0.005_v%.2lf_chi%.1lf_J50.txt","w");
	
	FILE *fd;
	fd=fopen (filename2,"w");
	
	FILE *fc;
	fc=fopen (filename3,"w");
	
	// Perturbation cycle
	for (k =0; k<kmax;k++){
		shift = phi + angle;
		shift = atan2(sin(shift), cos(shift));
		
		pol=polarization();
		printf("%lf \n", pol.delta);
		
		for (i=0;i<imax;i++){
			if (i%100==0){
				
				mean =0.0;
				meanperp =0.0;
				
				//printf("At %d angle for p is %lf \n", i, *(op+p)); 
				for (j=0;j<N;j++) { 	
					
					mean += *(op + j);
					meanperp += sin(*(op + j)-phi);
					
					//if (i<100000) {  				//set angle to save on disk
					fprintf(fc,"%lf ", *(op+j));
					fprintf(fd,"%lf %lf %lf %lf \n", (pos+j)->x, (pos+j)->y, cos(*(op+j)), sin(*(op+j)));	
				}
				// all velocities are rescaled
				time = (k *(imax + ifree)+ i)*factor[0];
				mean /= N;
				meanperp /=N;
				
				if (i==0) {
					sub=mean;
					}
					
				//if (sub>=1.5*M_PI && sub <=2*M_PI) sub =0.0;
				
				fprintf(ft, "%.8lf %.8lf %.8lf \n", time, mean-sub,meanperp);
				
				//listvic(p,mosaic);
				fprintf(fc, "\n");	
			}
			
			// Updates all lattice sites
			evolve_field(factor, mosaic, i*factor[0], shift);
		}
		
		fprintf(ft, "\n");
		fprintf(fc, "\n");

		// Relaxation cycle
		count = 0.0;
		phi = 0.0;
		cosines = 0.0;
		
		for (i=0;i<ifree; i++){
			if(i==0){
				initUP();
				startlattice();
			}
	
			if (i%100==0){
				pol = polarization ();
				phi_loc = atan2 (pol.y, pol.x);
				phi += phi_loc ;
				cosines += cos (phi_loc - shift);
				count += 1;
				for (j=0;j<N;j++) {
					fprintf(fd,"%lf %lf %lf %lf \n", (pos+j)->x, (pos+j)->y, 1*cos(*(op+j)), 1*sin(*(op+j)));
				}
				//listvic(p,mosaic);
			}
			
		// Updates all lattice sites
		evolve(factor, mosaic);
		}
		
		/*for (i=0;i<(int)(ifree/factor[0]); i++){

			if (i%100==0){
				pol = polarization ();
				phi_loc = atan2 (pol.y, pol.x);
				phi += phi_loc ;
				cosines += cos (phi_loc - shift);
				count += 1;
				for (j=0;j<N;j++) {
					fprintf(fd,"%lf %lf %lf %lf \n", (pos+j)->x, (pos+j)->y, 1*cos(*(op+j)), 1*sin(*(op+j)));}
				listvic(p,mosaic);
			}
			
		// Updates all lattice sites
		evolve(factor, mosaic);
		}*/
		
		phi /= count;
		cosines /= count;
		printf ("Phi %.8lf \t and %.8lf \n" ,phi ,1-cosines);	
	}
	
	fclose(ft);
	fclose(fd);
	fclose(fc);
}

// Finds bird in the middle of the lattice
void find_p(cell * mosaic, double v){
	double y = L/2;							// Center of the mosaic, it will move following the (very polarized) flock
	double x= L/2;							// Center of the mosaic
	//double x=L;							// Side of the mosaic
	
	y += v*T_term;							// After thermalization time
	y = (int) reshift(y);
	p = mosaic[indice(x,y)].list[0];		// Update index of the central bird
}

// Random number generation
void boxmuller(double *address){
	double u, theta;
	u=drand48();
	theta=drand48();
	
	*(address)=sqrt(-2*log(u))*cos(2*M_PI*theta);
	*(address +1)=sqrt(-2*log(u))*sin(2*M_PI*theta);	
}

// Index to vector converters
int x(int k){
	return k%L;}
	
int y(int k){
	return (int)(k/L);}
	
int indice(int x, int y){
	x= (x+(1+(int)(abs(x)/L))*L)%L;					//This line maps the real position into (0,L)
	y= (y+(1+(int)(abs(y)/L))*L)%L;				
	return x+L*y;}
	
double reshift(double posit){						//Bring the position into (0,L)
	int k;
	k=(int)(fabs(posit)/L) + 1;
	posit = fmod(posit + k*L,L); 
	return posit;}
	
void help(){
	printf ("Invalid input of parameters. Insert values separated by a single space: \n\
	temperature \t flock speed \t inertia  \n.");
	exit(EXIT_FAILURE);
}

// Use it to print the mosaic or the sorted lists
void printCells(cell * mosaic){
	int i,j;
	int sum=0;
	for(i=0; i<N; i++){
		printf("# On tile/bird %d : \t",i);
		for(j=0; j<mosaic[i].numData; j++) printf("%d\t", mosaic[i].list[j]);
		printf("\n# Total birds %d\n# Distances from %d-th bird: \t", mosaic[i].numData, i);
		for(j=0; j<mosaic[i].numData; j++) printf("%f\t", dist(mosaic[i].list[j],i));
		printf("\n");
		sum += mosaic[i].numData;
	}
	printf("Total birds: %d\n", sum);
}
