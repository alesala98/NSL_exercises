//The comments in this file are not detailed. See qmc1d.cpp for a better explanation.

void readInput();   	//Reads input from the file "input.dat"
void deleteMemory(); 	//Handles the dynamic allocation of memory
void initialize();  	//Initializes the variables
void consoleOutput(); 	//Writes the output on the screen
                                                                                                                                                                                                                                  
double potential_density_matrix(double val, double val_next);
double u_prime(double x, int m);
double u_sec(double x, int m);


//Potential_density_matrix returns only the potential part of the correlation between two adjacent timeslices.

double external_potential(double);  		//This is the external potential definition
double external_potential_prime(double); 	//...and here goes its first derivative
double external_potential_second(double); 	//... and its second derivative 


//The derivatives are necessary for the evaluation of the kinetic estimator, because it contains the laplacian operator! 
                                                                                                                                                                                                                     
void translation(); 		//Performs a rigid translation
void brownianBridge();  	//Reconstructs a segment of the polymer with a free particle propagation. 
void brownianMotion(int);  	//Reconstructs a segment at the extremities of the polymer with a free particle propagation. 
                                                                                                                 
double variationalWaveFunction(double);  		//VariationalWaveFunction is the variational wave function that is projected in a PIGS simulation.
double variationalWaveFunction_second(double);		//As for the potential, you have to specify its first and second derivative for the evaluation of the kinetic local energy.
double variationalLocalEnergy(double val);
                                                                                                               

int index_mask(int); 	//Index_mask is just a compatibility function that takes into account whether the polymer is open (PIGS) or closed in periodic boundary contitions (PIMC-ring polymer).


void upgradeAverages(); 	//At every MCSTEP accumulates the estimators values.

void upgradeHistogram(); 	//Fills the histogram of positions foreach MCSTEP
void endBlock(); 		//Finalizes the averages at the end of each block

double kineticEstimator(double,double);  	//Evaluates the kinetic energy along the polymer
void finalizePotentialEstimator();		//The last three functions are called at the end of the simulation, basically they average over each
void finalizeKineticEstimator();		//block and evaluate the error on the block average. This is an application of the central limit
void finalizeHistogram();			//theorem.
