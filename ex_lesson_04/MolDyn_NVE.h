//Parameters, observables
const int m_props=500;
unsigned int n_props, igofr;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp, bin_size, nbins, blk_norm;
double blk_av[m_props], glob_av[m_props], glob_av2[m_props];

//Averages
double acc,att;

//Configuration
const unsigned int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

//Thermodynamical state
unsigned int npart;
double energy,temp,vol,rho,box,rcut;

//Simulation
unsigned int nstep; 
unsigned int nblock;
int seed;
double delta;
bool restart;

//Functions
void Input(void);
void Reset(unsigned int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
void Averages(unsigned int);
void Equilibration();
double Force(int, int);
double Pbc(double);
double Error(double,double,unsigned int);
