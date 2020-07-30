#ifndef __NVT__
#define __NVT__

#include "random.h"

//Random numbers
int seed[4];
Random rnd;

//Parameters, observables
const int m_props=1000;
int n_props, iv, iw, igofr;
double vtail,ptail,bin_size,nbins,sd;
double walker[m_props];

//Averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_pot,stima_pres,err_pot,err_press,err_gdir;

//Configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part];

//Thermodynamical state
int npart;
double beta,temp,vol,rho,box,rcut;

//Simulation
int nstep, nblk;
double delta;

//Pi
const double pi=M_PI;

//Functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);

#endif
