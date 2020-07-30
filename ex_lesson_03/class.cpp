#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <vector>

#include "random.h"
#include "class.h"

using namespace std;


//Class Methods

Statistics :: Statistics(unsigned int M, unsigned int N){

	m_M=M;
	m_N=N;
	m_L=M/N;

}	
	
Statistics :: ~Statistics(){};     	
						
void Statistics :: Set_Iterations(unsigned int n, char c){	

	switch(c){
		case 'm':
			m_M=n;
			m_L=m_M/m_N;
			break;
		case 'n':
			m_N=n;
			m_L=m_M/m_N;
			break;
		default:
			cerr << endl << "Error: no such parameter exists" << endl;
			break;
	}

}

void Statistics :: Set_Parameters(double n, char c){

	switch(c){
		case 'S':
			m_S0=n;
			break;
		case 'T':
			m_T=n;
			break;
		case 'K':
			m_K=n;
			break;
		case 'r':
			m_r=n;
		case 's':
			m_sigma=n;
			break;
		default:
			cerr << endl << "Error: no such parameter exists" << endl;
			break;
	}

}
 
vector<vector<double>> Statistics :: European_Option(Random rnd, vector<vector<double>> support1, unsigned int nsteps, char a){

	unsigned int k;
	double sum, s;
	
	vector<vector<double>> support;
	vector<double> N, x(m_N,0.), ave(m_N,0.), ave2(m_N,0.), sum_prog(m_N,0.), sum2_prog(m_N,0.), err_prog(m_N,0.);   
 
	if(support1.size()!=0)						//Check if the vector support1 given as input has a finite dimension to preserve its information
		for(unsigned int i=0; i<support1.size(); i++)
			support.push_back(support1[i]);

	iota(x.begin(),x.end(),0);

	for(unsigned int i=0; i<nsteps*m_M; i++){

		if(i<m_N){

			x[i]*=m_L;
			N.push_back(rnd.Gauss(0,1));

		} 
		else
			N.push_back(rnd.Gauss(0,1));

	}

	if(a=='c')
		support.push_back(x);

	for(unsigned int i=0; i<m_N; i++){

		sum=0;
		
		for(unsigned int j=0; j<m_L; j++){

			s=m_S0;

			for(unsigned int l=0; l<nsteps; l++){				//Discretization over time. If nsteps==1, then the simulation is carried on without  
											//any kind of discretization, i.e. the assignment of exercise 3.1.1
				k=j+l+i*m_L;
				s*=exp((m_r-0.5*pow(m_sigma,2))*(m_T/nsteps)+(m_sigma*N[k]*sqrt(m_T/nsteps)));

			}

			if(a=='c')							//Call option
				sum+=exp(-m_r*m_T)*max(0.,s-m_K);

			if(a=='p')							//Put option
				sum+=exp(-m_r*m_T)*max(0.,m_K-s);

		}

		ave[i]=sum/m_L;			
		ave2[i]=pow(ave[i],2);

	}

	for(unsigned int i=0; i<m_N; i++){

		for(unsigned int j=0; j<i+1; j++){
			sum_prog[i]+=ave[j];
			sum2_prog[i]+=ave2[j];
		}

		sum_prog[i]/=(i+1);             //Cumulative average				
    		sum2_prog[i]/=(i+1);		//Cumulative square average
				
		err_prog[i]=Error(sum2_prog,sum_prog,i);
				
	}

	support.push_back(sum_prog);
	support.push_back(err_prog);

	return support;

}

double Statistics :: Error(vector<double> sum2, vector<double> sum_prog, double n){

		if(n==0)	
    			return 0;
		else
			return sqrt(abs((sum2[n]-pow(sum_prog[n],2)))/n);		
		
}
	
void Statistics :: Print(vector<vector<double>> vec, string filename){

	ofstream WriteData;
   	WriteData.open(filename);

	if(WriteData.is_open()){
		
		for(unsigned int j=0; j<vec[0].size(); j++){
		
			for(unsigned int i=0; i<vec.size(); i++)	
      				WriteData << vec[i][j] << " " ;
			
			WriteData << endl ;

		}
			
   	} else cerr << "Unable to open " << filename << endl;

  	WriteData.close();

}


