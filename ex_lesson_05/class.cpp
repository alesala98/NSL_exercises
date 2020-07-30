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
 
vector<vector<double>> Statistics :: Metropolis(Random rnd, vector<vector<double>> support1, char a, char b, double delta){

	double sum;
	
	vector<vector<double>> support;
	vector<double> ave(m_N,0.), ave2(m_N,0.), sum_prog(m_N,0.), sum2_prog(m_N,0.), err_prog(m_N,0.), coord(3,0.), newcoord(3,0.);   
 
	if(support1.size()!=0)						//Check if the vector support1 given as input has a finite dimension to preserve its information
		for(unsigned int i=0; i<support1.size(); i++)
			support.push_back(support1[i]);
	
	if(a=='s'){

		vector<double> x(m_N,0.);
		iota(x.begin(),x.end(),0);
		support.push_back(x);

	}

	for(unsigned int i=0; i<m_N; i++){

		sum=0;
		
		fill(coord.begin(), coord.end(), 0.);			//Starting point is (0,0,0), but you could change it to any value you like. However, you should be 
									//cautious, since with the same number of simulations and blocks Metropolis algorithm gives back wrong
		for(unsigned int j=0; j<m_L; j++){			//extimations. For example: starting point (20,20,20) -> r_unif_extimation = 1.9

			
			for(unsigned int k=0; k<3; k++){

				if(a=='s'){
					if(b=='u')
						newcoord[k]=coord[k]+(rnd.Rannyu()-0.5)*delta;
					if(b=='g')
						newcoord[k]=rnd.Gauss(coord[k],delta);
				}

				if(a=='p'){
					if(b=='u')
						newcoord[k]=coord[k]+(rnd.Rannyu()-0.5)*delta;
					if(b=='g')
						newcoord[k]=rnd.Gauss(coord[k],delta);
				}
					
			}

			double A=min(1., ProbWaveFunc(newcoord,a)/ProbWaveFunc(coord,a));

			if(A>=rnd.Rannyu())
				for(unsigned int k=0; k<3; k++)
					coord[k]=newcoord[k];
	
			sum+=sqrt(pow(coord[0],2)+pow(coord[1],2)+pow(coord[2],2));

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

double Statistics :: ProbWaveFunc(vector<double> coord, char a){

	double r;

	r=sqrt(pow(coord[0],2)+pow(coord[1],2)+pow(coord[2],2));

	if(a=='s')
		return (1./M_PI)*exp(-2*r);					//S wave function
	
	else
		return (1./(32.*M_PI))*pow(r*exp(-r/2.)*(coord[2]/r),2);	//P wave function

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


