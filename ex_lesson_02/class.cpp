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
						
void Statistics :: SetParameters(unsigned int n, char c){	

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
 
vector<vector<double>> Statistics :: Ave_StDev(Random rnd, vector<vector<double>> support1, char a){

	unsigned int k;
	double sum;
	
	vector<vector<double>> support;
	vector<double> r, x(m_N,0.), ave(m_N,0.), ave2(m_N,0.), sum_prog(m_N,0.), sum2_prog(m_N,0.), err_prog(m_N,0.);   
 
	if(support1.size()!=0)						//Check if the vector support1 given as input has a finite dimension to preserve its information
		for(unsigned int i=0; i<support1.size(); i++)
			support.push_back(support1[i]);

	iota(x.begin(),x.end(),0);

	for(unsigned int i=0; i<m_M; i++){

		if(i<m_N){

			x[i]*=m_L;

			if(a=='s')				//Uniform
				r.push_back(rnd.Rannyu());

			if(a=='i')				//Importance sampling
				r.push_back(rnd.Sampling());


		} 
		else{

			if(a=='s')				//Uniform
				r.push_back(rnd.Rannyu());

			if(a=='i')				//Importance sampling
				r.push_back(rnd.Sampling());

		}

	}

	if(a=='s')
		support.push_back(x);

	for(unsigned int i=0; i<m_N; i++){

		sum=0;
		
		for(unsigned int j=0; j<m_L; j++){

			k=j+i*m_L;

			if(a=='s')
				sum+=M_PI*cos(r[k]*M_PI/2.)/2.;

			if(a=='i')
				sum+=(M_PI*cos(r[k]*M_PI/2.)/2.)/(2*(1-r[k]));

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

		sum_prog[i]-=1.; 
				
	}

	support.push_back(sum_prog);
	support.push_back(err_prog);

	return support;

}

vector<vector<double>> Statistics :: RandomWalk(Random rnd, vector<vector<double>> support1, unsigned int nsteps, char a){

	vector<vector<double>> support;
	vector<double> ave(nsteps+1,0.), ave2(nsteps+1,0.), r(nsteps+1,0), x(nsteps+1,0), accu(3,0), err; 

	if(support1.size()!=0)							//Check if the vector support1 given as input has a finite dimension to preserve its information
		for(unsigned int i=0; i<support1.size(); i++)
			support.push_back(support1[i]);
			
		
	iota(x.begin(),x.end(),0);

	if(a=='d')
		support.push_back(x);

	for(unsigned int i=0; i<m_N; i++){					//Random walk iterations (m_L*m_N)

		fill(r.begin(), r.end(), 0);

		for(unsigned int j=0; j<m_L; j++){				//Random walk iterations in each block

			fill(accu.begin(), accu.end(), 0);

			if(a=='d'){

				for(unsigned int k=0; k<nsteps+1; k++){		//Random steps in each block

					unsigned int ran_dir=rnd.Rannyu(0,3);
					
					if(rnd.Rannyu()<0.5)
						accu[ran_dir]++;
					else
						accu[ran_dir]--;
				
					for(unsigned int l=0; l<3; l++)
						r[k]+=pow(accu[l],2);

				}

			}

			if(a=='c'){

				for(unsigned int k=0; k<nsteps+1; k++){		//Random steps in each block

					double theta=2*M_PI*rnd.Rannyu(), phi=acos(1-2*rnd.Rannyu());			//Sampling theta and phi using this method avoids generating random numbers
					accu[0]+=sin(phi)*cos(theta);							//exibiting a tendency to be dense near the poles and scarce near the equator
        				accu[1]+=sin(phi)*sin(theta);
       					accu[2]+=cos(phi);
						
					for(unsigned int l=0; l<3; l++)
						r[k]+=pow(accu[l],2);

				}

			}
			
		}

		for(unsigned int j=0; j<nsteps+1; j++){
		
			ave[j]+=(r[j]/m_L);
			ave2[j]=pow(ave[j],2);

		}

	}

	for(unsigned int j=0; j<nsteps+1; j++){

		ave[j]/=m_N;
		ave2[j]/=m_N;
											//The error in each step is computed through partial derivatives:
		err.push_back(Err_RW(ave2,ave,j));					//err_step(x)=(1/2)*x^(-1/2)*err(x)

		ave[j]=sqrt(ave[j]);							//Requested result sqrt((r_n)^2)
	
	}

	support.push_back(ave);
	support.push_back(err);

	return support;
	
}

double Statistics :: Error(vector<double> sum2, vector<double> sum_prog, double n){

		if(n==0)	
    			return 0;
		else
			return sqrt(abs((sum2[n]-pow(sum_prog[n],2)))/n);		
		
}

double Statistics :: Err_RW(vector<double> sum2, vector<double> sum, double n){

		if(n==0)	
    			return 0;
		else
			return (0.5/sqrt(sum[n]))*sqrt(abs((sum2[n]-pow(sum[n],2)))/m_N);
		
		
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


