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


//Function generating random angles between 0 and pi using accept/reject method

double gen_pi(Random rnd){

	double x=rnd.Rannyu(-1.,1.), y=rnd.Rannyu(0.,1.);

	while((pow(x,2)+pow(y,2))>=1){
		x=rnd.Rannyu(-1.,1.);
		y=rnd.Rannyu(0.,1.);
	}
	
	return acos(x/sqrt(pow(x,2)+pow(y,2)));	
	
}


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

	//1.Integral mean value & 2.Integral standard deviation

	iota(x.begin(),x.end(),0);

	for(unsigned int i=0; i<m_M; i++){

		if(i<m_N){
			r.push_back(rnd.Rannyu());
			x[i]*=m_L;
		} 
		else 
			r.push_back(rnd.Rannyu());

	}

	if(a=='m')
		support.push_back(x);

	for(unsigned int i=0; i<m_N; i++){

		sum=0;
		
		for(unsigned int j=0; j<m_L; j++){

			k=j+i*m_L;

			if(a=='m')
				sum+=r[k];

			if(a=='v')
				sum+=pow(r[k]-0.5,2);

		}

		ave[i]=sum/m_L;			//Montecarlo integration
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

		if(a=='m')			//Statistical uncertainty: mean
			sum_prog[i]-=0.5; 
			
		
		if(a=='v')			//Statistical uncertainty: standard deviation
			sum_prog[i]=(12.*sum_prog[i]-1)/12.;		

	}

	support.push_back(sum_prog);
	support.push_back(err_prog);

	return support;

}

vector<vector<double>> Statistics :: Pearson(Random rnd, vector<vector<double>> support1){		//3.Pearson's cumulative test

	vector<vector<double>> support;

	if(support1.size()!=0)						//Check if the vector support1 given as input has a finite dimension to preserve its information
		for(unsigned int i=0; i<support1.size(); i++)
			support.push_back(support1[i]);

	vector<double> chij(m_N,0.), accu(m_N,0.), rand;

	for(unsigned int i=0; i<m_N; i++){				

		rand.clear();
		fill(accu.begin(),accu.end(),0.);
		
		for(unsigned int k=0; k<m_M; k++)
			rand.push_back(rnd.Rannyu());

		sort(rand.begin(),rand.end());

		unsigned int k=0;
		for(unsigned int j=0; j<m_N; j++)
			while(rand[k]<((j+1.)/m_N) and k<m_M){
				accu[j]++;
				k++;
			}
		
		for(unsigned int l=0; l<m_N; l++)
			chij[i]+=pow((accu[l])-(m_M/m_N),2)/(m_M/m_N);
			
	}

	support.push_back(chij);

	return support;

}

vector<vector<double>> Statistics :: Central_Lim(Random rnd, vector<vector<double>> support1, unsigned int n, char x){

	vector<double> sum(n,0.), sum2(n,0.), sum10(n,0.), sum100(n,0.);		
																			
	vector<vector<double>> support;						

	if(support1.size()!=0)							//Check if the vector support1 given as input has a finite dimension to preserve its information
		for(unsigned int i=0; i<support1.size(); i++)
			support.push_back(support1[i]);									

	switch(x){
		case 'u':							//Uniform
			for(unsigned int i=0; i<n; i++){

				sum[i]=rnd.Rannyu(1.,6.);			//No sum

				for(unsigned int j=0; j<2; j++)			//Sum over two random numbers
					sum2[i]+=rnd.Rannyu(1.,6.)/2.;

				for(unsigned int j=0; j<10; j++)		//Sum over 10 random numbers	
					sum10[i]+=rnd.Rannyu(1,6)/10.;

				for(unsigned int j=0; j<100; j++)		//Sum over 100 random numbers	
					sum100[i]+=rnd.Rannyu(1,6)/100.;

			}
		
			support.push_back(sum);
			support.push_back(sum2);
			support.push_back(sum10);
			support.push_back(sum100);

			return support;

		case 'e':							//Exponential
			for(unsigned int i=0; i<n; i++){

				sum[i]=rnd.Exp(1.);				//No sum

				for(unsigned int j=0; j<2; j++)			//Sum over two random numbers
					sum2[i]+=rnd.Exp(1.)/2.;

				for(unsigned int j=0; j<10; j++)		//Sum over 10 random numbers	
					sum10[i]+=rnd.Exp(1.)/10.;

				for(unsigned int j=0; j<100; j++)		//Sum over 100 random numbers	
					sum100[i]+=rnd.Exp(1.)/100.;

			}
		
			support.push_back(sum);
			support.push_back(sum2);
			support.push_back(sum10);
			support.push_back(sum100);
	
			return support;

		case 'l':								//Cauchy-Lorentz
			for(unsigned int i=0; i<n; i++){

				sum[i]=rnd.Lorentz(1.,0.);				//No sum

				for(unsigned int j=0; j<2; j++)				//Sum over two random numbers
					sum2[i]+=rnd.Lorentz(1.,0.)/2.;

				for(unsigned int j=0; j<10; j++)			//Sum over 10 random numbers	
					sum10[i]+=rnd.Lorentz(1.,0.)/10.;

				for(unsigned int j=0; j<100; j++)			//Sum over 100 random numbers	
					sum100[i]+=rnd.Lorentz(1.,0.)/100.;;

			}
		
			support.push_back(sum);
			support.push_back(sum2);
			support.push_back(sum10);
			support.push_back(sum100);

			return support;

		default:
			cerr << endl << "Error: no such parameter exists" << endl;
			break;
	}

	return support;

}

vector<vector<double>> Statistics :: Buffon(Random rnd, vector<vector<double>> support, double d, double L){				//d -> distance between lines
																	//L -> needle length
	unsigned int n_hit;

	vector<double> y, theta, pi(m_N,0), pi_err(m_N,0), pi_prog(m_N,0.), pi_prog2(m_N,0.), x(m_N,0.);

	iota(x.begin(),x.end(),0);

	for(unsigned int i=0; i<m_M; i++){
		y.push_back(rnd.Rannyu(0,d));
		theta.push_back(gen_pi(rnd));
	}
	
	for(unsigned int i=0; i<m_N; i++){

		n_hit=0;
		x[i]*=m_L;

		for(unsigned int j=0; j<m_L; j++){

			unsigned int k=j+i*m_L;

			if(y[k]-(L*sin(theta[k])/2.)<=0 or y[k]+(L*sin(theta[k])/2.)>=d)		 
				n_hit++;									

		}

		pi[i]=2*m_L*L/(n_hit*d);

	}

	for(unsigned int i=0; i<m_N; i++){

		for(unsigned int j=0; j<i+1; j++){

			pi_prog[i]+=pi[j];
			pi_prog2[i]+=pow(pi[j],2);

		}

		pi_prog[i]/=(i+1);
		pi_prog2[i]/=(i+1);
		pi_err[i]=Error(pi_prog2,pi_prog,i);
	
	}
		
	support.push_back(x);
	support.push_back(pi_prog);
	support.push_back(pi_err);
		
	return support;

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

double Statistics :: Error(vector<double> sum2, vector<double> sum_prog, double n){

		if(n==0)	
    			return 0;
		else
			return sqrt(abs((sum2[n]-pow(sum_prog[n],2)))/n);		
		
}
