#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <vector>

#include "random.h"
#include "class.h"

using namespace std;


//Class Methods

Statistics :: Statistics(unsigned int M, unsigned int N, Random r){

	m_M=M;
	m_N=N;
	m_L=M/N;
	
	m_rnd=r;

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

void Statistics :: Optimization(bool opt){

	if(opt){

		cout << "Computing optimal parameters... " << endl << endl;

		double delta=4.;				//Metropolis step to obtain 50% acceptance rate
		double ene_true=10, mu_true, sigma_true;

        	for(double mu=0.7; mu<=0.9; mu+=0.005){

           		for(double sigma=0.5; sigma<=0.7; sigma+=0.005){

               			double sum=0, x=0., xnew=0.;
                		unsigned int attempted=0, accepted=0;

                		for(unsigned int i=0; i<m_M; i++){

                    			xnew=x+(m_rnd.Rannyu()-0.5)*delta;

                    			double A=min(1., pow(ProbWaveFunc(xnew,mu,sigma,'p'),2.)/pow(ProbWaveFunc(x,mu,sigma,'p'),2.));
                    		
                   			if(A >= m_rnd.Rannyu()){
                        			x=xnew;
                        			accepted++;
                   			}

                   			attempted++;
                    			sum+=(-ProbWaveFunc(x,mu,sigma,'h')*0.5+V(x)*ProbWaveFunc(x,mu,sigma,'p'))/(ProbWaveFunc(x,mu,sigma,'p'));

                		}

				cout << "Acceptance rate: " << (double)(accepted)/(double)(attempted) << "  "<< endl;

                		sum/=m_M;
                		if(sum <= ene_true){
					cout << "New minimal energy extimation: " << sum << endl;
                  			ene_true=sum;
                    			mu_true=mu; 
					sigma_true=sigma;
                		}
	
            		}

        	}

        	cout << endl << "Optimized parameters " << endl << "Mu: " << setprecision(6) << mu_true << endl << "Sigma: " << setprecision(6) << sigma_true << endl;

	}

}
 
vector<vector<double>> Statistics :: Metropolis(vector<vector<double>> support){

	double sum, x, xnew, delta=4.;
	double mu=0.81, sigma=0.62;		//Obtained from previous optimization

	vector<double> block(m_N,0.), ave(m_N,0.), ave2(m_N,0.), sum_prog(m_N,0.), sum2_prog(m_N,0.), err_prog(m_N,0.); 
	vector<vector<int>> histogram(m_N, vector<int>(m_N,0.));  

	iota(block.begin(),block.end(),0);
	support.push_back(block);

	cout << endl << endl << "Computing energy..." << endl;

	for(unsigned int i=0; i<m_N; i++){

		sum=0;
		x=0;
			
		for(unsigned int j=0; j<m_L; j++){

			xnew=x+(m_rnd.Rannyu()-0.5)*delta;

                    	double A=min(1., pow(ProbWaveFunc(xnew,mu,sigma,'p'),2.)/pow(ProbWaveFunc(x,mu,sigma,'p'),2.));
                    		
                   	if(A >= m_rnd.Rannyu())
                        	x=xnew;

                    	sum+=(-ProbWaveFunc(x,mu,sigma,'h')*0.5+V(x)*ProbWaveFunc(x,mu,sigma,'p'))/(ProbWaveFunc(x,mu,sigma,'p'));
	
			int bin=((x+3.)/0.06);		
			histogram[i][bin]++;		//Filling histogram

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

	cout << endl << endl << "Computing probability wave function..." << endl << endl;

	ofstream WriteOutput;
	WriteOutput.open("../Datas/histogram.dat");

	for(unsigned int ibin=0; ibin<100; ibin++){
        	
           	fill(sum_prog.begin(), sum_prog.end(), 0.);
           	fill(sum2_prog.begin(), sum2_prog.end(), 0.);
            	fill(err_prog.begin(), err_prog.end(), 0.);
        
        	for(unsigned int i=0; i<m_N; i++){

           		for(unsigned int j=0; j<i+1; j++){
                		sum_prog[i]+=histogram[j][ibin];
                		sum2_prog[i]+=pow(histogram[j][ibin],2.);
            		}
            
            		sum_prog[i]/=(i+1);
            		sum2_prog[i]/=(i+1);
            		err_prog[i]=Error(sum2_prog,sum_prog,i);

            		if(i==m_N-1)
                		WriteOutput << (-3+0.03)+(ibin*0.06) << " " << sum_prog[i]/(m_L*0.06) << " " << err_prog[i]/(m_L*0.06) << endl;

            	}
        }
	
	WriteOutput.close();

	return support;

}

double Statistics :: ProbWaveFunc(double x, double mu, double sigma, char a){

	if(a=='p')
		return exp(-pow(x-mu,2.)/(2.*sigma*sigma)) + exp(-pow(x+mu,2.)/(2.*sigma*sigma));					//p --> Psi(x)
	
	else
		return exp(-pow(x-mu,2.)/(2.*sigma*sigma))*(pow((x-mu)/(sigma*sigma),2.)-(1./(sigma*sigma))) +				//h --> H*Psi(x)
		       exp(-pow(x+mu,2.)/(2.*sigma*sigma))*(pow((x+mu)/(sigma*sigma),2.)-(1./(sigma*sigma)));

}

double Statistics :: V(double x){

	return pow(x,4)-2.5*pow(x,2);

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


