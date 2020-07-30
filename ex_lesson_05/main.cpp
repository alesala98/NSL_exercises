#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

#include "random.h"
#include "class.h"

using namespace std;

int main (int argc, char *argv[]){


	//Initialization Random Generator

   	Random rnd;
	rnd.SaveSeed();


	//EXERCISE 5.1
	
	unsigned int M=1E6, N=100;

	vector<vector<double>> vec;

	Statistics s(M,N);

	vec=s.Metropolis(rnd,vec,'s','u',2);			//First char select the wave function type: 'c'-> s orbital , 'p'-> p orbital
	vec=s.Metropolis(rnd,vec,'p','u',6);			//Second char select the transition probability distribution: 'u' -> uniform , 'g' -> gaussian
								//Double parameter sets delta-step used in Metropolis algorithm to obtain an acceptance near 50% 
								//Further details are available in the attached code files "class.cpp" and "class.h"

	s.Print(vec, "Datas/data_unif.dat");			//Saving datas on file "data_unif.dat" (uniform transition probability)
	
	vec.clear();

	vec=s.Metropolis(rnd,vec,'s','g',0.75);
	vec=s.Metropolis(rnd,vec,'p','g',1.9);

	s.Print(vec, "Datas/data_gauss.dat");			//Saving datas on file "data_gauss.dat" (gaussian transition probability)
	
	return 0;

}
