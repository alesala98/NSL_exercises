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


	//EXERCISE 8.2
	
	unsigned int M=1E6, N=100;
	bool opt=false;
	vector<vector<double>> vec;

	Statistics s(M,N,rnd);

	s.Optimization(opt);			//Once the optimal parameters are found, you should change opt to false

	vec=s.Metropolis(vec); 			//Compute the hamiltonian of the system with its uncertainty and the probability wave function

	s.Print(vec,"../Datas/hamiltonian.dat");

	return 0;

}
