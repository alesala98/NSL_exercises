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


	//EXERCISE 2.1
	
	//1.Integral using uniform distribution & 2.Integral using importance sampling method

	unsigned int M=1E6, N=100;

	vector<vector<double>> vec;

	Statistics s(M,N);

	vec=s.Ave_StDev(rnd,vec,'s');		//'s' -> standard integration, 'i' -> importance sampling integration
						//The code uses the same class to calculate the mean value and the standard deviation; the character "a" is used to distinguish different 
	vec=s.Ave_StDev(rnd,vec,'i');		//calculations performed by the method Ave_StDev. Further details are available in the attached code files "class.cpp" and "class.h"
	
	s.Print(vec, "data1.dat");		//Saving datas on file "data1.dat"


	//EXERCISE 2.2

	//1.Discrete random walk & 2.Continuum random walk

	vec.clear();

	s.SetParameters(10E5,'m');
	s.SetParameters(100,'n');

	vec=s.RandomWalk(rnd,vec,100,'d');	//'d' -> discrete case, 'i' -> continuum case
						//The code uses the same class to simulate both required random walks; the character "a" is used to distinguish different 
	vec=s.RandomWalk(rnd,vec,100,'c');	//calculations performed by the method RandomWalk. Further details are available in the attached code files "class.cpp" and "class.h"

	s.Print(vec, "data2.dat");		//Saving datas on file "data2.dat"
	
	return 0;

}
