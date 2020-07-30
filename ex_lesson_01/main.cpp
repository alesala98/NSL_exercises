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


	//EXERCISE 1.1

	//1.Integral mean value & 2.Integral standard deviation

	unsigned int M=1E5, N=100;

	vector<vector<double>> vec;

	Statistics s(M,N);

	vec=s.Ave_StDev(rnd,vec,'m');		//'m' -> mean, 'v' -> standard deviation
						//The code uses the same class to calculate the mean value and the standard deviation; the character "a" is used to distinguish different 
	vec=s.Ave_StDev(rnd,vec,'v');		//calculations performed by the method Ave_StDev. Further details are available in the attached code files "class.cpp" and "class.h"
	
	//3.Pearson's cumulative test

	s.SetParameters(1E4,'m');
	s.SetParameters(100,'n');

	vec=s.Pearson(rnd,vec);

	s.Print(vec, "data1.dat");		//Saving datas on file "data1.dat"

	
	//EXERCISE 1.2

	vec.clear();

	vec=s.Central_Lim(rnd,vec,1E4,'u');
	vec=s.Central_Lim(rnd,vec,1E4,'e');
	vec=s.Central_Lim(rnd,vec,1E4,'l');
	
	s.Print(vec, "data2.dat");		//Saving datas on file "data2.dat"
		

	//EXERCISE 1.3

	vec.clear();

	s.SetParameters(1E5,'m');
	s.SetParameters(100,'n');

	vec=s.Buffon(rnd,vec,30,10);		//s.Buffon(Random,vector<vector<double>>,d,L)	d -> distance between lines, L-> needle length
						
	s.Print(vec, "data3.dat");		//Saving datas on file data3.dat


	return 0;

}
