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


	//EXERCISE 3.1.1
	
	unsigned int M=1E5, N=100;

	vector<vector<double>> vec;

	Statistics s(M,N);

	s.Set_Parameters(100.,'S');
	s.Set_Parameters(1.,'T');
	s.Set_Parameters(100.,'K');
	s.Set_Parameters(0.1,'r');
	s.Set_Parameters(0.25,'s');

	vec=s.European_Option(rnd,vec,1,'c');		//The char select the option type: 'c'-> call option , 'p'-> put option 
							//Further details are available in the attached code files "class.cpp" and "class.h"
	vec=s.European_Option(rnd,vec,1,'p');		
			
	s.Print(vec, "Datas/data1.dat");			//Saving datas on file "data1.dat"

	
	//EXERCISE 3.1.2

	vec.clear();

	vec=s.European_Option(rnd,vec,100,'c');	
	vec=s.European_Option(rnd,vec,100,'p');

	s.Print(vec, "Datas/data2.dat");			//Saving datas on file "data1.dat"


	return 0;

}
