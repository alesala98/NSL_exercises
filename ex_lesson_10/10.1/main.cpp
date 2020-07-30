#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

#include "random.h"
#include "class.h"

using namespace std;

int main (int argc, char *argv[]){

	// Initialization Random Generator
   	Random rnd;
	rnd.SaveSeed();

	// EXERCISE 10.1
	unsigned int ncit = 32, t_gen = 2000, t_ind = 300;	// Number of cities, temperature generations and steps per generation
	double t_start = 4., t_end = 0.001;			// Starting - ending temperatures
	
	cout << endl << "Traveling salesman problem" << endl << endl;
	cout << "Start temperature:     " << t_start << endl;
	cout << "End temperature:       " << t_end << endl;
	cout << "Number of cities:      " << ncit << endl;
	
	// Initialization: creating world and allocating needed classes/variables
	cout << endl << "Creating world..." << endl;

	Generation g(ncit,t_gen,rnd);

	vector<City> cities = Generate_Cities(32,'s',rnd);		// Choose 'c' if you want cities on a circumference or 's' if you want them inside a 
									// square. Further details can be found in the attached 'class.cpp' and 'class.h' files
	vector<int> path(ncit);
	path = g.Gen_Individual();
	g.Check(path);

	cout << endl << "World created, starting simulation (annealing method) ... " << endl;

	double temp = t_start;

	// Simulation
	for(unsigned int i=0; i<t_gen; i++){
		
		cout << endl << "Temperature of the system: " << temp;

		double beta = 1./temp;

		for(unsigned int j=0; j<t_ind; j++){
			path = g.Mutation(path, cities, beta);
			g.Check(path);
		}

		g.Print_L(path, cities, i);
		temp -= ((t_start-t_end)/t_gen);

	}

	g.Print_Path(path, cities);
					
	return 0;

}
