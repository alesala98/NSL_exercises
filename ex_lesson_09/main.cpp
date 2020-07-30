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

	// EXERCISE 9.1
	unsigned int ngen = 2000, nind = 300, ncit = 32;		// Number of generations, individuals and cities
	double p = 2.;							// Rigged roulette probability exponent
	
	cout << endl << "Traveling salesman problem" << endl << endl;
	cout << "Number of generations:     " << ngen << endl;
	cout << "Number of individuals:     " << nind << endl;
	cout << "Number of cities:          " << ncit << endl;
	
	// Initialization: creating world and allocating needed classes/variables
	cout << endl << "Creating world..." << endl;

	Generation g(ncit,nind,rnd);

	vector<City> cities = Generate_Cities(32,'s',rnd);		// Choose 'c' if you want cities on a circumference or 's' if you want them inside a 
									// square. Further details can be found in the attached 'class.cpp' and 'class.h' files
	vector<vector<int>> individuals(nind);
	vector<vector<int>> new_individuals(nind);
	vector<vector<int>> new_born(2);

	for(unsigned int i=0; i<nind; i++){
		individuals[i] = g.Gen_Individual();
		g.Check(individuals[i]);
	}

	cout << endl << "World created, starting simulation... " << endl;

	// Simulation
	for(unsigned int i=0; i<ngen; i++){
		
		cout << endl << "Generation: " << i;
		individuals = g.Sort(individuals, cities);

		for(unsigned int j=0; j<nind; j+=2){
			
			new_born = g.Crossover(individuals[g.Roulette(p)], individuals[g.Roulette(p)]);
			new_born = g.Mutation(new_born[0], new_born[1]);
			
			g.Check(new_born[0]);
			g.Check(new_born[1]);

			new_individuals[j] = new_born[0];
			new_individuals[j+1] = new_born[1];

		}

		individuals = new_individuals;
		g.Print_L(individuals, cities, i);

	}

	g.Print_Path(individuals, cities);
					
	return 0;

}
