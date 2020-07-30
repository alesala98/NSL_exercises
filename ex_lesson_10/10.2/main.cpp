#include <iostream>
#include <mpi.h>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

#include "random.h"
#include "class.h"

using namespace std;

int main (int argc, char *argv[]){

	int size, rank;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Status stat1, stat2, stat3, stat4;
	MPI_Request req;
	
	// Initialization Random Generator
   	Random rnd;
	Init(rnd,rank);
	rnd.SaveSeed();

	// EXERCISE 10.2
	unsigned int ngen = 2000, nind = 500, ncit = 32, nmigr = 10;		// Number of generations, individuals, cities and migration
	double p = 2.;								// Rigged roulette probability exponent
	
	cout << endl << "Parallel traveling salesman problem " << rank << endl << endl;
	cout << "Number of generations:     " << ngen << endl;
	cout << "Number of individuals:     " << nind << endl;
	cout << "Number of cities:          " << ncit << endl;
	
	// Initialization: creating world and allocating needed classes/variables
	cout << endl << "Creating world..." << endl;

	Generation g(ncit,nind,rank,rnd);

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

	double tstart = MPI_Wtime();

	// Simulation
	for(unsigned int i=0; i<ngen; i++){
		
		cout << endl << "Node: " << rank << " generation: " << i;
		individuals = g.Sort(individuals, cities);

		// Exchange best individuals of all different nodes
		if((i+1)%nmigr == 0){	

			cout << endl << "Exchanging best individuals ... " << endl;

			int dim = individuals[0].size(), isend;
			vector<int> v(3);
			iota(v.begin(),v.end(),1);

			if(rank == 0)										// 'root' node has rank == 0. It chooses randomly which node ('isend') 
				isend = rnd.Rannyu(1,3);							// to communicate with. Then, using MPI_Bcast, root sends the information  
														// to all the other nodes (*).  Therefore, the communication between root 
			MPI_Bcast(&isend, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);	// (*)				// and isend is established, as well as for the other two nodes. In particular, 
														// these other nodes see each other using a vector v which is filled with the 
			cout << endl << "Node " << rank << " sending to node " << isend << endl;		// ranks of all the nodes that are neither root nor isend (**). Note that this 
			if(rank == 0){										// implementation works quite well for 4 different nodes, but if you increase 
														// the number of nodes you have to make some changes in the for loop (**) in order 
				vector<int> tosend = individuals[0];						// to assure that the nodes that are neither root nor isend chose one another randomly
				vector<int> receive(dim);

				MPI_Isend(&tosend[0], ncit, MPI_INTEGER, isend, 0, MPI_COMM_WORLD, &req);
				MPI_Recv(&receive[0], ncit, MPI_INTEGER, isend, 1, MPI_COMM_WORLD, &stat2);
		
				individuals[0] = receive;

			} else if(rank == isend){

				cout << endl << "Node " << rank << " sending to node " << 0 << endl;

				vector<int> tosend = individuals[0];
				vector<int> receive(dim);

				MPI_Send(&tosend[0], ncit, MPI_INTEGER, 0, 1, MPI_COMM_WORLD);
				MPI_Recv(&receive[0], ncit, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, &stat1);
		
				individuals[0] = receive;
		
			}

			v.erase(find(v.begin(),v.end(),isend)); 		
			
			for(unsigned int irank=0; irank<v.size(); irank++){	// (**)		// Communication between the nodes that are neither root nor isend
				
				if(rank == v[irank]){
		
					cout << endl << "Node " << rank << " sending to node " << v[irank+1] << endl;

					vector<int> tosend = individuals[0];
					vector<int> receive(dim);

					MPI_Isend(&tosend[0], ncit, MPI_INTEGER, v[irank+1], irank+1, MPI_COMM_WORLD, &req);
					MPI_Recv(&receive[0], ncit, MPI_INTEGER, v[irank+1], irank, MPI_COMM_WORLD, &stat4);
		
					individuals[0] = receive;
					
				} else if(rank == v[irank+1]){
					
					cout << endl << "Node " << rank << " sending to node " << v[irank] << endl;

					vector<int> tosend = individuals[0];
					vector<int> receive(dim);

					MPI_Send(&tosend[0], ncit, MPI_INTEGER, v[irank], irank, MPI_COMM_WORLD);
					MPI_Recv(&receive[0], ncit, MPI_INTEGER, v[irank], irank+1, MPI_COMM_WORLD, &stat3);
		
					individuals[0] = receive;
					
				}

			v.clear();

			}
			
		} 	// Exchange terminated

		// Mutations and crossover
		for(unsigned int j=0; j<nind; j+=2){
			
			new_born = g.Crossover(individuals[g.Roulette(p)], individuals[g.Roulette(p)]);
			new_born = g.Mutation(new_born[0], new_born[1]);
			
			g.Check(new_born[0]);
			g.Check(new_born[1]);

			new_individuals[j] = new_born[0];
			new_individuals[j+1] = new_born[1];

		}

		individuals = new_individuals;
		g.Print_L(individuals, cities, i, rank);

	}

	double tend = MPI_Wtime();

	cout << endl << "Node " << rank << " simulation time: " << tend-tstart << endl;

	g.Print_Path(individuals, cities, rank);
		
	MPI_Finalize();	
		
	return 0;

}
