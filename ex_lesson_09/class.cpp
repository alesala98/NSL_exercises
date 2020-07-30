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

// Class City Methods

City :: City(double x, double y){

	m_x = x;
	m_y = y;

}	
				
void City :: Set_Param(double n, char c){	

	switch(c){
		case 'x':
			m_x = n;
			break;
		case 'y':
			m_y = n;
			break;
		default:
			cerr << endl << "Error: no such parameter exists" << endl;
			break;
	}

}

double City :: Get_Param(char c){

	switch(c){
		case 'x':
			return m_x;
		case 'y':
			return m_y;
		default:
			cerr << endl << "Error: no such parameter exists" << endl;	
			exit(0);
	}

}

// Class Generation methods
 
Generation :: Generation(unsigned int nc, unsigned int nind, Random r){

	m_ncity = nc;
	m_individual = nind;

	m_rnd = r;

}

vector<int> Generation :: Gen_Individual(){

	vector<int> ind(m_ncity);

    	for(unsigned int i=0; i<m_ncity; i++)
        	ind[i] = i;

    	random_shuffle(ind.begin()+1, ind.end());

	return ind;

}

void Generation :: Check(vector<int> ind){
	
	for(unsigned int i=0; i<ind.size(); i++)
		if(ind[i]<0 or ind[i]>(int)ind.size()){
			cerr << "Wrong value ... Computation aborted" << endl << endl;			
			exit(0);
		}

	for(unsigned int i=0; i<ind.size(); i++)
		for(unsigned int j=i+1; j<ind.size(); j++)
			if(ind[i] == ind[j]){
				cerr << "Duplicated value ... Computation aborted" << endl << endl;
				exit(1);
			}

}

vector<int> Generation :: Pair_Permutation(vector<int> ind){

	unsigned int i = m_rnd.Rannyu(1, ind.size());
	unsigned int j;

	do{
		 j = m_rnd.Rannyu(1, ind.size());
	}while(j == i);	
	
	swap(ind[i],ind[j]);

	return ind;

}

vector<int> Generation :: Block_Permutation(vector<int> ind){

	unsigned int mid = ind.size()/2;
	unsigned int k = m_rnd.Rannyu(1, mid);
	unsigned int l = m_rnd.Rannyu(k, mid);

	vector<int> newind = ind;

	for(unsigned int i=k; i<l; i++){
		newind[i] = ind[i+mid];
		newind[i+mid] = ind[i];
	}

	return newind;

}

vector<int> Generation :: Shift(vector<int> ind){

	unsigned int k = m_rnd.Rannyu(1, ind.size());

	rotate(ind.begin()+1, ind.begin() + k, ind.end());

	return ind;
	
}

vector<int> Generation :: Inversion(vector<int> ind){

	unsigned int k = m_rnd.Rannyu(1, ind.size());
	unsigned int l = m_rnd.Rannyu(k, ind.size());

	vector<int> newind = ind;

	for(unsigned int i=0; i<(l-k); i++)
		newind[i+k] = ind[l-i-1];

	return newind;

}

unsigned int Generation :: Roulette(double p){

		unsigned int i;

		do{
			i = (m_individual*pow(m_rnd.Rannyu(), p)) + 1;
		}while(i > m_individual-1);

		return i;

}

vector<vector<int>> Generation :: Mutation(vector<int> ind1, vector<int> ind2){

	vector<vector<int>> inds(2);

	inds[0] = ind1;
	inds[1] = ind2;

	for(unsigned int i=0; i<2; i++){			// Mutation probability must be < 10%

		if(m_rnd.Rannyu() < 0.1)			// Pair permutation
			inds[i] = Pair_Permutation(inds[i]);

		if(m_rnd.Rannyu() < 0.1)			// Shift
			inds[i] = Shift(inds[i]);

		if(m_rnd.Rannyu() < 0.1)			// Block permutation
			inds[i] = Block_Permutation(inds[i]);

		if(m_rnd.Rannyu() < 0.1)			// Inversion
			inds[i] = Inversion(inds[i]);
	
	}
	
	return inds;		
	
}

vector<vector<int>> Generation :: Crossover(vector<int> parent1, vector<int> parent2){

	vector<vector<int>> inds(2);

	inds[0] = parent1;		// Copy all elements of parent1 in one of the two sons. Then procede to crossover.
	inds[1] = parent2;

	if(m_rnd.Rannyu() > 0.3){					// Crossover probability must be at least > 50%

		unsigned int j = m_rnd.Rannyu(1, parent1.size());	// Cut position

		vector<int> support1 = parent2;				// Copy all elements of parent_x in vector support_y
        	vector<int> support2 = parent1;

        	for(unsigned int i=0; i<support1.size(); i++){						// Function 'find' searches in vector 'parent1' the elements contained in 'parent2' 
													// between the cut position j and the end of vector 'parent1'. If an element is not 
            		if(find(parent1.begin() + j, parent1.end(), support1[i]) == parent1.end()){	// found (i.e. 'find' function gives back the iterator to last element of 'parent1')
                		support1.erase(support1.begin()+i);					// it means that such element is already stored in vector 'parent1' before the cut j. 
                		i--;									// Therefore, there is no need to preserve it in vector 'support1' and the script erases it 
            		}										// changing the size of vector 'support1'. Otherwise, once the script reaches the third for 
													// loop (*), it duplicates this value which causes the Check method to abort the simulation,
		}											// which is not a good thing at all.

		for(unsigned int i=0; i<support2.size(); i++){						//Same thing as above, this time done for parent2

            		if(find(parent2.begin() + j, parent2.end(), support2[i]) == parent2.end()){
                		support2.erase(support2.begin()+i);
                		i--;
            		}

        	}

       		for(unsigned int i=0; i<support1.size(); i++){		// (*) 
            		inds[0][i+j]=support1[i];			// Cut parents at same position, conserve first part and complete second 
            		inds[1][i+j]=support2[i];			// one with missing cities of the consort stored in vector support_x
		}

        }

	return inds;

}

double Generation :: Norm(vector<int> ind, vector<City> C, int n){

	double L = 0.;

	if(n == 1){		//Module-Norm

		for(unsigned int i=0; i<C.size()-1; i++)
			L += sqrt(pow(C[ind[i]].Get_Param('x') - C[ind[i+1]].Get_Param('x'), 2.) + pow(C[ind[i]].Get_Param('y') - C[ind[i+1]].Get_Param('y'), 2.));

		L += sqrt(pow(C[C.size()-1].Get_Param('x') - C[0].Get_Param('x'), 2.) + pow(C[C.size()-1].Get_Param('y') - C[0].Get_Param('y'), 2.));	// Adding path from N_City to 1_City 

	}

	if(n == 2){		//Square-Norm

		for(unsigned int i=0; i<C.size()-1; i++)
			L += (pow(C[ind[i]].Get_Param('x') - C[ind[i+1]].Get_Param('x'), 2.) + pow(C[ind[i]].Get_Param('y') - C[ind[i+1]].Get_Param('y'), 2.));

		L += (pow(C[C.size()-1].Get_Param('x') - C[0].Get_Param('x'), 2.) + pow(C[C.size()-1].Get_Param('y') - C[0].Get_Param('y'), 2.));	// Adding path from N_City to 1_City 
	
	}

	return L;

}

vector<vector<int>> Generation :: Sort(vector<vector<int>> inds, vector<City> cities){		// Sort algorithm first computes the distances for all the individuals (vector dist) (*). Then,
												// it creates a vector of indexes (support), which is sorted comparing the previously computed 
	vector<double> dist(inds.size());							// distances (**). Finally, input vector inds is sorted following the order given by the indexes 
												// vector: this  is done copying inds in newinds (***).
    	for(unsigned int i=0; i<inds.size(); i++)		// (*)
        	dist[i] = Norm(inds[i], cities, 1);		// Choose '1' to select 1-Norm
								// Choose '2' to select 2-Norm
	vector<int> support(dist.size());
	iota(support.begin(), support.end(), 0);
   	
	sort(support.begin(), support.end(), [&dist](size_t i1, size_t i2) {return dist[i1] < dist[i2];});	// (**)

    	vector<vector<int>> newinds(inds.size());

    	for(unsigned int i=0; i<inds.size(); i++)		// (***)
      		newinds[i] = inds[support[i]];
    	
    	return newinds;						// The output of this algorithm is a vector newinds that stores in the first position
								// the individual representing the best path, i.e. the path that minimizes L^1 (or L^2)
}

void Generation :: Print_L(vector<vector<int>> inds, vector<City> cities, unsigned int ngen){

	ofstream Ave_L1, Best_L1;
	Ave_L1.open("Datas/ave_L1.dat", ios::app);
	Best_L1.open("Datas/best_L1.dat", ios::app);

	double dist1 = 0.;

	vector<vector<int>> support = Sort(inds, cities);

	Best_L1 << ngen << " " << Norm(support[0],cities,1) << endl;

	for(unsigned int i=0; i<support.size()/2; i++)
		dist1 += Norm(support[i],cities,1);
	
	Ave_L1 << ngen << " " << 2*dist1/support.size() << endl;

	Ave_L1.close();
	Best_L1.close();


}

void Generation :: Print_Path(vector<vector<int>> inds, vector<City> cities){

	ofstream Path;
	Path.open("Datas/path.dat", ios::app);

	vector<vector<int>> support = Sort(inds, cities);

    	for(unsigned int i=0; i<m_ncity; i++)
        	Path << cities[support[0][i]].Get_Param('x') << " " << cities[support[0][i]].Get_Param('y') << endl;

    	Path << cities[support[0][0]].Get_Param('x') << " " << cities[support[0][0]].Get_Param('y') << endl;		// Print again first city in oder to 
															// draw the path from N_City to 1_City
	Path.close();

}


// Functions

vector<City> Generate_Cities(unsigned int ncities, char c, Random rnd){

	vector<City> cities;

	if(c == 'c'){			// Cities generated on a circumference

		for(unsigned int i=0; i<ncities; i++){
			double alfa = rnd.Rannyu(0.,2*M_PI);
        		City C = City(cos(alfa), sin(alfa));
        		cities.push_back(C);
		}

	}
	else if(c == 's'){		// Cities generated inside a square

		for(unsigned int i=0; i<ncities; i++){
        		City C = City(rnd.Rannyu(), rnd.Rannyu());
        		cities.push_back(C);
		}

	}

	return cities;

}


