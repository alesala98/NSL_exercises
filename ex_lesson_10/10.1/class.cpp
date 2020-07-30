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

vector<int> Generation :: Mutation(vector<int> path, vector<City> cities, double beta){

	vector<int> support;

	unsigned int rand = m_rnd.Rannyu(0,4);

	if(rand == 0)			// Pair permutation
		support = Pair_Permutation(path);

	if(rand == 1)			// Shift
		support = Shift(path);

	if(rand == 2)			// Block permutation
		support = Block_Permutation(path);

	if(rand == 3)			// Inversion
		support = Inversion(path);
	
	double A = min(1., exp(beta*(Norm(support, cities, 1) - Norm(path, cities, 1))));

	if(A >= m_rnd.Rannyu())
		support = path;

	return support;		
	
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

void Generation :: Print_L(vector<int> path, vector<City> cities, unsigned int t_gen){

	ofstream Best_L1;
	Best_L1.open("../Datas/best_L1.dat", ios::app);

	Best_L1 << t_gen << " " << Norm(path,cities,1) << endl;

	Best_L1.close();


}

void Generation :: Print_Path(vector<int> path, vector<City> cities){

	ofstream Path;
	Path.open("../Datas/path.dat", ios::app);

    	for(unsigned int i=0; i<m_ncity; i++)
        	Path << cities[path[i]].Get_Param('x') << " " << cities[path[i]].Get_Param('y') << endl;

    	Path << cities[path[0]].Get_Param('x') << " " << cities[path[0]].Get_Param('y') << endl;			// Print again first city in oder to 
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


