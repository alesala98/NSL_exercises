#ifndef __Class__
#define __Class__

using namespace std;

// Classes

class City{

	private:
		double m_x, m_y;

	public:
  		City(double, double);		//Constructor
        	~City(){};     			//Destructor
						
 		void Set_Param(double, char);	//Methods
		double Get_Param(char);
  		
};

class Generation{

	private:
		unsigned int m_ncity, m_individual;
		Random m_rnd;
		
	public:
		Generation(unsigned int, unsigned int, Random);		//Constructor
		~Generation(){};					//Destructor

		vector<int> Gen_Individual();				//Methods
		void Check(vector<int>);

		vector<int> Pair_Permutation(vector<int>);
		vector<int> Block_Permutation(vector<int>);
		vector<int> Shift(vector<int>);
		vector<int> Inversion(vector<int>);

		vector<vector<int>> Mutation(vector<int>, vector<int>);
		vector<vector<int>> Crossover(vector<int>, vector<int>);
		vector<vector<int>> Sort(vector<vector<int>>, vector<City>);

		double Norm(vector<int>, vector<City>, int);

		unsigned int Roulette(double);

		void Print_L(vector<vector<int>>, vector<City>, unsigned int);
		void Print_Path(vector<vector<int>>, vector<City>);
		
};

// Functions

vector<City> Generate_Cities(unsigned int, char, Random);


#endif // __Class__

