#ifndef __Class__
#define __Class__

using namespace std;

class Statistics {

	private:
		unsigned int m_M, m_N, m_L;
		Random m_rnd;

	public:
  		Statistics(unsigned int, unsigned int, Random);		//Constructor
        	~Statistics();     					//Destructor
						
 		void Set_Iterations(unsigned int, char);		//Methods
		void Optimization(bool);
  		vector<vector<double>> Metropolis(vector<vector<double>>);
		double ProbWaveFunc(double, double, double, char);
		double V(double);
		double Error(vector<double>, vector<double>, double);
		void Print(vector<vector<double>>, string);
		
};	

#endif // __Class__

