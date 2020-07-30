#ifndef __Class__
#define __Class__

using namespace std;

class Statistics {

	private:
		unsigned int m_M, m_N, m_L;

	public:
  		Statistics(unsigned int, unsigned int);		// constructor
        	~Statistics();     				// destructor
						
 		void Set_Iterations(unsigned int, char);						//methods
  		vector<vector<double>> Metropolis(Random, vector<vector<double>>, char, char, double);
		double ProbWaveFunc(vector<double>, char);
		double Error(vector<double>, vector<double>, double);
		void Print(vector<vector<double>>, string);
		
};	

#endif // __Class__

