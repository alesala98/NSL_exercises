#ifndef __Class__
#define __Class__

using namespace std;

class Statistics {

	private:
		unsigned int m_M, m_N, m_L;
		double m_S0, m_T, m_K, m_r, m_sigma;

	public:
  		Statistics(unsigned int, unsigned int);		// constructor
        	~Statistics();     				// destructor
						
 		void Set_Iterations(unsigned int, char);						//methods
		void Set_Parameters(double, char);
  		vector<vector<double>> European_Option(Random, vector<vector<double>>, unsigned int, char);
		double Error(vector<double>, vector<double>, double);
		void Print(vector<vector<double>>, string);
		
};	

#endif // __Class__

