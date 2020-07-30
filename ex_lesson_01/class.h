#ifndef __Class__
#define __Class__

using namespace std;

class Statistics {

	private:
		unsigned int m_M, m_N, m_L;

	public:
  		Statistics(unsigned int, unsigned int);		// constructor
        	~Statistics();     				// destructor
						
 		void SetParameters(unsigned int, char);		//methods
  		vector<vector<double>> Ave_StDev(Random, vector<vector<double>>,char);
		vector<vector<double>> Pearson(Random, vector<vector<double>>);
		vector<vector<double>> Central_Lim(Random, vector<vector<double>>, unsigned int, char);
		vector<vector<double>> Buffon(Random, vector<vector<double>>, double, double);
		double Error(vector<double>, vector<double>, double);

		void Print(vector<vector<double>>, string);
		
};	

#endif // __Class__

