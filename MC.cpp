#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
#include <deque>
#include <sys/time.h>

using namespace std;
class MC{
public:
	double minpot1( vector<int> &v, vector<double> &r ){
		double m;
		m = -r[v.back()];
		v.pop_back();

		while(v.empty() != true){
			m = min(m, -r[v.back()]);
			v.pop_back();
		}
		return m;
	}
	double minpot2( vector<int> &v, vector<double> &r ){
		double m;
		
		m = r[v.back()];
		v.pop_back();

		while(v.empty() != true){
			m = min(m, r[v.back()]);
			v.pop_back();
		}
		return m;
	}


	void initflow( vector<int> &x, vector<double> const &r, vector<int> const &c, int A){
		int i;
		for(i=0; i < A; i++){
			if( r[i] >= 0){
			x[i] = 0;
			}else{
			x[i] = c[i];
			}
		}
	}

	double cost( vector<int> const &x, vector<int> const &d, int A){
		int i;
		double totalcost = 0;
		for(i = 0; i < A; i++){
			totalcost += double(x[i])*double(d[i]);
		}
		return totalcost;
	}

	bool checkfeas( vector<int> const &x, vector<int> const &c, vector<int> const &s, int A, int N){
		int i;
		for(i=0; i < A; i++){
			if( (x[i] > c[i]) || (x[i]<0) ){ return false; }
		}
		for(i=0; i < N; i++){
			if(s[i] != 0){ return false; }
		}
		return true;
	}

	void printflow( vector<int> &x, int A){
		int i;
		for(i=0; i < A; i++){
			cout<<"x["<<i<<"] = "<<x[i]<<endl;
		}
	}

	double diff_sec(timeval t1, timeval t2){
		return (t2.tv_sec - t1.tv_sec +  (t2.tv_usec - t1.tv_usec)/1000000.0 );
	}


	void readdata(const char *filename, int &N, int &A, vector<int> &startnode, vector<int> &endnode, vector<int> &d, vector<int> &c, vector<int> &b, vector<int> &degree){
		ifstream inputdata;
		int databuffer; 
		int i;
		inputdata.open(filename);
		if(!inputdata){
			cerr << "Can't open input file " << filename << endl;
			exit(1);
			}
		inputdata >> N;
		inputdata >> A;

		startnode.resize(A);
		endnode.resize(A);
		d.resize(A);
		c.resize(A);
		b.resize(N);
		degree.resize(N);

		for(i = 0; i < A; i++){
			inputdata >> databuffer;
			startnode[i] = databuffer - 1;
			degree[databuffer-1] += 1;
			
			inputdata >> databuffer;
			endnode[i] = databuffer - 1;
			degree[databuffer-1] += 1;
			
			inputdata >> d[i];
			inputdata >> c[i];
		}
		
		for(i=0; i<N; i++){
			inputdata >> b[i];
		}
		inputdata.close();
	}

	void solve( int argc, char *argv[] ) {

		int i; int j;
		int m; int n;
		int N; 
		int A; 
		const char *filename = argv[1];
		timeval t1;
		timeval t2;
		double dif;

		vector<int> startnode;
		vector<int> endnode;
		vector<int> d;
		double d_max; 
		vector<int> c; 
		vector<int> b; 
		vector<int> degree;
		int maxdeg; 
		int maxdegnode;

		gettimeofday(&t1, NULL);

		readdata(filename, N, A, startnode, endnode, d, c, b, degree);

		maxdeg=0;
		for(j=0; j<N; j++){ 
			if(maxdeg < degree[j]){ 
				maxdegnode = j;
				maxdeg = degree[j];
			}
		}

		vector< vector<int> > arcout(N, vector<int>(maxdeg));  
		vector< vector<int> > arcin(N, vector<int>(maxdeg));

		vector<int> numarcout(N,0);
		vector<int> numarcin(N,0); 

		for(i=0; i < A; i++){
			m = startnode[i];
			n = endnode[i];
			arcin[n][numarcin[n]] = i;
			numarcin[n] += 1;
			arcout[m][numarcout[m]] = i;
			numarcout[m] += 1;
		}


		d_max = fabs(double(d[0]));
		for(i=0; i < A; i++){
			if(fabs(double(d[i]))> d_max){ 
				d_max = fabs(double(d[i])); 
			}
		}

		double epsilon = d_max;
		vector<double> u(N, 0); 
		vector<int> x(A, 0); 
		vector<double> r(A); 
		vector<int> s(N); 
		deque<int> nodeq;
		int ibar;

		for(i=0; i < A; i++){
			r[i] = double(d[i]) + u[startnode[i]] - u[endnode[i]];
		}

		bool proceed;
		int k;
		double alpha;
		double alpha1;
		double alpha2;
		int beta;
		vector<int> temp;
		int arc;

		while(epsilon >= 1/( double(N) )){

			initflow(x,r,c,A);

			for(i=0; i < N; i++){
				s[i] = -b[i];
				for(j=0; j < numarcout[i]; j++){
					k = arcout[i][j];
					s[i] += x[k];
				}
				for(j=0; j < numarcin[i]; j++){
					k = arcin[i][j];
					s[i] -= x[k];
				} 
				if(s[i]>0){
					nodeq.push_back(i);
				}
			}

			while(nodeq.empty() == false){
				ibar = nodeq.back();
				nodeq.pop_back();
				while(s[ibar] > 0){
					proceed = true;
					temp.clear();

					for(j=0; j<numarcout[ibar] && proceed == true; j++){
						k = arcout[ibar][j];
						if( (r[k] <= epsilon) && (r[k] >= (epsilon/2)) && (x[k] > 0) ){ 
							arc = k; 
							proceed = false;
						}
					}

					if(proceed == false ){
						k = arc;
						beta = min( x[k], s[ibar] );
						x[k] -= beta;
						s[ibar] -= beta;
						s[endnode[k]] += beta;

						if( (s[endnode[k]] > 0) && (s[endnode[k]] - beta <= 0) ){ 
							nodeq.push_front(endnode[k]);
						}
					}else{
						for( j=0; j<numarcin[ibar] && proceed==true; j++ ){
							k = arcin[ibar][j];
							if( (r[k] >= (-epsilon)) && (r[k] <= (-epsilon/2)) && (x[k] < c[k]) ){
								arc = k;
								proceed = false;
							}
						}
				
						if(proceed == false ){
							k = arc;
							beta = min( c[k] - x[k], s[ibar]);
							x[k] += beta;
							s[ibar] -= beta;
							s[startnode[k]] += beta;

							if( (s[startnode[k]] > 0) && (s[startnode[k]] - beta) <= 0 ){
								nodeq.push_front(startnode[k]);
							}
						}
					}
					if( proceed==true ){  
						for(j=0; j < numarcout[ibar]; j++){
							k = arcout[ibar][j];
							if( x[k] > 0 ){ temp.push_back(k); }
						}
						if(temp.empty() == false){
							alpha1 = minpot1( temp, r );
							temp.clear();

							for(j=0; j < numarcin[ibar]; j++){
								k = arcin[ibar][j];
								if( x[k] < c[k] ){ temp.push_back(k); }
							}
							if(temp.empty() == false){
								alpha2 = minpot2( temp, r );
								temp.clear();
								alpha = min(alpha1, alpha2) + epsilon;
							}else{
								alpha = alpha1 + epsilon;
							}
						}else{
							for(j=0; j < numarcin[ibar]; j++){
								k = arcin[ibar][j];
								if( x[k] < c[k] ){ temp.push_back(k); }
							}
							alpha = minpot2( temp, r) + epsilon;
							temp.clear();
						}

						u[ibar] = u[ibar] + alpha;

						for(j=0; j < numarcout[ibar]; j++){
							k = arcout[ibar][j];
							r[k] += alpha;
						}
						for(j=0; j < numarcin[ibar]; j++){
							k = arcin[ibar][j];
							r[k] -= alpha;
						}
					}     
				}
			}
			epsilon = epsilon/2;
		}
		gettimeofday(&t2, NULL); 
		dif = diff_sec(t1, t2);
		cout<<"Run time: "<<dif<<" seconds"<<endl;
	
		if( checkfeas(x,c,s,A,N) ){ 
			cout<<"Current flow is feasible"<<endl;
			cout<<"Minimum cost is: "<<cost(x,d,A)<<endl;
		}else{
			cout<<"Current flow not feasible"<<endl;
		}
	}
};

/*
int main( int argc, char *argv[] ){
	MC ac;
	ac.solve(argc, argv);
	return 0;
}
*/    
