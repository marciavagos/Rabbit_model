#include "MersenneTwister.h"
#include <time.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#ifndef GENPURPOSE_HRDBARLC_H
#define GENPURPOSE_HRDBARLC_H


using namespace std;

class GenPurpose
{	
		
	public:	

	static double sum(double* vals, int l)
	{
		double sumval = 0;
		for(int i=0; i<l; i++)
		{
			sumval += vals[i];
		}
		
		return sumval;
	}
	
	static double average(double* vals, int l)
	{
		return sum(vals, l) / l;
	}
	
	static double std(double* vals, int le)
	{
		double avg = average(vals, le);
		double sumval = 0;
		for(int i=0; i<le; i++)
		{
			sumval += ((vals[i] - avg) * (vals[i] - avg));
		}
		return safeSqrt(sumval / (le-1));
	}
	
	static bool occursInList(int* list, int listlength, int val)
	{
		for(int i=0; i<listlength; i++)
		{
			if(list[i] == val) return true;
		}
		return false;
	}
	
	static bool occursInList(double* list, int listlength, double val)
	{
		for(int i=0; i<listlength; i++)
		{
			if(abs(list[i] - val) < 1E-10) return true;
		}
		return false;
	}
	
	static bool parseBoolVal(string val)
	{
		if(sCompI(val,"True") || val == "1") return true;
		return false;
	}
	
	static bool sCompI(string sone, string stwo)
	{
		if(sone.length() != stwo.length()) return false;
	
		for(int i=0; i<sone.length(); i++)
		{
			if(toupper(sone[i]) != toupper(stwo[i])) return false;
		}
		return true;
	}

	static int round(double val)
	{
		//slightly biased towards positive infinity since 10.5 is rounded to 11, etc.
		return (int)floor(val + 0.5);
	}
	
	static double safeSqrt(double val)
	{
		if(val < 0) return 0;
		return sqrt(val);
	}

	static void wait( double seconds )
	{
		clock_t endwait;
		endwait = clock () + (long)(seconds * CLOCKS_PER_SEC) ;
		while (clock() < endwait) {}
	}

	static void trim(string& StringToModify)
	{   
		if(StringToModify.empty()) return;   
	
		int startIndex = StringToModify.find_first_not_of(" ");   
		int endIndex = StringToModify.find_last_not_of(" ");   
		string tempString = StringToModify;   
		StringToModify.erase();   
		StringToModify = tempString.substr(startIndex, (endIndex-startIndex+ 1) );
	}

	static int countOccurences(string search, string target)
	{
		string temp = search;
		int c = 0;
		bool done = false;

		while(!done)
		{
			int pos = temp.find(target);
			if(pos < 0)
			{
				done = true;
			}
			else
			{
			c++;
				temp = temp.substr(pos+target.length(), temp.length() - pos - target.length() + 1);		
			}
		}

		return c;
	}

	static string expandArrayRanges(string source)
	{
		string result = source;
		bool done = false;
		
		while(!done)
		{
			//cout << "AAA: Result = " << result << endl;
			int pos = result.find("-");
			if(pos < 0)
			{
				done = true;
			}
			else
			{
				int startnum, endnum;
				string leftside = result.substr(0,pos);
				string rightside = result.substr(pos + 1, result.length() - pos);
				
				int leftopt = leftside.rfind(",");
				if(leftopt < 0)
				{
					startnum = (int)atof(leftside.c_str());
					leftside = "";
				}
				else
				{
					startnum = (int)atof((leftside.substr(leftopt + 1, leftside.length() - leftopt)).c_str());
					leftside = leftside.substr(0, leftopt);
					leftside.append(",");
				}
			
				int rightopt = rightside.find(",");
				if(rightopt < 0)
				{
					endnum = (int)atof(rightside.c_str());
					rightside = "";
				}
				else
				{
					endnum = (int)atof((rightside.substr(0,rightopt)).c_str());
					rightside = rightside.substr(rightopt + 1, rightside.length() - rightopt);
					string temp = ",";
					temp.append(rightside);
					rightside = temp;
				}
				
				//cout << "AAA: Startnum = " << startnum << ", Endnum = " << endnum << endl;
				string replacewith = "";
				for(int iii=startnum; iii<=endnum; iii++)
				{
					char buffer[7];
					sprintf(buffer, "%d,", iii);
					replacewith.append(buffer);
				}
				replacewith = replacewith.substr(0, replacewith.length() - 1);
				
				result = leftside.append(replacewith.append(rightside));
			}
		}
		
		return result;
	}
	
	static double* parseArray(string source, string seperat)
	{
		int c = countOccurences(source, seperat) + 1;
		double* res = new double[c];

		string temp = source;
		for(int i =0; i < c - 1; i++)
		{
			int pos = temp.find(seperat);
			string item = temp.substr(0, pos);
			res[i] = atof(item.c_str());
			temp = temp.substr(pos + 1, temp.length() - pos);		
		}

		res[c-1] = atof(temp.c_str());

		return res;
	}

	static int* parseIntArray(string source, string seperat)
	{
		int c = countOccurences(source, seperat) + 1;
		int* res = new int[c];

		string temp = source;
		for(int i =0; i < c - 1; i++)
		{
			int pos = temp.find(seperat);
			string item = temp.substr(0, pos);
			res[i] = atoi(item.c_str());
			temp = temp.substr(pos + 1, temp.length() - pos);		
		}

		res[c-1] = atoi(temp.c_str());

		return res;
	}
	
	static string* parseStringArray(string source, string seperat)
	{
		int c = countOccurences(source, seperat) + 1;
		string* res = new string[c];

		string temp = source;
		for(int i =0; i < c - 1; i++)
		{
			int pos = temp.find(seperat);
			string item = temp.substr(0, pos);
			res[i] = item.c_str();
			temp = temp.substr(pos + 1, temp.length() - pos);		
		}

		res[c-1] = temp.c_str();

		return res;
	}
	
	static double** parseMatrix(string source)
	{
		int numrows = countOccurences(source, ";") + 1;
		double** res = new double*[numrows];
		
		string temp = source;
		for(int i =0; i < numrows - 1; i++)
		{
			int pos = temp.find(";");
			string item = temp.substr(0, pos);
			
			res[i] = parseArray(item);
			temp = temp.substr(pos + 1, temp.length() - pos);		
		}

		res[numrows-1] = parseArray(temp.c_str());
		
		return res;
	}
	
	static double* parseArray(string source)
	{
		return parseArray(source, ",");
	}
	
	static int* parseIntArray(string source)
	{
		return parseIntArray(source, ",");
	}
	
	static int applyBinomialInverseTransform(MTRand *r, int n, double p)
	{
		if(abs(1.0/sqrt(n) * (sqrt((1-p)/p) - sqrt(p / (1 - p)))) > 0.3) // Box, Hunter and Hunter (1978). Statistics for experimenters. Wiley. p. 130.
		{																 // Criterion for approximation with normal dist according to Wikipedia.
			//cout << "1";
			double rv;
			#pragma omp critical
			{
				rv = r->rand();
			}
			double p_k = 1.0;
			double n_over_k = 1;
			double mp_nmk = pow(1.0 - p, n);
			double probsum = 0;
		
			//cout << " n = " << n << ", p = " << p << endl;
			for(int k=0; k<=n; k++)
			{			
				//cout << " k = " << k << "; P(k) = " << n_over_k << " * " << p_k << " * " << mp_nmk << " = " << n_over_k * p_k * mp_nmk << endl;			
				probsum += n_over_k * p_k * mp_nmk;
				if(probsum > rv) return k;
			
				p_k = p_k * p;
				mp_nmk = mp_nmk / (1.0 - p);
				n_over_k = (n - k) / (k+1.0) * n_over_k;
			}
		
			return n;
		}
		else
		{
			// Use normal approximation
			//cout << "2";
			int res;
			#pragma omp critical
			{
				res = round(r->randNorm(n*p, sqrt(n*p*(1-p))));
			}			
			
			if(res < 0 || res > n)
			{
				//cout << "3";
				double rv = r->rand();
				double p_k = 1.0;
				double n_over_k = 1;
				double mp_nmk = pow(1.0 - p, n);
				double probsum = 0;
			
				//cout << " n = " << n << ", p = " << p << endl;
				for(int k=0; k<=n; k++)
				{			
					//cout << " k = " << k << "; P(k) = " << n_over_k << " * " << p_k << " * " << mp_nmk << " = " << n_over_k * p_k * mp_nmk << endl;			
					probsum += n_over_k * p_k * mp_nmk;
					if(probsum > rv) return k;
				
					p_k = p_k * p;
					mp_nmk = mp_nmk / (1.0 - p);
					n_over_k = (n - k) / (k+1.0) * n_over_k;
				}
				
				return n;
			}
			else
			{
				return res;
			}
		}
	}	
};

#endif //GENPURPOSE_HRDBARLC_H