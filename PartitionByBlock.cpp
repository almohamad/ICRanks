#include<cstdlib>
#include<algorithm>
#include<cmath>
#include<string>
#include<ctime>
#include<vector>
#include<Rcpp.h>

using namespace Rcpp;

double Summation(const NumericVector& x, const NumericVector& sigma, const int& LowI, const int& UppI)
{
	// LowI should be by default equal to 0
	// UppI should be by default equal to SampleSize - 1
	double SumRes = x[LowI]/(sigma[LowI]*sigma[LowI]);
	double SumSigma = 1/(sigma[LowI]*sigma[LowI]);
	for (int i = LowI + 1; i <= UppI; i++)
	{
		SumRes += x[i]/(sigma[i]*sigma[i]);
		SumSigma += 1/(sigma[i]*sigma[i]);
	}
	return SumRes/SumSigma;
}

double LogLikelihood(const NumericVector& y, const NumericVector& sigma, const int& LowI, const int& UppI) // Calculate the LogLiklihood for points with indices between LowI and UppI
{
	// Calculate the average
	double MeanY = Summation(y, sigma, LowI, UppI); // division by the number of elements is included
	double LogLikhood = 0;
	for (int i = LowI; i <= UppI; i++)
	{
		LogLikhood += (y[i] - MeanY)*(y[i] - MeanY) / (sigma[i] * sigma[i]);
	}
	return LogLikhood;
}

void BinaryConfig(unsigned long long int c, int*& Config, int& l, const int& Shift, const int& start)
{
	unsigned long long int residu = c;
	int counter = 0; l = 0;
	while (residu > 1)
	{
		if (residu % 2 == 1)
		{
			Config[l+start] = counter+Shift;
			l++;
		}

		residu /= 2;
		counter++;
	}
	if (residu == 1)
	{
		Config[l+start] = counter+Shift;
		l++;
	}
}
void RankUpdate(IntegerVector& Lower, IntegerVector& Upper, const int* InqPosi, const int& l, const int& n)
{
	for (int i = 0; i <= InqPosi[0]; i++)
	{
		Lower[i] = 0;
		if (InqPosi[0] > Upper[i]) {
			Upper[i] = InqPosi[0];
		}
	}
	int j = 0;
	
	while (j <= (l - 2))
	{
		for (int i = InqPosi[j] + 1; i <= InqPosi[j + 1]; i++)
		{
			if (InqPosi[j] + 1 < Lower[i]) {
				Lower[i] = InqPosi[j] + 1;
			}
			if (InqPosi[j + 1] > Upper[i]) {
				Upper[i] = InqPosi[j + 1];
			}
		}
		j++;
	}
	for (int i = InqPosi[l - 1] + 1; i < n; i++)
	{
		if (InqPosi[l - 1] + 1 < Lower[i]) {
			Lower[i] = InqPosi[l - 1] + 1;
		}
		Upper[i] = n - 1;
	}
}
unsigned long long int binomialCoeff(int n, int k)
{
	if (k>n) return 0;
	unsigned long long int res = 1;

	// Since C(n, k) = C(n, n-k)
	if (k > n - k){
		k = n - k;
	}

	// Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
	for (int i = 0; i < k; ++i)
	{
		res *= (n - i);
		res /= (i + 1);
	}

	return res;
}
bool PAVACheck(const NumericVector& y_temp, const NumericVector& sigma_temp, const int& l, const int* InqPosi, const int& n)
{
	if(l==0) return false;
	double MeanBlock1 = Summation(y_temp,sigma_temp,0,InqPosi[0]);
	double MeanBlock2;
	bool CheckPAVA = false;
	int k = 0;
	while (k <= (l - 2))
	{
		MeanBlock2 = Summation(y_temp,sigma_temp,InqPosi[k] + 1,InqPosi[k + 1]);
		if(MeanBlock1>MeanBlock2)
		{
			CheckPAVA = true;
		}
		else
		{
			MeanBlock1 = MeanBlock2;
		}
		k++;
	}
	if(CheckPAVA == false)
	{
		MeanBlock2 = Summation(y_temp,sigma_temp,InqPosi[l - 1] + 1,n-1);
		if(MeanBlock1>MeanBlock2)
		{
			CheckPAVA = true;
		}
	}
	return CheckPAVA;
}
void CorrectPermutationsF(const NumericVector& y, const NumericVector& sigma, const NumericVector& crit, IntegerVector& Lower, IntegerVector& Upper, const int* InqPosi, const int& l, const int& n, const bool& EqSigma)
{
	NumericVector y_temp(n);
	NumericVector sigma_temp(n);
	IntegerVector Lower_temp(n);
	IntegerVector Upper_temp(n);
	
	double LR=0;
	int k = 0;
	bool CheckPAVA = false;
	for(int i=1; i<=n-1; i++)
	{
		int j;
		for(j=1; j<=n-i; j++)
		{
			// Set things into order
			for (int u = 0; u<n; u++)
			{
				Lower_temp[u] = u;
				Upper_temp[u] = u;
				y_temp[u] = y[u];
				sigma_temp[u] = sigma[u];
			}
			// Apply the permutation
			for(int s=j; s<=j+i-1; s++)
			{
	    		y_temp[s] = y[s-1];
	  			sigma_temp[s] = sigma[s-1];
	  			Lower_temp[s] = s-1;
	  			Upper_temp[s] = s-1;
			}
			y_temp[j-1] = y[j+i-1];
			sigma_temp[j-1] = sigma[j+i-1];
			Lower_temp[j-1] = j+i-1;
			Upper_temp[j-1] = j+i-1;
			
			// Test if a PAVA is needed. If yes then break.
			// Could improve this part by using the means in the calculation of the LR later
			CheckPAVA = PAVACheck(y_temp, sigma_temp, l, InqPosi, n);
			if(CheckPAVA == true)
			{
				continue;
			}
			// Calculate the LR corresponding to the permuted vector
			LR = LogLikelihood(y_temp,sigma_temp,0,InqPosi[0]);
			k = 0;
			while (k <= (l - 2))
			{
				LR += LogLikelihood(y_temp,sigma_temp,InqPosi[k] + 1,InqPosi[k + 1]);
				k++;
			}
			LR += LogLikelihood(y_temp,sigma_temp,InqPosi[l - 1] + 1,n-1);		
			if(LR < crit[l])
			{
				RankUpdate(Lower_temp, Upper_temp, InqPosi, l, n);
				// Adjust the CI for the permuted centers
				for(int s=j; s<=j+i-1; s++)
				{
		  			Lower[s-1] = fmin(Lower[s-1], Lower_temp[s]);	
		  			Upper[s-1] = fmax(Upper[s-1], Upper_temp[s]);
				}
				Lower[j+i-1] = fmin(Lower[j+i-1], Lower_temp[j-1]);
				Upper[j+i-1] = fmax(Upper[j+i-1], Upper_temp[j-1]);
				
			}
			else
			{
				if(EqSigma) break;
			}
			
		}
		if(j<n-i && EqSigma) break;
	}

}

void CorrectPermutationsB(const NumericVector& y, const NumericVector& sigma, const NumericVector& crit, IntegerVector& Lower, IntegerVector& Upper, const int* InqPosi, const int& l, const int& n, const bool& EqSigma)
{
	NumericVector y_temp(n);
	NumericVector sigma_temp(n);
	IntegerVector Lower_temp(n);
	IntegerVector Upper_temp(n);
	
	double LR=0;
	int k = 0;
	bool CheckPAVA = false;
	for(int i=1; i<=n-1; i++)
	{
		int j;
		for(j=1; j<=n-i; j++)
		{
			// Set things into order
			for (int u = 0; u<n; u++)
			{
				Lower_temp[u] = u;
				Upper_temp[u] = u;
				y_temp[u] = y[u];
				sigma_temp[u] = sigma[u];
			}
			// Apply the permutation
			for(int s=j; s<=j+i-1; s++)
			{
	    		y_temp[s-1] = y[s];
	  			sigma_temp[s-1] = sigma[s];
	  			Lower_temp[s-1] = s;
	  			Upper_temp[s-1] = s;
			}
			y_temp[j+i-1] = y[j-1];
			sigma_temp[j+i-1] = sigma[j-1];
			Lower_temp[j+i-1] = j-1;
			Upper_temp[j+i-1] = j-1;
			
			// Test if a PAVA is needed. If yes then break.
			// Could improve this part by using the means in the calculation of the LR later
			CheckPAVA = PAVACheck( y_temp, sigma_temp, l, InqPosi, n);
			if(CheckPAVA == true)
			{
				continue;
			}
			// Calculate the LR corresponding to the permuted vector
			LR = LogLikelihood(y_temp,sigma_temp,0,InqPosi[0]);
			k = 0;
			while (k <= (l - 2))
			{
				LR += LogLikelihood(y_temp,sigma_temp,InqPosi[k] + 1,InqPosi[k + 1]);
				k++;
			}
			LR += LogLikelihood(y_temp,sigma_temp,InqPosi[l - 1] + 1,n-1);		
			
			if(LR < crit[l])
			{
				RankUpdate(Lower_temp, Upper_temp, InqPosi, l, n);
				// Adjust the CI for the permuted centers
				for(int s=j; s<=j+i-1; s++)
				{
		  			Lower[s] = fmin(Lower[s], Lower_temp[s-1]);	
		  			Upper[s] = fmax(Upper[s], Upper_temp[s-1]);
				}
				Lower[j-1] = fmin(Lower[j-1], Lower_temp[j+i-1]);
				Upper[j-1] = fmax(Upper[j-1], Upper_temp[j+i-1]);
				/*if(Upper[3]==8) {
					Rcout<<"\n I = "<<i<<", J = "<<j<<"\n InqPosi = (";
					for(int ss = 0; ss<l; ss++)
						Rcout<<InqPosi[ss]<<",";
					Rcout<<")\n";
				}*/
			}
			else
			{
				if(EqSigma) break;
			}
		}
		if(j<n-i && EqSigma) break;
	}
}
void UnrankCombin(int*& S, unsigned long long int m, int k, unsigned long long int**& CnkMat)
{
	int i = k - 1;
	while (i >= 0)
	{
		int l = i;
		while (CnkMat[l][i + 1] <= m)
		{
			l++;
		}
		S[i]= l - 1;
		m = m - CnkMat[l - 1][i + 1];
		i--;
	}

}
// An auxiliary function which does not permute anything. It is expected that the y's are not ordered and that a PAVA check has to be made.
void PartitioningRankingLevel(const NumericVector& y, const NumericVector& sigma, IntegerVector& Lower, IntegerVector& Upper, const NumericVector& crit, unsigned long long int**& CnkMat, const int& n, const bool& trace)
{
	bool CheckPAVA = false;
	double Likelihood0 = 0;
	// Calculate the Likelihood matrix of the blocks.
	double** LikelihoodMat = new double*[n];
	for (int i = 0; i<n; i++)
	{
		LikelihoodMat[i] = new double[n];
		for (int j = i; j<n; j++)
		{
			LikelihoodMat[i][j] = LogLikelihood(y, sigma, i, j);
		}
	}
	
	int* InqPosi = new int[n];
	//if(trace == true) {
	//		Rcout<<"Processed levels:";
	//}
	for (int l = 1; l <= n - 2; l++)
	{
		//			int l = 3; int i = 8;
	//	if(trace == true) Rcout<<l<<".";
		unsigned long long int m = CnkMat[n - 1][l];
		for (unsigned long long int c = 0; c<m; c++)
		{
			UnrankCombin(InqPosi, c, l, CnkMat);
			// Test if a PAVA is needed. If yes then break.
			// Could improve this part by using the means in the calculation of the LR later
			CheckPAVA = PAVACheck(y, sigma, l, InqPosi, n);
			if(CheckPAVA == true)
			{
				continue;
			}
			// Inside each configuration, calculate each group's share in the likelihood
			Likelihood0 = LikelihoodMat[0][InqPosi[0]];
			int j = 0;
			while (j <= (l - 2))
			{
				Likelihood0 += LikelihoodMat[InqPosi[j] + 1][InqPosi[j + 1]];
				j++;
			}
			Likelihood0 += LikelihoodMat[InqPosi[l - 1] + 1][n - 1];
			// Update the ranking
			if (Likelihood0<crit[l])
			{
				RankUpdate(Lower, Upper, InqPosi, l, n);				
			}
		}
		
	}
	delete[] InqPosi;

	/*NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i]+1;
		CIs[n+i] = Upper[i]+1;
	}*/
	
	// Free some space.
	for (int i = 0; i<n; i++)
	{
		delete[] LikelihoodMat[i];
	}
	delete[] LikelihoodMat;
	
	//return CIs;
}


// [[Rcpp::export]]
NumericMatrix PartitioningRankingLevelEqSig(NumericVector y, NumericVector sigma, NumericVector crit, int n, bool trace)
{
	
	// Calculate the Likelihood matrix of the blocks.
	double** LikelihoodMat = new double*[n];
	for (int i = 0; i<n; i++)
	{
		LikelihoodMat[i] = new double[n];
		for (int j = i; j<n; j++)
		{
			LikelihoodMat[i][j] = LogLikelihood(y, sigma, i, j);
		}
	}
	// Calculate the C_n^k matrix in order to use it in the UnrankCombin function
	unsigned long long int ** CnkMat = new unsigned long long int*[n];
	for (int i = 0; i<n; i++)
	{
		CnkMat[i] = new unsigned long long int[n];
		CnkMat[i][i] = 1;
		for (int j = 0; j<i; j++)
		{
			CnkMat[i][j] = binomialCoeff(i, j);
			CnkMat[j][i] = 0;
		}
	}

	// Definition and initialization of the vector of ranks
	IntegerVector Lower(n);
	IntegerVector Upper(n);
	for (int i = 0; i<n; i++)
	{
		Lower[i] = i;
		Upper[i] = i;
	}
	// Test the upper level with all equalities
	double Likelihood0 = LikelihoodMat[0][n - 1];
	if (Likelihood0<crit[0])
	{
		for (int i = 0; i<n; i++)
		{
			Lower[i] = 0;
			Upper[i] = n - 1;
		}
		if(trace == true) {
			Rcout << "Process ended with trivial confidence intervals.\n";
		}
	}
	else
	{
		int* InqPosi = new int[n];
		if(trace == true) {
				Rcout<<"Processed levels:";
		}
		for (int l = 1; l <= n - 2; l++)
		{
			//			int l = 3; int i = 8;
			if(trace == true) Rcout<<l<<".";
			unsigned long long int m = CnkMat[n - 1][l];
			for (unsigned long long int c = 0; c<m; c++)
			{
				UnrankCombin(InqPosi, c, l, CnkMat);
				// Inside each configuration, calculate each group's share in the likelihood
				Likelihood0 = LikelihoodMat[0][InqPosi[0]];
				int j = 0;
				while (j <= (l - 2))
				{
					Likelihood0 += LikelihoodMat[InqPosi[j] + 1][InqPosi[j + 1]];
					j++;
				}
				Likelihood0 += LikelihoodMat[InqPosi[l - 1] + 1][n - 1];
				// Update the ranking
				if (Likelihood0<crit[l])
				{
					RankUpdate(Lower, Upper, InqPosi, l, n);
					CorrectPermutationsF(y, sigma, crit, Lower, Upper, InqPosi, l, n, true);
					CorrectPermutationsB(y, sigma, crit, Lower, Upper, InqPosi, l, n, true);
					
				}
				//	CorrectPermutationsF(y, sigma, crit, Lower, Upper, InqPosi, l, n, EqSigma, NbOfPermut);
				//	CorrectPermutationsB(y, sigma, crit, Lower, Upper, InqPosi, l, n, EqSigma);

			}
			
		}
		delete[] InqPosi;
	}
	NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i]+1;
		CIs[n+i] = Upper[i]+1;
	}
	
	// Free some space.
	for (int i = 0; i<n; i++)
	{
		delete[] LikelihoodMat[i];
	}
	delete[] LikelihoodMat;
	for (int i = 0; i<n; i++)
	{
		delete[] CnkMat[i];
	}
	delete[] CnkMat;
	return CIs;
}


double PartitionCoverage(int***& ResCIsMat, const IntegerVector& Lower_temp, const IntegerVector& Upper_temp, const int& n, const int& MM, const int& a)
{
	int coveragePartition = MM;
	bool coverage;
	for(int j = 0; j<MM; j++)
	{
		coverage = true;
		for(int i = 0; i<n; i++)
		{
			//if(InqPosi[0]==2 && InqPosi[1]==4 && InqPosi[2]==8 && l==3 && trace == true)Rcout<<"["<<ResCIs(i,0)<<","<<ResCIs(i,1)<<"] contains ["<<Lower_temp[i]<<","<<Upper_temp[i]<<"]\n";
			if(ResCIsMat[a][i][2*j]>Lower_temp[i] || ResCIsMat[a][i][2*j+1]<Upper_temp[i])
			{
				coverage = false;
				break;
			}
		}
		//Rcout<<"....\n";
		if(coverage == false)
		{
			coveragePartition--;			
		}
	}
	return 1.0*coveragePartition/(MM*1.0);
}
void PartitioningRankingGeneralProcInit(int***& ResCIsMat, int***& ResCIsGridMat, int*& AlphaRescaled, IntegerVector& Lower, IntegerVector& Upper, const IntegerVector& EmpOrderInit, unsigned long long int**& CnkMat, const NumericMatrix& crit, const int& n, const int& MM, const int& gridSize, const double& alpha, const bool& trace)
{
	double* coveragePartition = new double[gridSize];
	int* InqPosi = new int[n];
	IntegerVector Lower_temp(n), Upper_temp(n), EmpOrder = seq(0,n-1);
	bool coverage;
	double minCoverage;
	int alphaRescaled;
	for (int l = 1; l <= n - 2; l++)
	{
		unsigned long long int m = CnkMat[n - 1][l];
		for (unsigned long long int c = 0; c<m; c++)
		{
			UnrankCombin(InqPosi, c, l, CnkMat);
			Lower_temp = clone(EmpOrder);
			Upper_temp = clone(EmpOrder);
			RankUpdate(Lower_temp, Upper_temp, InqPosi, l, n);

			coveragePartition[0] = fabs(PartitionCoverage(ResCIsMat, Lower_temp, Upper_temp, n, MM, 0)-(1-alpha));
			minCoverage = coveragePartition[0]; alphaRescaled = 0;
			for(int a = 1; a<gridSize; a++)
			{
				coveragePartition[a] = fabs(PartitionCoverage(ResCIsMat, Lower_temp, Upper_temp, n, MM, a)-(1-alpha));
				if(coveragePartition[a]<minCoverage)
				{
					alphaRescaled = a;
					minCoverage = coveragePartition[a];
				}
			}
			AlphaRescaled[(l-1)*CnkMat[n - 1][l-1] + c] = alphaRescaled;
			coverage = true;
			for(int i = 0; i<n; i++)
			{
				if(ResCIsGridMat[alphaRescaled][EmpOrderInit[i]][0]>Lower_temp[i] || ResCIsGridMat[alphaRescaled][EmpOrderInit[i]][1]<Upper_temp[i])
				{
					coverage = false;
					break;
				}
			}
			if(coverage == true)
			{
				for(int s = 0 ; s<n; s++)
				{
					
					Lower[s] = fmin(Lower[s], Lower_temp[s]);
					Upper[s] = fmax(Upper[s], Upper_temp[s]);
					
				}
			}
		}
		
	}
	delete[] InqPosi;
	delete[] coveragePartition;
}
void PartitioningRankingGeneralProc(int***& ResCIsMat, int***& ResCIsGridMat, int*& AlphaRescaled, IntegerVector& Lower, IntegerVector& Upper, const IntegerVector& EmpOrderInit, unsigned long long int**& CnkMat, const NumericMatrix& crit, const int& n, const int& MM, const int& gridSize, const double& alpha, const bool& trace)
{
	int* InqPosi = new int[n];
	IntegerVector Lower_temp(n), Upper_temp(n), EmpOrder = seq(0,n-1);
	bool coverage;
	int alphaRescaled;
	for (int l = 1; l <= n - 2; l++)
	{
		unsigned long long int m = CnkMat[n - 1][l];
		for (unsigned long long int c = 0; c<m; c++)
		{
			UnrankCombin(InqPosi, c, l, CnkMat);
			Lower_temp = clone(EmpOrder);
			Upper_temp = clone(EmpOrder);
			RankUpdate(Lower_temp, Upper_temp, InqPosi, l, n);

			alphaRescaled = AlphaRescaled[(l-1)*CnkMat[n - 1][l-1] + c];
			coverage = true;
			for(int i = 0; i<n; i++)
			{
				if(ResCIsGridMat[alphaRescaled][EmpOrderInit[i]][0]>Lower_temp[i] || ResCIsGridMat[alphaRescaled][EmpOrderInit[i]][1]<Upper_temp[i])
				{
					coverage = false;
					break;
				}
			}
			if(coverage == true)
			{
				for(int s = 0 ; s<n; s++)
				{
					
					Lower[s] = fmin(Lower[s], Lower_temp[s]);
					Upper[s] = fmax(Upper[s], Upper_temp[s]);
					
				}
			}
		}
		
	}
	delete[] InqPosi;
}

void RescalTopHpo(int***& ResCIsMat, int***& ResCIsGridMat, IntegerVector& Lower, IntegerVector& Upper, const IntegerVector& EmpOrderInit, const int& n, const int& MM, const int& gridSize, const double& alpha)
{
	double* coveragePartition = new double[gridSize];
	IntegerVector Lower_temp(n), Upper_temp(n);
	bool coverage;
	double minCoverage;
	int alphaRescaled;

	for(int i = 0; i < n; i++)
	{
		Lower_temp[i] = 0;
		Upper_temp[i] = n-1;
	}
	
	
			coveragePartition[0] = fabs(PartitionCoverage(ResCIsMat, Lower_temp, Upper_temp, n, MM, 0)-(1-alpha));
			minCoverage = coveragePartition[0]; alphaRescaled = 0;
			for(int a = 1; a<gridSize; a++)
			{
				coveragePartition[a] = fabs(PartitionCoverage(ResCIsMat, Lower_temp, Upper_temp, n, MM, a)-(1-alpha));
				if(coveragePartition[a]<minCoverage)
				{
					alphaRescaled = a;
					minCoverage = coveragePartition[a];
				}
			}
			
			coverage = true;
			for(int i = 0; i<n; i++)
			{
				if(ResCIsGridMat[alphaRescaled][EmpOrderInit[i]][0]>Lower_temp[i] || ResCIsGridMat[alphaRescaled][EmpOrderInit[i]][1]<Upper_temp[i])
				{
					coverage = false;
					break;
				}
			}
			if(coverage == true)
			{
				for(int s = 0 ; s<n; s++)
				{
					
					Lower[s] = 0;
					Upper[s] = n-1;
					
				}
			}
		
		
	
	//delete[] InqPosi;
	delete[] coveragePartition;
}

// [[Rcpp::export]]
NumericMatrix PartitioningRankingLevelEqSigRescaled(NumericVector y, NumericVector sigma, NumericMatrix crit, NumericMatrix SampleCoverage, int MM, int n, int NbOfPermut, double alpha, int gridSize, bool trace)
{
	int powN2 = 1;
	for(int i = 0; i<n; i++)
	{
		powN2 *=2; 
	}
	//Rcout<<"powN2 = "<<powN2<<"\n";
	//Rcout<<"\n Start : gridSize = "<<gridSize<<", MM = "<<MM<<"\n";
	if(trace == true) Rcout<<"Initializing variables...\n";
	int*** ResCIsMat = new int**[gridSize]; // alpha values range in 0.05, ..., 0.4
	//std::vector<std::vector<std::vector<int> > > ResCIsMat; ResCIsMat.resize(gridSize);
	int*** ResCIsMat_temp = new int**[gridSize];
	//std::vector<std::vector<std::vector<int> > > ResCIsMat_temp; ResCIsMat_temp.resize(gridSize);
	int*** ResCIsGridMat = new int**[gridSize];
	//std::vector<std::vector<std::vector<int> > > ResCIsGridMat; ResCIsGridMat.resize(gridSize);
	int* AlphaRescaled = new int[powN2];
	//std::vector<int> AlphaRescaled; AlphaRescaled.resize(powN2);
	for(int i = 0; i<gridSize; i++)
	{
		ResCIsMat[i] = new int*[n];
		//ResCIsMat[i].resize(n);
		ResCIsMat_temp[i] = new int*[n];
		//ResCIsMat_temp[i].resize(n);
		ResCIsGridMat[i] = new int*[n];
		//ResCIsGridMat[i].resize(n);
		for(int j = 0; j<n; j++)
		{
			ResCIsMat[i][j] = new int[2*MM];
			//ResCIsMat[i][j].resize(2*MM);
			ResCIsMat_temp[i][j] = new int[2*MM];
			//ResCIsMat_temp[i][j].resize(2*MM);
			ResCIsGridMat[i][j] = new int[2];
			//ResCIsGridMat[i][j].resize(2);
		}
	}
	
	NumericMatrix ResCIs(n,2);
	if(trace == true) Rcout<<"Calculating SCIs for the ranks of a sample.\n";
	for(int a = 0; a<gridSize; a++)
	{
		for(int s = 0; s<MM; s++)
		{					
			// Caclulate simultaneous CIs for ranks using sample i
			ResCIs = PartitioningRankingLevelEqSig(SampleCoverage.column(s), sigma, crit.column(a), n, false);
			for(int i = 0; i<n; i++)
			{
				ResCIsMat[a][i][2*s] = ResCIs(i,0)-1; // Even for Lower CIs
				ResCIsMat[a][i][2*s+1] = ResCIs(i,1)-1; // Odd for Upper CIs	
			}
		}
		ResCIs = PartitioningRankingLevelEqSig(y, sigma, crit.column(a), n, false);
		for(int i = 0; i<n; i++)
		{
			ResCIsGridMat[a][i][0] = ResCIs(i,0)-1;
			ResCIsGridMat[a][i][1] = ResCIs(i,1)-1;
		}
		
	}
		
	IntegerVector RandPermutInit, RandPermut;
	RandPermutInit = seq(0,n-1);
	// Calculate the C_n^k matrix in order to use it in the UnrankCombin function
	unsigned long long int ** CnkMat = new unsigned long long int*[n];
	for (int i = 0; i<n; i++)
	{
		CnkMat[i] = new unsigned long long int[n];
		CnkMat[i][i] = 1;
		for (int j = 0; j<i; j++)
		{
			CnkMat[i][j] = binomialCoeff(i, j);
			CnkMat[j][i] = 0;
		}
	}

	// Definition and initialization of the vector of ranks
	IntegerVector Lower(n), Lower_temp(n), EmpOrderInit(n);
	IntegerVector Upper(n), Upper_temp(n);
	for (int i = 0; i<n; i++)
	{
		Lower[i] = i;
		Upper[i] = i;
		EmpOrderInit[i] = i;
	}
	// Calculate first the adjustment for the hypothesis mu_1=...=mu_n
	RescalTopHpo(ResCIsMat, ResCIsGridMat, Lower, Upper, EmpOrderInit, n, MM, gridSize, alpha);
	
	// Calculate CIs for ranks using the general proc without any permutations (correctly ordered hypotheses)
	if(trace == true) Rcout<<"Perform the correctly ordered hypotheses...\n";
	PartitioningRankingGeneralProcInit(ResCIsMat, ResCIsGridMat, AlphaRescaled, Lower, Upper, EmpOrderInit, CnkMat, crit, n, MM, gridSize, alpha, trace);
	bool CheckTrivialCIs = true;
	for(int i = 0; i<n; i++)
	{
		if(Lower[i] != 0 || Upper[i] != n-1)
		{
			CheckTrivialCIs = false;
			break;
		}
	}
	if(CheckTrivialCIs == false){
	if(trace == true) {Rcout<<"\n Forward permutations ("<<n-1<<" permutations)\n";}
	//### part1: permute (i,i+1,i+2,...)
	for(int I = 1; I<=n-1; I++)//#row index
	{
		if(trace == true) {Rcout<<I<<".";}
		for(int J = 1; J <= n-I; J++)//#column index
		{
			CheckTrivialCIs = true;
			for(int i = 0; i<n; i++)
			{
				if(Lower[i] != 0 || Upper[i] != n-1)
				{
					CheckTrivialCIs = false;
					break;
				}
			}
			if(CheckTrivialCIs == true) break;
			for(int s = 0; s<n; s++)
			{
				Lower_temp[s] = EmpOrderInit[s];
				Upper_temp[s] = EmpOrderInit[s];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][s][2*j];
						ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][s][2*j+1];
					}
				}			
			}
			
			//# Permute
			for(int s = J; s<=J+I-1; s++)
			{
				Lower_temp[s] = EmpOrderInit[s-1]; 
				Upper_temp[s] = EmpOrderInit[s-1];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][s-1][2*j];
						ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][s-1][2*j+1];
					}
				}			
			}
			Lower_temp[J-1] = EmpOrderInit[J+I-1];	
			Upper_temp[J-1] = EmpOrderInit[J+I-1];
			for(int i = 0; i<gridSize; i++)
			{
				for(int j = 0; j<MM; j++)
				{
					ResCIsMat_temp[i][J-1][2*j] = ResCIsMat[i][J+I-1][2*j];
					ResCIsMat_temp[i][J-1][2*j+1] = ResCIsMat[i][J+I-1][2*j+1];
				}
			}
		
			//# Calculate the SCI for ranks
			PartitioningRankingGeneralProc(ResCIsMat_temp, ResCIsGridMat, AlphaRescaled, Lower_temp, Upper_temp, clone(Lower_temp), CnkMat, crit, n, MM, gridSize, alpha, trace);
			
			//# Permute the ranks
			for(int s = J; s<=J+I-1; s++)
			{
				Lower[s-1] = fmin(Lower[s-1], Lower_temp[s]);
				Upper[s-1] = fmax(Upper[s-1], Upper_temp[s]);
			}
			Lower[J+I-1] = fmin(Lower[J+I-1], Lower_temp[J-1]);
			Upper[J+I-1] = fmax(Upper[J+I-1], Upper_temp[J-1]);
		}
	}
	
	if(trace == true) {Rcout<<"\n Applying backward permutations.\n";}
	int I;
	//### Part 2: permute (i+4,i+3,..)
	for(int J = 1; J<=n-1; J++)//#column index
	{
		if(trace == true) {Rcout<<J<<".";}
		I = 2;
		while(I<=(n-J))//#row index
		{
			CheckTrivialCIs = true;
			for(int i = 0; i<n; i++)
			{
				if(Lower[i] != 0 || Upper[i] != n-1)
				{
					CheckTrivialCIs = false;
					break;
				}
			}
			if(CheckTrivialCIs == true) break;

			for(int s = 0; s<n; s++)
			{
				Lower_temp[s] = EmpOrderInit[s];
				Upper_temp[s] = EmpOrderInit[s];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][s][2*j];
						ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][s][2*j+1];
					}
				}
			}
			//# Permute
			for(int s = J; s<= J+I-1; s++)
			{
			  	Lower_temp[s-1] = EmpOrderInit[s];
				Upper_temp[s-1] = EmpOrderInit[s];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s-1][2*j] = ResCIsMat[i][s][2*j];
						ResCIsMat_temp[i][s-1][2*j+1] = ResCIsMat[i][s][2*j+1];
					}
				}
			}
			Lower_temp[J+I-1] = EmpOrderInit[J-1];
			Upper_temp[J+I-1] = EmpOrderInit[J-1];
			for(int i = 0; i<gridSize; i++)
			{
				for(int j = 0; j<MM; j++)
				{
					ResCIsMat_temp[i][J+I-1][2*j] = ResCIsMat[i][J-1][2*j];
					ResCIsMat_temp[i][J+I-1][2*j+1] = ResCIsMat[i][J-1][2*j+1];
				}
			}
						
			PartitioningRankingGeneralProc(ResCIsMat_temp, ResCIsGridMat, AlphaRescaled, Lower_temp, Upper_temp, clone(Lower_temp), CnkMat, crit, n, MM, gridSize, alpha, trace);
			
			//# permute the ranks
			for(int s = J; s<=J+I-1; s++)
			{
			  Lower[s] = fmin(Lower[s], Lower_temp[s-1]);
			  Upper[s] = fmax(Upper[s], Upper_temp[s-1]);
			}
			Lower[J-1] = fmin(Lower[J-1], Lower_temp[J+I-1]);
			Upper[J-1] = fmax(Upper[J-1], Upper_temp[J+I-1]);
		 	I = I+1;
		}
	}
	
	if(trace == true) {Rcout<<"\n Applying "<<NbOfPermut<<" Random permutations.\n";}
	for(int I = 1; I<=NbOfPermut; I++)//#row index
	{
		CheckTrivialCIs = true;
		for(int i = 0; i<n; i++)
		{
			if(Lower[i] != 0 || Upper[i] != n-1)
			{
				CheckTrivialCIs = false;
				break;
			}
		}
		if(CheckTrivialCIs == true) break;

		if((trace == true) && (I % 500000 == 0)) {Rcout<<I<<".";}
		//# Permute
		RandPermut = sample(RandPermutInit,n);//clone(RandPermutInit);
		//std::random_shuffle(RandPermut.begin(), RandPermut.end());
		for(int s = 0; s<n; s++)
		{
			Lower_temp[s] = EmpOrderInit[RandPermut[s]];
			Upper_temp[s] = EmpOrderInit[RandPermut[s]];
			for(int i = 0; i<gridSize; i++)
			{
				for(int j = 0; j<MM; j++)
				{
					ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][RandPermut[s]][2*j];
					ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][RandPermut[s]][2*j+1];
				}
			}
		}
		PartitioningRankingGeneralProc(ResCIsMat_temp, ResCIsGridMat, AlphaRescaled, Lower_temp, Upper_temp, clone(Lower_temp), CnkMat, crit, n, MM, gridSize, alpha, trace);

		//# Permute the ranks
		for(int s = 0; s<n; s++)
		{
			Lower[RandPermut[s]] = fmin(Lower[RandPermut[s]], Lower_temp[s]);
			Upper[RandPermut[s]] = fmax(Upper[RandPermut[s]], Upper_temp[s]);
		}
	}
	}
	Rcout<<"\n";
	NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i]+1;
		CIs[n+i] = Upper[i]+1;
	}
	
	// Free some space.
	for (int i = 0; i<n; i++)
	{
		delete[] CnkMat[i];
	}
	delete[] CnkMat;
	for(int a = 0; a<gridSize; a++)
	{
		for(int i = 0; i<n; i++)
		{
			delete[] ResCIsMat[a][i];
			delete[] ResCIsMat_temp[a][i];
			delete[] ResCIsGridMat[a][i];
		}
		delete[] ResCIsMat[a];
		delete[] ResCIsMat_temp[a];
		delete[] ResCIsGridMat[a];
	}
	delete[] ResCIsMat;
	delete[] ResCIsMat_temp;
	delete[] ResCIsGridMat;
	delete[] AlphaRescaled;

	return CIs;
}


// [[Rcpp::export]]
NumericMatrix PartitioningRankingLevelUneqSig(NumericVector y, NumericVector sigma, NumericVector crit, int n, bool trace, const int NbOfPermut, bool SwapPerm)
{
	return 0;
}
// [[Rcpp::export]]
NumericMatrix PartitioningRankingBlockCorrectOrder(NumericVector y, NumericVector sigma, NumericVector crit, NumericVector MinBlock, NumericVector MaxBlock, IntegerVector Lower, IntegerVector Upper, int n, bool trace)
{
	return 0;
}
// An internal function for an instance of observations.
NumericMatrix OnlyBlockRanking_instance(const NumericVector& y, const NumericVector& sigma, const NumericVector& crit, IntegerVector& Lower, IntegerVector& Upper,  int n)
{
	return 0;
}
// ----------------------------------------------------------
// .......... To do Function ................................
// NumericMatrix OnlyBlockRanking(NumericVector y, NumericVector sigma, NumericVector crit, int n, bool trace, const int NbOfPermut)
// ---------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix OnlyBlockRanking(NumericVector y, NumericVector sigma, NumericVector crit, int n, bool trace, const int NbOfPermut, bool SwapPerm)
{
	return 0;
}
// ----------------------------------------------------------
// ...........................................
//................ Bracketing  ...................
// ...........................................
// ----------------------------------------------------------
/*double critFun(const int& x, const double& Slop)
{
	return Slop*x;
}*/
void WhichBounds(const NumericVector& y, const int& I, const int& J, int& minInd, int& maxInd)
{ // Determine the position of the min and max values
	minInd = 1; maxInd = J-I+1;
	int minY = y[I], maxY = y[J];
	for(int i = 1; i<J-I+1; i++)
	{
		if(y[i+I]<minY)
		{
			minY = y[i+I];
			minInd = i+1;
		}
		if(y[i+I]>maxY)
		{
			maxY = y[i+I];
			maxInd = i+1;
		}
	}
}

void IndividContribs(const NumericVector& y_temp, const NumericVector& sigma_temp, std::vector<std::vector<double> > & LogL, const int& K, const int& L, const double& Binf, const double& Bsup, std::vector<std::vector<double> >& IndividContribBlock, std::vector<std::vector<std::vector<double> > >& AverageBlock, const double& Slop, const double& Intercept, const int& n)
{
}

void ApproximatePartitionWrongOrderFast(const NumericVector& y_temp, const NumericVector& sigma_temp, const IntegerVector& EmpOrder_temp, std::vector<int> &Lower, std::vector<int>& Upper, std::vector<std::vector<double> > & LogL, std::vector<std::vector<std::vector<double> > > & AverageBlock, std::vector<std::vector<std::vector<double> > > & AverageBlock_temp, std::vector<std::vector<double> > & IndividContribBlock, std::vector<std::vector<double> > & IndividContribBlock_temp, const double& Slop, const double& Intercept, const int& n, const double& maxY, const double& minY)
{
}
// [[Rcpp::export]]
NumericMatrix ApproximatePartitionPermutations(NumericVector yInit, NumericVector sigmaInit, IntegerVector LowerInit, IntegerVector UpperInit, int n, double Slop, double Intercept, double minY, double maxY, bool trace, const bool SwapPerm, const int NbOfPermut)
{
	return 0;
}


// This function is the same as the one written in the R code, but I'm using it here only to make things run faster.
NumericMatrix TukeyProc(const NumericVector& y, const NumericVector& sigma, const double& qq, const int& n)
{
	return 0;
}
// [[Rcpp::export]]
NumericMatrix TukeyRankingLevelEqSigRescaled(NumericVector y, NumericVector sigma, NumericMatrix crit, NumericMatrix SampleCoverage, int MM, int n, int NbOfPermut, double alpha, int gridSize, bool trace)
{
	return 0;
}
void PartitioningRankingGeneralProcTuk(std::vector<std::vector<std::vector<int> > >& ResCIsMat, std::vector<std::vector<std::vector<int> > >& ResCIsGridMat, IntegerVector& Lower, IntegerVector& Upper, const IntegerVector& EmpOrderInit, unsigned long long int**& CnkMat, const NumericMatrix& crit, const int& n, const int& MM, const int& gridSize, const double& alpha, const bool& trace)
{
}
// [[Rcpp::export]]
NumericMatrix TukeyRankingLevelUneqSigRescaled(NumericVector y, NumericVector sigma, NumericMatrix crit, NumericMatrix SampleCoverage, int MM, int n, int NbOfPermut, double alpha, int gridSize, bool trace)
{
	return 0;
}


