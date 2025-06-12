//
//
#ifndef __UEBDAFUNCTIONS_HPP
#define __UEBDAFUNCTIONS_HPP

#include <ctime>
#include "Eigen/Dense"
#include<iomanip>
#include<iostream>
#include <fstream>
#include <cstring>
#include<sstream>
//#include <stdio.h>
#include "mpi.h"
#include <netcdf.h>
//#include <netcdf_par.h>
#include <cmath>
#include <set>
#include <vector>
#include <algorithm>
#include <random>

using namespace Eigen;

#include "eigenmvn.h"
//#include <cuda_runtime.h>
//#include "device_launch_parameters.h"
#pragma warning(disable : 4996)
//#include <algorithm>
//using namespace std;
//error handling for netcdf
#define  ERR(e) {std::cout<<"Error: "<< nc_strerror(e)<<std::endl; return 2; }
//error handling for cuda
//define  cuda_checkERR(err) { if (err != cudaSuccess) std::cout << "Error: "<<cudaGetErrorString(err)<< std::endl; exit(EXIT_FAILURE); }
//__device__ __host__
//void cuda_checkERR(cudaError_t err);

struct davar {
	char stateVar[256];   //state to assimilate
	float obsErrStdev;		  //observed (state or equivalent) var Standard deviation	
	float obsCorLength;           //correlation length of observation
	int  infType;            //type of obs state file (text 0 or nc 1)
	char stateInpFile[256];
	char infvarName[256];
	char inftimeVar[256];
	int numNcfiles;              // if infType = 1	
	int m_numObs;                //we need this only if infType = 0
	int nRecs;
};

class uebEnKF {  
	public:		
		uebEnKF(const char* daconFile, bool iuseHmatrix, std::vector<std::pair<int, int> > icellCoordinates, 
			int iYearDA, int iMonthDA, int iDayDA, double iHourDA, float iy0, float ix0) : useHmatrix(iuseHmatrix), cellCoordinates(icellCoordinates),
			YearDA(iYearDA), MonthDA(iMonthDA), DayDA(iDayDA), HourDA(iHourDA), y0(iy0), x0(ix0)
		{
			readDaContr(daconFile);
			//initDAMatrices();		
			debugOutputFile.open("debugOutput.txt", std::ios::out);
			debugOutputFile.close();
		}
		//uebEnKF();
		//uebEnKF(uebEnKF& uebEnKFCell0);
		//uebEnKF& operator= (uebEnKF& uebEnKFCell0);
		~uebEnKF(){
			//clean up
		}
		
		//da setttings
		int ng_gridSize, ns_statSize, mo_obseSize, es_enseSize;          //  matrix dimensions: grid dim, number of states, ensemble size, num. observations
		int nDaStates;             //number of assimiated obs
		float forcEnStdev;    //forcing ensemble standard deviation except temperature
		float tempEnStdev;    //temperature forcing ensemble standard deviation 
		float forcCorLength;     //correlation length for forcing
		static const int danumAssmnStates = 6;
		//float daStatesStdev[danumAssmnStates];  //states' standard deviation 

		float stateCorLength;     //correlation length
		float dxC;
		float dyC;
		std::vector<davar> daContArr;
		std::vector<std::pair<int, int> > cellCoordinates;
		Eigen::VectorXi Hc_hgVector;  //Hc= grid cell indices with observation
		Eigen::VectorXi Hc_hsVector;  //Hc= state index of observation
		std::vector<int> stateIndex;    //the state to assimilate 8.8.16 7th state for snow surface temp 
		//coordinate values of the bottom left of watershed grid file (from watershed.nc)
		float y0, x0;
		std::vector<int> startIndexDA, ncReadStartDA; // , tEndDA;										 
		double daTime;
		int YearDA, MonthDA, DayDA;
		double HourDA;
		char* uebdaStates[danumAssmnStates] = { "Us", "SWE", "tausn", "refDepth", "totalRefDepth", "TSURFs" };// "Wc", , "Tave"

		//debug file
		std::ofstream debugOutputFile;
		/*__host__ __device__*/
		void readDaContr(const char* daconFile);
		/*__host__ __device__*/
		void initDAMatrices();  // std::vector<std::pair<int, int> > cellCoordinates);
		void setHMatrices(float **daYcorrArr, float **daXcorrArr);
		void setHcMatrices(float **daYcorrArr, float **daXcorrArr);
		// mean of 1s and std.dev 
		//void getStdNorm_Samples_Default(int numSamples); // , double* &stdNormArr);
		/*__host__ __device__*/
		void getMultivarStdNorm_Samples_Forc_Default(float** &stdNormArr);
		/*__host__ __device__*/
		void getMultivarStdNorm_Samples_Forc_Tempr(float** &stdNormArr);
		/*__host__ __device__*/
		void getMultivarNorm_Samples_Obs_Default(float** &stdNormArr);
		
		/*__host__ __device__*/
		void getDaArr(int ida, int &numNc, float** RegArray, float &tcorVar, MPI::Intracomm inpComm, MPI::Info inpInfo);
		/*__host__ __device__*/
		void updateDaArr(float*** RegArray, float tcorVar, double** tvarArr);     //, int nDataPoints);
		/*__host__ __device__*/
		void getEnKFArrays(float** Xh_stateObsSpace, float** &X1_o_ensAnomalyObsSpace, float **&y_obsStateResidual, float **&Pzz_i_ObsStateCov_inv);
        	/*__host__ __device__*/
		void runEnKF(float** X1_o_ensAnomalyObsSpace, float **y_obsStateResidual, float **Pzz_i_ObsStateCov_inv, float** stateOutputArr, float** &stateOutputUpdate);
		void runEnKF(float** Xh_stateObsSpace, float** ensObservationErr, float** stateOutputArr, float** &stateOutputUpdate);       // , bool NormalDist);  //float* &ensembleUpdateArr,

		// flag to write warnings,...etc 
		static const int uebda_debugout1 = 1;  // print major debug information
		static const int uebda_debugout2 = 0;  // more detailed print 
		static const int uebda_debugout3 = 0;  // print in 'main' 
		int uebda_outflag;    //output print

	    /*__host__ __device__*/
		void  printDebugOutputs();
		void  readTStextFile_multiVal(const char* inforcFile, int ida, float* &yCorArr, float* &xCorArr, double* &tvarArr, float** &varvalArr);   // int &nrecords);
		void  readTStextFileTimeValPair(const char* inforcFile, std::vector<std::pair<double, float>> &tvar_in, int &nrecords);
		//time related funcs //snowxv.cpp							  
		// __host__ __device__  
		void  UPDATEtime(int &YEAR, int &MONTH, int &DAY, double &HOUR, double DT);
		//    function to return number of days in February checking for leap years
		//// __host__ __device__ 
		int lyear(int year);
		//To convert the real date to julian date
		//// __host__ __device__  
		int julian(int yy, int mm, int dd);
		//these were copied from functions.f90
		//COMPUTES JULIAN DATE, GIVEN CALENDAR DATE AND time.  INPUT CALENDAR DATE MUST BE GREGORIAN.  
		//// __host__ __device__ 
		double julian(int I, int M, int K, double H);
		//COMPUTES CALENDAR DATE AND time, GIVEN JULIAN DATE.  INPUT JULIAN DATE CAN BE BASED ON ANY UT-LIKE time SCALE
		//// __host__ __device__  
		void  calendardate(double TJD, int &I, int &M, int &K, double &H);

		//stdnormalSamplestxt << stdnormSamples.transpose() << std::endl;

    private:
        	
		Eigen::EigenMultivariateNormal<double> std_norm_dist_Forc_Default;
		Eigen::EigenMultivariateNormal<double> norm_dist_0Mean_Tempr;     //for temperature
		Eigen::EigenMultivariateNormal<double> norm_dist_0Mean_Default;		

		Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> Q_modelErrCov,
			H_hMaxtirx; //use this if not using function hX(); 
		bool useHmatrix;
		Eigen::VectorXf Z_obs;        
	    // covariance matrix based on distance between grid cells
		Matrix<double, Dynamic, Dynamic, RowMajor> covarF, R_obsErrCov;
		VectorXd meanF, meanT, meanM;

		//read multiple slubs (y,x arrays) along the time dim for data with t, y, x config---time as slowely varying array; and pvar_in already allocated
		// __host__ __device__ 
		int readNC_yxSlub_givenT(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int &tStart, float** &pvar_in, float &tcorvar, int &numNc, MPI::Intracomm inpComm, MPI::Info inpInfo);

};

#endif

