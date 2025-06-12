//#include "uebpg.h"
//#include "nctest.h"
//#include"mpi.h"
#include "uebpgdecls.h"
#include <time.h>
#include "uebdafunctions.h"
//#include <queue>
#pragma warning(disable : 4996)
using namespace std;

/*__global__ void callUEBRun(uebCell *uebCellArray, int nCells)
{
	int indx = blockIdx.x*blockDim.x + threadIdx.x;
	if (indx < nCells)
	{
		//float *dev_loOutArr = NULL;
		//cudaMalloc(&dev_loOutArr, 700000 * sizeof(float));		
		uebCellArray[indx].runUEB();
	}
}
__device__ __host__ void cuda_checkERR(cudaError_t err)
{
	if (err != cudaSuccess){
	   std::cout << "Error: " << cudaGetErrorString(err) << std::endl;
	   exit(EXIT_FAILURE);
    }
}

__host__ __device__ void checkDeviceMemory()
{
	//mem check on device
	size_t freeM, totalM;
	float freeMB, totalMB, allocMB;
	cudaMemGetInfo((size_t*)&freeM, (size_t*)&totalM);
	freeMB = (size_t)freeM / (1024*1024);
	totalMB = (size_t)totalM / (1024*1024);
	allocMB = totalMB - freeMB;
	printf(" %f  MB of   %f   MB total available device memory allocated. Remaining memory =   %f MB\n", allocMB, totalMB, freeMB);
}

__host__ __device__ void estimateThroughput(size_t dataSize, clock_t beginTime, clock_t endTime)
{
	double bandWidth = dataSize * 2.0;
	double GFLOPs = (double)(bandWidth * CLOCKS_PER_SEC) / (double)(endTime - beginTime);
	printf(" Estimated throughput =  %lf GFLOPs\n", GFLOPs);
}*/

int main(int argc, char* argv[])
{
	//gpu control	
	int threadsPerBlock = 255, blocksPerGrid = 1;
	//time
	float timeControl = 0.0, timeWS = 0.0, timeSitestate = 0.0, timeTSArrays = 0.0, timeParam = 0.0, timeParamSiteInptcontrol = 0.0, timeModelRun = 0.0;
	float* OutVarValues; //= new float*[70]; //[70];
	float*** aggoutvarArray = NULL;
	float ***ncoutArray = NULL;
		
	char conFile[256], paramFile1[256], sitevarFile1[256], inputconFile1[256], outputconFile1[256], watershedFile1[256], aggoutputconFile1[256], aggoutputFile1[256], daconFile1[256];
	char wsvarName1[256], wsycorName1[256], wsxcorName1[256];
	int **wsArray = NULL;
	int dimlen1 = 0, dimlen2 = 0;
	int wsfillVal = -9999;	
	float SiteState[32];
	float *wsxcorArray = NULL, *wsycorArray = NULL;
	int daAssimlate = 0;   //whether to run DAssimilation 1=yes; 0 (default)=no

	//params ParamVAlues;
	sitevar *strsvArray = new sitevar[32];
	char * svFile[32];
	char * svVarName[32];
	for (int i = 0; i < 32; i++){
		svFile[i] = new char[256];
		svVarName[i] = new char[256];
	}
	int svType[32];
	//inpforcvar strinpforcArray[13];
	//outputs
	pointOutput *pOut = NULL;
	ncOutput *ncOut = NULL;
	aggOutput *aggOut = NULL;
	ncOutput *daOut = NULL;
	int npout = 0, nncout = 0, naggout = 0, nZones = 0, ndaout=0;
	const char * zName = "Outletlocations"; //12.24.14 watershed zonning for aggregation--	
	float *z_ycor = NULL;
	float *z_xcor = NULL;
	int zoneid = 0;
	int *ZonesArr = NULL;
	//inptimeseries *strintsArray[11];
	float **RegArray[13] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };  //3.19.15  [5]; // [xstride]; //assuming max nc files for a variable =5
																											   //float *tcorvar[13], *tsvarArray[13],
	float *ycorArr = NULL, *xcorArr = NULL;
	int ncTotaltimestep[13] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	int numNc[13] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };				//6.24.14

	size_t ntimesteps = 0, nysteps = 0, nxsteps = 0;
	int tinitTime = 0;
	int npar = 32;
	/*int NUMtimeSTEP,NREFYR,NREFMO,NREFDAY,NumOP;*/
	int ModelStartDate[3], ModelEndDate[3]; //check this 
	double ModelStartHour, ModelEndHour, ModelDt, ModelUTCOffset;
	double modelSpan;	
	int inpDailyorSubdaily;
	//int numTotalTimeSteps;
	const int numOut = 71 ;
	char headerLine[256];
	int retvalue = 0;
	int numgrid = 0;

	const char* tNameout = "time";
	int outtSteps = 0;
	int outtStride = 1, outyStep = 1, outxStep = 1;
	float* t_out;
	float out_fillVal = -9999.0;
	int outDimord = 0, aggoutDimord = 1, daoutDimord = 1;
	int *yIndxArr = NULL, *xIndxArr = NULL;
	const char* tlong_name = "time";
	const char* tcalendar = "standard";
	char* uebVars[numOut] = { "Year", "Month", "Day", "dHour", "atff", "HRI", "Eacl", "Ema", "cosZen", "Ta", "P", "V", "RH", "Qsi", "Qli", "Qnet",
		"Us", "SWE", "tausn", "Pr", "Ps", "Alb", "QHs", "QEs", "Es", "SWIT", "QMs", "Q", "FM", "Tave", "TSURFs", "cump", "cumes",
		"cumMr", "NetRads", "smelt", "refDepth", "totalRefDepth", "cf", "Taufb", "Taufd", "Qsib", "Qsid", "Taub", "Taud",
		"Qsns", "Qsnc", "Qlns", "Qlnc", "Vz", "Rkinsc", "Rkinc", "Inmax", "intc", "ieff", "Ur", "Wc", "Tc", "Tac", "QHc",
		"QEc", "Ec", "Qpc", "Qmc", "Mc", "FMc", "SWIGM", "SWISM", "SWIR", "errMB", "Trange" };
	int outvarindx = 17, aggoutvarindx = 17;
	int size =1, rank = 0, irank = 0, jrank;
	double intermStart_Time = 0.0, startTimeT = 0.0, TotalTime = 0.0, paramSite_Time = 0.0, inputTS_Time = 0.0, computeRun_Time = 0.0, outputWrite_Time = 0.0, dataCopy_Time = 0.0;
	double TsReadTime = 0.0, TSStartTime, ComputeStartTime, ComputeTime = 0.0, OutWriteTime = 0.0;
	clock_t beginTime, endTime;
	//beginTime = clock();
	MPI::Init(argc, argv);
	//how many processes
	size = MPI::COMM_WORLD.Get_size(); //	MPI_Comm_size(MPI_COMM_WORLD,&size);
	//which rank is yours? 
	rank = MPI::COMM_WORLD.Get_rank(); //_Comm_rank(MPI_COMM_WORLD,&rank);
	//std::cout << "\n rank "<< rank << " of "<< size << " processes has started\n" << std::endl;	
	MPI::Intracomm worldComm = MPI::COMM_WORLD;
	MPI::Info worldInfo = MPI::INFO_NULL;
	if (rank == 0)
	{
		//microsecond wall time: to time block of work
		startTimeT = MPI::Wtime();
		intermStart_Time = MPI::Wtime();
		TsReadTime = 0.0;
		ComputeTime = 0.0;
	}
	//  Input Arguments		
	if (argc > 1)
	{
		//conFile = new char[sizeof(argv[0])];
		strcpy(conFile, argv[1]);
	}
	else
	{
		if (rank == 0)
			std::cout << "file not found exiting" << std::endl;
		MPI::Finalize();
		return 1;
		//cin >> conFile;
	}
	FILE* pconFile = fopen(conFile, "rt");
	fgets(headerLine, 256, pconFile);
	fscanf(pconFile, "%s\n %s\n %s\n %s\n %s\n %s\n %s\n", paramFile1, sitevarFile1, inputconFile1, outputconFile1, daconFile1, aggoutputFile1, watershedFile1);
	fscanf(pconFile, "%s %s %s\n", wsvarName1, wsycorName1, wsxcorName1);
	//new vs2012 appears to have issues with passing char[256] for const char*
	const char *paramFile = paramFile1, *sitevarFile = sitevarFile1, *inputconFile = inputconFile1, *outputconFile = outputconFile1, *aggoutputFile = aggoutputFile1,
		*watershedFile = watershedFile1, *wsvarName = wsvarName1, *wsycorName = wsycorName1, *wsxcorName = wsxcorName1, *daconFile = daconFile1;
	//read simulation related parameters including start and end datetimes, and model time step dt
	fscanf(pconFile, "%d %d %d %lf\n", &ModelStartDate[0], &ModelStartDate[1], &ModelStartDate[2], &ModelStartHour);
	fscanf(pconFile, "%d %d %d %lf\n", &ModelEndDate[0], &ModelEndDate[1], &ModelEndDate[2], &ModelEndHour);
	fscanf(pconFile, "%lf\n %lf\n %d\n %d %d %d\n %d %d %d\n %d\n %d\n", &ModelDt, &ModelUTCOffset, &inpDailyorSubdaily, &outtStride, &outyStep, &outxStep, &outDimord, &aggoutDimord, &daoutDimord, &threadsPerBlock, &daAssimlate);
	//close control file
	fclose(pconFile);
	//time units
	char tunits[256];
	int hhMod = (int)floor(ModelStartHour);
	int mmMod = (int)(remainder(ModelStartHour, 1.0) * 60);
	sprintf(tunits, "hours since %d-%d-%d %d:%d:00 UTC", ModelStartDate[0], ModelStartDate[1], ModelStartDate[2], hhMod, mmMod);
	const char* tUnitsout = tunits;
	//read watershed (model domain) netcdf file	
	retvalue = readwsncFile(watershedFile, wsvarName, wsycorName, wsxcorName, wsycorArray, wsxcorArray, wsArray, dimlen1, dimlen2, wsfillVal, worldComm, worldInfo);
	//std::cout<<"dim1 = "<<dimlen1<<" dim2 = "<< dimlen2<<std::endl;
	/*printf("fillvalue= %d ",wsfillVal);
	for(int i=0;i<dimlen1;i++){
	for(int j=0;j<dimlen2;j++)
	std::cout<<wsArray[i][j];
	std::cout<<"\n";
	}*/	
	//aggregation zone info
	float * wsArray1D = new float[dimlen1*dimlen2];
	for (int i = 0; i < dimlen1; i++)
		for (int j = 0; j < dimlen2; j++)
			wsArray1D[i*dimlen2 + j] = wsArray[i][j];
	//set contains unique id values
	std::set<int> zValues(wsArray1D, wsArray1D + (dimlen1*dimlen2));
	//std::cout << zValues.size() << std::endl;
	//std::remove_if(zValues.begin(), zValues.end(), [&wsfillVal](int a){ return a == wsfillVal; });
	std::set<int> fillSet;
    fillSet.insert (wsfillVal);
	//std::cout << "fill: " << fillSet.size() << " value: " << *(fillSet.begin())<<std::endl;
	std::vector<int> zVal(zValues.size());
	std::vector<int>::iterator it = std::set_difference(zValues.begin(), zValues.end(), fillSet.begin(), fillSet.end(), zVal.begin());  // exclude _FillValue
	zVal.resize(it - zVal.begin());  //now zVal contains unique watershed ids excluding fill value
	//std::cout << zVal.size()<<std::endl;
	z_ycor = new float[zVal.size()];
	z_xcor = new float[zVal.size()];
	//std::cout << zValues.size() << std::endl;
	nZones = zVal.size();
	for (int iz = 0; iz < zVal.size(); iz++)
	{
		//#_12.24.14 change these with actual outlet locations coordinates
		z_ycor[iz] = 0.0;
		z_xcor[iz] = 0.0;
		//std::cout << zValues[iz];		
	}
	//read site vars 
	//std::cout<<"Reading site variable ");
	readSiteVars(sitevarFile, strsvArray); //svDefaults,svFile,svVarName,svType);
	/*std::cout<<"\n site variables read \n");
	for(int i=0;i<32;i++)
	std::cout<<"%f ",strsvArray[i].svdefValue);
	std::cout<<"\n");*/
	for (int i = 0; i < 32; i++)
		if (strsvArray[i].svType == 1)
		{
			//std::cout<<"%d %s %s\n",i, strsvArray[i].svFile,strsvArray[i].svVarName);
			retvalue = read2DNC(strsvArray[i].svFile, strsvArray[i].svVarName, strsvArray[i].svArrayValues, worldComm, worldInfo);
			/*for(int ih=0;ih<13;ih++)
			{
				for(int jv=0;jv<16;jv++)
					std::cout<<"%f ",strsvArray[i].svArrayValues[ih][jv]);
				std::cout<<"\n");
			}*/
	    }
	paramSite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();	
	//vector of active cells
	std::vector<std::pair<int, int> > activeCells;
	for (int iy = 0; iy < dimlen1; iy++)
		for (int jx = 0; jx < dimlen2; jx++)
			if (wsArray[iy][jx] != wsfillVal && strsvArray[16].svType != 3)  //compute cell && no accumulation zone //***tbc what happens if it is accumulation zone?
				activeCells.push_back(std::make_pair(iy, jx));	
	int remLength = activeCells.size() % size;
	int extraCell = 0;
	if (remLength != 0) {
		if (rank < remLength)              //       (rank/remLength) == 0
			extraCell = 1;
	}
	int numCells0 = activeCells.size() / size;
	int numCells = numCells0 + extraCell;
	std::cout << "proc " << rank << " num cells = " << numCells << std::endl;
	if (rank == 0) {
		std::cout << "total number of uebcells = " << activeCells.size() << std::endl;
		std::cout << " size of uebCell class =  " << sizeof(uebCell) / ((double)(1024 * 1024)) << " MB" << std::endl;
	}
	// create ueb model gridcell instance and copy to arrays of grid cells
	uebCell objCell0 (paramFile, ModelStartDate, ModelEndDate, ModelStartHour, ModelEndHour, ModelDt, ModelUTCOffset, inpDailyorSubdaily,outtStride);	
	//output control
	readOutputControl(outputconFile, pOut, npout, ncOut, nncout, aggOut, naggout, daOut, ndaout);
	//create output netcdf
	outtSteps = objCell0.numTotalTimeSteps;             //--7.20.16 in the future conside saving output every outstrid'th t-step
	if (rank == 0)
		std::cout << "number of time steps: " << " " << objCell0.numTotalTimeSteps << std::endl;
	t_out = new float[outtSteps];
	for (int it = 0; it < outtSteps; ++it)
		t_out[it] = it*outtStride*ModelDt;      //in hours since model start time 
	//output array written to netcdf files
	ncoutArray = new float**[nncout];
	for (int inc = 0; inc < nncout; inc++)
	{
		ncoutArray[inc] = new float*[numCells];
		for (int nindx = 0; nindx < numCells; nindx++) {
			ncoutArray[inc][nindx] = new float[outtSteps];
			for (int to = 0; to < outtSteps; to++)
				ncoutArray[inc][nindx][to] = out_fillVal;
		}
	}
	yIndxArr = new int[numCells];
	xIndxArr = new int[numCells];
	//aggregated output arrays # There should be better way than this
	aggoutvarArray = new float**[nZones];
	float * totalAgg = new float[outtSteps];
	ZonesArr = new int[nZones];
	for (int j = 0; j < nZones; j++)
	{
		ZonesArr[j] = 0;
		aggoutvarArray[j] = new float*[naggout];
		for (int i = 0; i < naggout; i++)
		{
			aggoutvarArray[j][i] = new float[outtSteps];
			for (int it = 0; it < outtSteps; it++)
				aggoutvarArray[j][i][it] = 0.0;
		}
	}
	//netcdf output files
	for (int icout = 0; icout < nncout; icout++)
		retvalue = create3DNC_uebOutputs(ncOut[icout].outfName, (const char*)ncOut[icout].symbol, (const char*)ncOut[icout].units, tNameout, tUnitsout,
		tlong_name, tcalendar, outtSteps, outDimord, t_out, &out_fillVal, watershedFile, wsvarName, wsycorName, wsxcorName, worldComm, worldInfo);
	//create aggregate ouput file
	retvalue = create3DNC_uebAggregatedOutputs(aggoutputFile, aggOut, naggout, tNameout, tUnitsout, tlong_name, tcalendar, outtSteps, aggoutDimord, t_out, &out_fillVal,
		watershedFile, wsvarName, wsycorName, wsxcorName, nZones, zName, z_ycor, z_xcor, worldComm, worldInfo);
	outputWrite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();
	
	//read input /forcing control file--all possible entries of input control have to be provided
	objCell0.readInputForContr(inputconFile);	
	
	for (int it = 0; it < 13; it++) {
		if (objCell0.infrContArr[it].infType == 0) {
			RegArray[it] = new float *[1];
			readTStextFile(objCell0.infrContArr[it].infFile, RegArray[it][0], ncTotaltimestep[it]);   //tsvarArray[it][0] 3.19.15      ntimesteps[0] 12.18.14
		}
		else if (objCell0.infrContArr[it].infType == 1) {
			RegArray[it] = create2DArray_Contiguous(dimlen1, dimlen2);
		}
	}

	//if (daAssimlate == 1) {
	bool dauseHmatrix = true;
	uebEnKF objEnKF(daconFile, dauseHmatrix, activeCells, ModelStartDate[0], ModelStartDate[1], ModelStartDate[2], ModelStartHour, wsycorArray[0], wsxcorArray[0]);
	//std::vector<std::pair<double, float>> daTSArray;
	int* danumNc = new int [objEnKF.nDaStates];
	float ***daRegArray = new float** [objEnKF.nDaStates];	
	float **daYcorrArr = new float* [objEnKF.nDaStates];
	float **daXcorrArr = new float* [objEnKF.nDaStates];
	double **daTcorrArr = new double* [objEnKF.nDaStates];
	for (int ida = 0; ida < objEnKF.nDaStates; ida++) {
		if (objEnKF.daContArr[ida].infType == 0) {
			objEnKF.readTStextFile_multiVal((const char*)objEnKF.daContArr[ida].stateInpFile, ida, daYcorrArr[ida], daXcorrArr[ida], daTcorrArr[ida], daRegArray[ida]);
		}
		else if (objEnKF.daContArr[ida].infType == 1) {
				daRegArray[ida] = create2DArray_Contiguous(dimlen1, dimlen2);
		}
		danumNc[ida] = 0;
	}
	objEnKF.initDAMatrices();
	objEnKF.setHcMatrices(daYcorrArr, daXcorrArr);
	/*std::cout << "size of daTS Vector: " << daTSArray.size() << std::endl;
	std::cout << " Time values " << std::endl;
	for (int tvt = 0; tvt < daTSArray.size(); tvt++)
		std::cout << std::setprecision(15)<< daTSArray[tvt].first << " ";
	std::cout << std::endl;
	std::cout << " Observed values " << std::endl;
	for (int tvt = 0; tvt < daTSArray.size(); tvt++)
		std::cout << std::setprecision(15)<<daTSArray[tvt].second << " ";
	std::cout << std::endl;*/
	TsReadTime += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();
	std::cout << "proc " << rank << " after reading time series" << std::endl;
	//}
	// data assimilation outputs
	const int danumStates = objEnKF.danumAssmnStates;
	float ***daoutArray[danumStates];         // = NULL;
	for (int is = 0; is < danumStates; is++) {
		daoutArray[is] = new float**[numCells];
		for (int nindx = 0; nindx < numCells; nindx++) {
			daoutArray[is][nindx] = new float*[objCell0.numTotalTimeSteps];
			for (int it = 0; it < objCell0.numTotalTimeSteps; it++) {
				daoutArray[is][nindx][it] = new float[objEnKF.es_enseSize + 1];
				for (int ie = 0; ie < (objEnKF.es_enseSize + 1); ie++)
					daoutArray[is][nindx][it][ie] = out_fillVal;                        // no da, output will be fill val 8.5.16
			}
		}
	}
	//data assimilation netcdf output files
	for (int icout = 0; icout < ndaout; icout++)
		retvalue = createMultiDnc_uebOutputs(daOut[icout].outfName, (const char*)daOut[icout].symbol, (const char*)daOut[icout].units, tNameout, tUnitsout, tlong_name,
			tcalendar, outtSteps, daoutDimord, t_out, &out_fillVal, watershedFile, wsvarName, wsycorName, wsxcorName, "ensembleNumber", objEnKF.es_enseSize + 1, worldComm, worldInfo);
	
	bool NormalDist = true;
	float* ensForcingMultiplier[7];
	for (int iforc = 0; iforc < 7;iforc++)
		ensForcingMultiplier[iforc] = new float[objEnKF.es_enseSize];
	//float* ensObservationMultiplier = new float[objCell0.daContArr.nEns]; 
	//if (rank = 0) {
	float **multivarNormalDistSamplesForc[7];
	for (int iforc = 0; iforc < 7; iforc++)
		multivarNormalDistSamplesForc[iforc] = create2DArray_Contiguous(activeCells.size(), objEnKF.es_enseSize);
	float **stateOutputUpdateArr = create2DArray_Contiguous(danumStates, objEnKF.es_enseSize + 1);
	float **stateOutputArr = create2DArray_Contiguous(danumStates, objEnKF.es_enseSize);

	float **Xh_stateObsSpace = create2DArray_Contiguous(objEnKF.mo_obseSize, objEnKF.es_enseSize);
	float** ensObservationErr = create2DArray_Contiguous(objEnKF.mo_obseSize, objEnKF.es_enseSize);
	/*float **X1_o_ensAnomalyObsSpace = create2DArray_Contiguous(objEnKF.mo_obseSize, objEnKF.nEns);
	float **y_obsStateResidual = create2DArray_Contiguous(objEnKF.mo_obseSize, objEnKF.nEns);
	float **Pzz_i_ObsStateCov_inv = create2DArray_Contiguous(objEnKF.mo_obseSize, objEnKF.mo_obseSize);	//float ***stateOutputArr = new float **[activeCells.size()];																														
    *///MPI::Datatype stateEnsembleMatrix = MPI::FLOAT.Create_contiguous(danumStates * objEnKF.daContArr.nEns);
	//stateEnsembleMatrix.Commit();

	//}  //if (rank = 0)	
	uebCell *uebCellArr = new uebCell[numCells];
	std::cout << "proc " << rank << " created uebCell arrays" << std::endl;	  
	//for (irank = rank; irank < activeCells.size() - remLength; irank += size)
	int cellIndx = 0;
	for (irank = rank; irank < activeCells.size(); irank +=size)
	{
		//track grid cell		
		yIndxArr[cellIndx] = activeCells[irank].first;
		xIndxArr[cellIndx] = activeCells[irank].second;
		uebCellArr[cellIndx] = objCell0;
		uebCellArr[cellIndx].uebCellY = activeCells[irank].first;
		uebCellArr[cellIndx].uebCellX = activeCells[irank].second;
		for (int is = 0; is < 32; is++)
		{
			if (strsvArray[is].svType == 1)
				SiteState[is] = strsvArray[is].svArrayValues[uebCellArr[cellIndx].uebCellY][uebCellArr[cellIndx].uebCellX];
			else
				SiteState[is] = strsvArray[is].svdefValue;
		}
		uebCellArr[cellIndx].setSiteVars_and_Initconds(SiteState);
		//intialize ensemble states
		uebCellArr[cellIndx].setInitialEnsembleStates(objEnKF.es_enseSize);                   //, (const char*) objEnKF.daContArr.forcName);
		cellIndx++;
	}
	std::cout << "proc " <<rank<< " done setting site vars" << std::endl;
	//write outut headers
	for (irank = 0; irank < numCells; irank++)
	{
		//point outputs
		for (int ipout = 0; ipout < npout; ipout++)
		{
			if (uebCellArr[irank].uebCellY == pOut[ipout].ycoord && uebCellArr[irank].uebCellX == pOut[ipout].xcoord)
			{
				FILE* pointoutFile = fopen((const char*)pOut[ipout].outfName, "w");
				for (int vnum = 0; vnum < 4; vnum++)
					fprintf(pointoutFile, "%8s", uebVars[vnum]);            // header
				for (int vnum = 4; vnum < numOut; vnum++)
					fprintf(pointoutFile, "%16s ", uebVars[vnum]);            // header
				fclose(pointoutFile);
			}
		}
		//debug outputs
		if (irank % (outyStep*dimlen2 + outxStep) == 0) {
			char testPrint[256];
			char ind[256];
			strcpy(testPrint, "ZTest");
			sprintf(ind, "%d", uebCellArr[irank].uebCellY);
			strcat(testPrint, ind);
			strcat(testPrint, "_");
			sprintf(ind, "%d", uebCellArr[irank].uebCellX);
			strcat(testPrint, ind);
			strcat(testPrint, ".txt");
			FILE* outFile = fopen(testPrint, "w");
			for (int vnum = 0; vnum < 4; vnum++)
				fprintf(outFile, "%8s", uebVars[vnum]);            // header
			for (int vnum = 4; vnum < numOut; vnum++)
				fprintf(outFile, "%16s ", uebVars[vnum]);            // header
			fclose(outFile);
		}
	}

	paramSite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();
	//MPI::COMM_WORLD.Barrier();
	int *remRanks = NULL;  // new int[remLength];
	MPI::Group worldGroup = MPI::COMM_WORLD.Get_group();
	MPI::Group remGroup = NULL;        // worldGroup.Incl(remLength, remRanks);
	MPI::Intracomm remComm = NULL;     // MPI::COMM_WORLD.Create(remGroup);
	int newSize = -1;  //remComm.Get_size();
	int newRank = -1;  //remComm.Get_rank();	
	if (remLength > 0)  //if there are remaining compute cells after even distribution 
	{
		remRanks = new int[remLength];
		for (int ir = 0; ir < remLength; ir++)           // = leftBorder; ir < activeCells.size(); ir++)
			remRanks[ir] = ir;                           // [ir - leftBorder] = ir;
		remGroup = worldGroup.Incl(remLength, remRanks);
		remComm = MPI::COMM_WORLD.Create(remGroup);	
		if (extraCell == 1)                                    //(rank < remLength)         //only for processes in the new comm group
		{
			newSize = remComm.Get_size();
			newRank = remComm.Get_rank();
			std::cout << " rank " << rank << " of WorldComm has " << newRank << " of new comm group of size " << newSize << " remaining cells " << remLength << std::endl;
		} // if(rank < remLength)
	}

	//comment the following out when running without gpu
	/* 
	//start gpu code
	int curDv, curStr; 
	//memory on device checkDeviceMemory()
	//cuda err check
	cudaError_t err = cudaSuccess;		
	err = cudaGetDevice(&curDv);
	cuda_checkERR(err);
	std::cout <<"proc "<<rank<< " current device  = " << curDv << std::endl;
    cudaStream_t oStream;
    //cudaSetDevice(cudIndx);
	err = cudaStreamCreate(&oStream);
	cuda_checkERR(err);
	std::cout << "proc " << rank << " current stream  = " << oStream << std::endl;
	uebCell *dev_uebCellArr = NULL;
	err = cudaMalloc(&dev_uebCellArr, numCells*sizeof(uebCell));
	cuda_checkERR(err);
    std::cout<<"proc "<<rank<<" device memory alloc"<<std::endl;
	dataCopy_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();
	//end gpu code 
	*/
	double EJD = uebCellArr[0].julian(uebCellArr[0].modelEndDate[0], uebCellArr[0].modelEndDate[1], uebCellArr[0].modelEndDate[2], uebCellArr[0].modelEndHour);
	double currentModelDateTime = uebCellArr[0].julian(uebCellArr[0].modelStartDate[0], uebCellArr[0].modelStartDate[1], uebCellArr[0].modelStartDate[2], uebCellArr[0].modelStartHour);
	
	if (rank ==0)
		printf(" start date = %lf  end date = %lf\n", currentModelDateTime, EJD);
	//7.20.16 this may not be needed anymore
	/*if (uebCellArr[0].inpDailyorSubdaily == 0)                     //the last 24 steps of forcing are read at once; adjust the EJD so that last time datetime is EndDate - 23*DT
		EJD -= (22 * uebCellArr[0].modelDT);
	else*/
	EJD -= 0.125;
//TODO---what if the assimilation obs. ends before simulation time end? 
	bool updataDaArray = true;             //when a stored DA data is used, then get the next value
	float tcorVarVal = 0.0, dAtcorVarVal = 0.0;
	int tsimStep = 0;
	int iranko, ranko;
	std::cout << "grid cells with observations: " << std::endl;
	for (irank = 0; irank < objEnKF.mo_obseSize; irank++)
		std::wcout << objEnKF.Hc_hgVector[irank] << "   yc = " << uebCellArr[objEnKF.Hc_hgVector[irank]].uebCellY << "  xc " << uebCellArr[objEnKF.Hc_hgVector[irank]].uebCellX << std::endl;
	
	while (EJD > currentModelDateTime)
	{				
//        printf("proc %d  current date =  %lf\n",rank, currentModelDateTime);		
		objCell0.getInpForcArr(numNc, RegArray, tcorVarVal, worldComm, worldInfo);
//        std::cout<<"proc "<<rank<<" copied forcing arrays"<<std::endl;
		for (irank =  0; irank < numCells; irank++)
			uebCellArr[irank].updateInpForcArr(RegArray);   // , ncTotaltimestep);
//		std::cout << "proc " << rank << " forcing array updated." << std::endl;       // Number of time steps to run = " << uebCellArr[0].numSimTimeSteps << std::endl;
		//DA 7.15.16
		//if (daAssimlate == 1) {   
		if (updataDaArray) {
			for (int ida = 0; ida < objEnKF.nDaStates; ida++)
			   objEnKF.getDaArr(ida, danumNc[ida], daRegArray[ida], dAtcorVarVal, worldComm, worldInfo);
//			std::cout << "proc " << rank << " copied da state arrays" << std::endl;
			//for (irank = 0; irank < numCells; irank++)
			objEnKF.updateDaArr(daRegArray, dAtcorVarVal, daTcorrArr);
//			std::cout << "proc " << rank << " DA observed state array updated" << std::endl;
			updataDaArray = false;   // wait until this data is used in assimilation 
		}
		//}
		inputTS_Time += (MPI::Wtime() - intermStart_Time);
		intermStart_Time = MPI::Wtime();	
		//comment the following out when running without gpu
		/*
		//start gpu
		//cudaSetDevice(cudIndx);
		//memory on device checkDeviceMemory()
		//if(numNc==0){
		err = cudaMemcpyAsync(dev_uebCellArr, uebCellArr, numCells*sizeof(uebCell), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
		cuda_checkERR(err);
		err = cudaStreamSynchronize(oStream);
		cuda_checkERR(err);
		//}
		dataCopy_Time += (MPI::Wtime() - intermStart_Time);
		intermStart_Time = MPI::Wtime();
		std::cout << "process " << rank << " data copied to device" << std::endl; 	//std::cout<<err<<std::endl;
		//memory on device    checkDeviceMemory()
		// Launch Kernel	
		blocksPerGrid = (numCells + threadsPerBlock - 1) / threadsPerBlock;
		//call device run function	
		callUEBRun << < blocksPerGrid, threadsPerBlock, 0, oStream >> >(dev_uebCellArr, numCells);
		//synchronization
		err = cudaStreamSynchronize(oStream);   /// cudaDeviceSynchronize();
		cuda_checkERR(err);
		computeRun_Time += (MPI::Wtime() - intermStart_Time);
		intermStart_Time = MPI::Wtime();
		std::cout << "process " << rank << " finished device compute tasks" << std::endl;
		//copy data back
		err = cudaMemcpyAsync(uebCellArr, dev_uebCellArr, numCells*sizeof(uebCell), cudaMemcpyDeviceToHost, oStream);// != cudaSuccess)	
		cuda_checkERR(err);
		err = cudaStreamSynchronize(oStream);/// cudaDeviceSynchronize();
		cuda_checkERR(err);
		dataCopy_Time += (MPI::Wtime() - intermStart_Time);
		intermStart_Time = MPI::Wtime();
		std::cout << "process " << rank << " data copied to host" << std::endl;	//std::cout << err << std::endl;	
		//end of gpu call
		*/
		/*for (irank = rank; irank < activeCells.size() - remLength; irank += size)
		     for (int it = 0; it < outtSteps; ++it)
			     std::cout << " " << OutVarValues[70 *it + 17];
		*/

		//control simulation ---regardless of whether filter run was called or not
		//this sets all variables not passed through filter
		//this avoids running the ensemble when there is no DA
		for (irank = 0; irank < numCells; irank++)
		{
			uebCellArr[irank].setForcingAndRadiationParamterization();
			//uebCellArr[irank].runUEB(objEnKF.es_enseSize);           //
			uebCellArr[irank].runUEB();           //no ensemble run, no da
		}
		// call Ensemble Kalman Filter
		if (rank == 0) {
			std::cout << std::setprecision(15) << std::endl << " current model datetime = " << currentModelDateTime << ";  data assimilation next time = " << objEnKF.daTime << "; 0.5DT = " << 0.5*ModelDt / 24;
			std::cout << std::setprecision(15) << "; CurrModelDateTime - NextDaDateTime = " << fabs(objEnKF.daTime - currentModelDateTime) << std::endl;
		}
		if (daAssimlate == 1)
		{
			if (fabs(objEnKF.daTime - currentModelDateTime) <= 0.5*ModelDt / 24)    //call only when there is assimilation data---else no data (o) will be printed 8.5.16  ---<= 0.5 DT assimilate to the closest time step
			{
				//generate multi-variate forcing and observation multipliers
				if (rank == 0)
				{
					objEnKF.getMultivarStdNorm_Samples_Forc_Tempr(multivarNormalDistSamplesForc[0]);
					for (int iforc = 1; iforc < 7; iforc++)
						objEnKF.getMultivarStdNorm_Samples_Forc_Default(multivarNormalDistSamplesForc[iforc]);
					//
					objEnKF.getMultivarNorm_Samples_Obs_Default(ensObservationErr);
				}
				//get states in observation space
				for (irank = 0; irank < objEnKF.mo_obseSize; irank++)
				{
					iranko = objEnKF.Hc_hgVector[irank] / size;
					ranko = objEnKF.Hc_hgVector[irank] % size;

					for (int iforc = 0; iforc < 7; iforc++)
					{
						if (rank == 0)
							MPI::COMM_WORLD.Isend(&multivarNormalDistSamplesForc[iforc][objEnKF.Hc_hgVector[irank]][0], objEnKF.es_enseSize, MPI::FLOAT, ranko, 1000+100*irank+10*ranko+iforc );
						if (rank == ranko)
							MPI::COMM_WORLD.Recv(&ensForcingMultiplier[iforc][0], objEnKF.es_enseSize, MPI::FLOAT, 0, 1000 + 100 * irank + 10 * ranko + iforc);
					}
					if (rank == ranko)
						uebCellArr[iranko].runUEBEnsembles(objEnKF.es_enseSize, ensForcingMultiplier, objEnKF.Hc_hsVector[irank], Xh_stateObsSpace[irank]);
					MPI::COMM_WORLD.Bcast(&Xh_stateObsSpace[irank][0], objEnKF.es_enseSize, MPI::FLOAT, ranko);
				}
				MPI::COMM_WORLD.Bcast(&ensObservationErr[0][0], objEnKF.mo_obseSize * objEnKF.es_enseSize, MPI::FLOAT, 0);
				if (objEnKF.uebda_debugout3 == 1)
				{
					std::cout << " Proc: " << rank << " state in obs. space: " << std::endl;
					for (irank = 0; irank < objEnKF.mo_obseSize; irank++) {
						for (int ie = 0; ie < objEnKF.es_enseSize; ie++)
							std::cout << Xh_stateObsSpace[irank][ie] << "  ";
						std::wcout << endl;
					}
					std::cout << " Proc: " << rank << " obs error: " << std::endl;
					for (irank = 0; irank < objEnKF.mo_obseSize; irank++) {
						for (int ie = 0; ie < objEnKF.es_enseSize; ie++)
							std::cout << ensObservationErr[irank][ie] << "  ";
						std::wcout << endl;
					}
				}
				//objEnKF.getEnKFArrays(Xh_stateObsSpace, X1_o_ensAnomalyObsSpace, y_obsStateResidual, Pzz_i_ObsStateCov_inv);
				for (irank = 0; irank < numCells0; irank++)
				{
					for (int iforc = 0; iforc < 7; iforc++)
						MPI::COMM_WORLD.Scatter(&multivarNormalDistSamplesForc[iforc][irank*size][0], objEnKF.es_enseSize, MPI::FLOAT, ensForcingMultiplier[iforc], objEnKF.es_enseSize, MPI::FLOAT, 0);
					//run ensemlbes			
					uebCellArr[irank].runUEBEnsembles(objEnKF.es_enseSize, ensForcingMultiplier, stateOutputArr);
					//if (irank == 0) {
						/*std::cout << std::endl << " cell 0 SWE ensemble outputs for time step: " << tsimStep << std::endl;
						for (int ie = 0; ie < objEnKF.daContArr.nEns; ie++)
							std::cout << stateOutputArr[1][ie] << "   ";
						std::cout << std::endl;*/
						//}			
					objEnKF.runEnKF(Xh_stateObsSpace, ensObservationErr, stateOutputArr, stateOutputUpdateArr);
					if (objEnKF.uebda_debugout3 == 1)
					{
						std::cout << std::endl << " sample updated state ensemble for time step: " << tsimStep <<" grid cell: "<<irank<< std::endl;
						for (int is = 0; is < danumStates; is++)
						{
							for (int ie = 0; ie < objEnKF.es_enseSize + 1; ie++)
								std::cout << stateOutputUpdateArr[is][ie] << " ";
							std::cout << std::endl;
						}
						std::cout << std::endl << " state ensemble mean for time step: " << tsimStep << " grid cell: " << irank << std::endl;
						for (int is = 0; is < danumStates; is++)
							std::cout << stateOutputUpdateArr[is][objEnKF.es_enseSize] << " ";            //is = danumStates*irank*size + is = danumStates*0*size + is
						std::cout << std::endl;
					}
					//update background state for next step		
					uebCellArr[irank].updateBackgroundStates(objEnKF.es_enseSize, stateOutputUpdateArr);
					for (int is = 0; is < danumStates; is++)						
						for (int ie = 0; ie < objEnKF.es_enseSize + 1; ie++)
							daoutArray[is][irank][tsimStep][ie] = stateOutputUpdateArr[is][ie];					
				}
				//MPI::COMM_WORLD.Barrier();	
				if (remLength > 0)  //if there are remaining compute cells after even distribution 
				{
					if (extraCell == 1)                                    //(rank < remLength)         //only for processes in the new comm group
					{
						for (int iforc = 0; iforc < 7; iforc++)
							remComm.Scatter(&multivarNormalDistSamplesForc[iforc][numCells0*size][0], objEnKF.es_enseSize, MPI::FLOAT, ensForcingMultiplier[iforc], objEnKF.es_enseSize, MPI::FLOAT, 0);
						//run ensemlbes			
						uebCellArr[numCells0].runUEBEnsembles(objEnKF.es_enseSize, ensForcingMultiplier, stateOutputArr);      // , ensAnomaly, ensMean[numCells0], NormalDist);
						/*if (irank == 0) {
						std::cout << std::endl << " cell 0 SWE ensemble outputs for time step: " << tsimStep << std::endl;
						for (int ie = 0; ie < objEnKF.daContArr.nEns; ie++)
						std::cout << stateOutputArr[1][ie] << "   ";
						std::cout << std::endl;
						}*/
						objEnKF.runEnKF(Xh_stateObsSpace, ensObservationErr, stateOutputArr, stateOutputUpdateArr);
						if (objEnKF.uebda_debugout3 == 1)
						{
							std::cout << std::endl << " sample updated state ensemble for time step: " << tsimStep << " grid cell: " << irank << std::endl;
							for (int is = 0; is < danumStates; is++)
							{
								for (int ie = 0; ie < objEnKF.es_enseSize + 1; ie++)
									std::cout << stateOutputUpdateArr[is][ie] << " ";
								std::cout << std::endl;
							}
							std::cout << std::endl << " state ensemble mean for time step: " << tsimStep << " grid cell: " << irank << std::endl;
							for (int is = 0; is < danumStates; is++)
								std::cout << stateOutputUpdateArr[is][objEnKF.es_enseSize] << " ";            //is = danumStates*irank*size + is = danumStates*0*size + is
							std::cout << std::endl;
						}
						//update background state for next step	
						uebCellArr[numCells0].updateBackgroundStates(objEnKF.es_enseSize, stateOutputUpdateArr);
						for (int is = 0; is < danumStates; is++)								
							for (int ie = 0; ie < objEnKF.es_enseSize + 1; ie++)
								daoutArray[is][numCells0][tsimStep][ie] = stateOutputUpdateArr[is][ie];						
					} // if(rank < remLength)
				}
				updataDaArray = true;            //copy the next observation 
			}
			else     // no da, make sure intial states are set for next time step
			{
				for (irank = 0; irank < numCells; irank++)
					uebCellArr[irank].setNextStepStates(objEnKF.es_enseSize);
			}
		}
    	//if(tsimStep == 394)
		//	std::getchar();
		computeRun_Time += (MPI::Wtime() - intermStart_Time);
		intermStart_Time = MPI::Wtime();
		//
		for (irank = 0; irank < numCells; irank++)
		{
			//write nc outputs				
			for (int icout = 0; icout < nncout; icout++)
			{
				for (int vindx = 0; vindx < 70; vindx++)     //TODO: try to do without looping?  
				{
					if (strcmp(ncOut[icout].symbol, uebVars[vindx]) == 0)
					{
						outvarindx = vindx;
						break;
					}
				}
				ncoutArray[icout][irank][tsimStep] = uebCellArr[irank].OutVarValues[outvarindx];        //t_out[it]3.20.15  //use timeStiride to sample outputs if it is dense (e.g hourly data for a year may be too big to save in one nc file)
				//write var values
				//retvalue = WriteTSto3DNC((const char*)ncOut[icout].outfName, (const char*)ncOut[icout].symbol, outDimord, uebCellY, uebCellX, outtSteps, t_out);                //, worldComm, worldInfo);
			}
			//#_??aggregated outputs 12.24.14
			zoneid = wsArray[uebCellArr[irank].uebCellY][uebCellArr[irank].uebCellX] - 1;
			//ZonesArr[zoneid] += 1;
			for (int iagout = 0; iagout < naggout; iagout++)
			{
				for (int vindx = 0; vindx < 70; vindx++)
				{
					if (strcmp(aggOut[iagout].symbol, uebVars[vindx]) == 0)
					{
						aggoutvarindx = vindx;
						break;
					}
				}
				aggoutvarArray[zoneid][iagout][tsimStep] += uebCellArr[irank].OutVarValues[aggoutvarindx];
			}
			//point outputs
			for (int ipout = 0; ipout < npout; ipout++)
			{
				if (uebCellArr[irank].uebCellY == pOut[ipout].ycoord && uebCellArr[irank].uebCellX == pOut[ipout].xcoord)
					uebCellArr[irank].printPointOutputs((const char*)pOut[ipout].outfName);
			}
			//debug outputs
			if (uebCellArr[irank].snowdgt_outflag == 1)
				uebCellArr[irank].printDebugOutputs();
		}
		outputWrite_Time += (MPI::Wtime() - intermStart_Time);
		intermStart_Time = MPI::Wtime();

		//if (tsimStep == 405) std::getchar();		

		//std::cout<<"proc "<<rank<<" before currDT compute"<<std::endl;
		//udate time   7.19.16 ---look for alternative later
		tsimStep++;
		for (irank = 0; irank < numCells; irank++)
			uebCellArr[irank].updateSimTime();
		currentModelDateTime = uebCellArr[0].julian(uebCellArr[0].modelStartDate[0], uebCellArr[0].modelStartDate[1], uebCellArr[0].modelStartDate[2], uebCellArr[0].modelStartHour);
		//progress is calculated and written here
		//numgrid + uebCellArr[0].numSimTimeSteps;
		if (rank == 0) {                               /// && numgrid % (dimlen1*uebCellArr[0].numSimTimeSteps) == 0)
			std::cout<<std::endl << " time step: " << tsimStep << std::endl;
			std::cout << "   percent completed: " << ((float)tsimStep / objCell0.numTotalTimeSteps)*100.0 << " %" << std::endl;
		}
		std::fflush(stdout);		
		//std::cout<<"proc "<<rank<<" at the end of time loop"<<std::endl;
	}
	//std::cout << "number of active cells = " << activeCells.size() << " number cellIndx = " << cellIndx << std::endl;
	std::cout << "process " << rank << " completed computation" << std::endl;
	MPI::COMM_WORLD.Barrier();
	//cellIndx++;
	//nc outputs
	for (int icout = 0; icout < nncout; icout++)
		retvalue = WriteTSto3DNC_Block((const char*)ncOut[icout].outfName, (const char*)ncOut[icout].symbol, outDimord, yIndxArr, xIndxArr, numCells-extraCell, outtSteps, ncoutArray[icout], worldComm, worldInfo);              //, worldComm, worldInfo);
	//da output
	for (int icout = 0; icout < ndaout; icout++)
	{    
		for (int sindx = 0; sindx < danumStates; sindx++)     //TODO: try to do without looping?  
		{
			if (strcmp(daOut[icout].symbol,objEnKF.uebdaStates[sindx]) == 0)
			{
				outvarindx = sindx;
				break;
			}
		}
		retvalue = WriteTStoMultiDnc_Block((const char*)daOut[icout].outfName, (const char*)daOut[icout].symbol, daoutDimord, yIndxArr, xIndxArr, numCells0, objCell0.numTotalTimeSteps, objEnKF.es_enseSize + 1, daoutArray[outvarindx], worldComm, worldInfo);              //, worldComm, worldInfo);
	}
	std::cout << "process " << rank << " wrote block outputs" << std::endl;
	//MPI::COMM_WORLD.Barrier();	
	if (remLength > 0)  //if there are remaining compute cells after even distribution 
	{
		if (extraCell == 1){                                    //(rank < remLength)         //only for processes in the new comm group
			//write var values
			for (int icout = 0; icout < nncout; icout++)
				retvalue = WriteTSto3DNC((const char*)ncOut[icout].outfName, (const char*)ncOut[icout].symbol, outDimord, yIndxArr[numCells - 1], xIndxArr[numCells - 1], outtSteps, ncoutArray[icout][numCells - 1], remComm, worldInfo);
			//da output
			for (int icout = 0; icout < ndaout; icout++) {  //
				for (int sindx = 0; sindx < danumStates; sindx++){     //TODO: try to do without looping?  
					if (strcmp(daOut[icout].symbol,objEnKF.uebdaStates[sindx]) == 0){
						outvarindx = sindx;
						break;
					}
				}
				retvalue = WriteTStoMultiDnc_Block((const char*)daOut[icout].outfName, (const char*)daOut[icout].symbol, daoutDimord, &yIndxArr[numCells0], &xIndxArr[numCells0], 1, objCell0.numTotalTimeSteps, objEnKF.es_enseSize + 1, &daoutArray[outvarindx][numCells0], remComm, worldInfo);              //, worldComm, worldInfo);
			}
			std::cout << "process " << rank << " wrote single outputs" << std::endl;
		} 
	}     // if(rank < remLength)		
	outputWrite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();	
	//MPI::COMM_WORLD.Barrier();
	//aggregation/ reduction 
	for (int it = 0; it < outtSteps; it++)
		totalAgg[it] = 0.0;
	//
	for (irank = 0; irank < numCells; irank++)
	{
		//#_??aggregated outputs 12.24.14
		zoneid = wsArray[uebCellArr[irank].uebCellY][uebCellArr[irank].uebCellX] - 1;
		ZonesArr[zoneid] += 1;
	}
	//MPI::COMM_WORLD.Barrier();
	int rankrec = 0;               //receiver rank 
	int totalZonecells = 1, zonValue = 0;
	for (int izone = 0; izone < nZones; izone++)
	{
		rankrec =  izone*size / nZones;
		//std::cout << "process " << rank << " before first reduce to rank: " << rankrec << std::endl;
		zonValue = ZonesArr[izone];
		MPI::COMM_WORLD.Reduce(&zonValue, &totalZonecells, 1, MPI::INT, MPI::SUM, rankrec);
		//std::cout<<"process "<<rank<<" total zone cells "<<totalZonecells<<std::endl;
		if (totalZonecells < 1)
			totalZonecells = 1;		
		for (int iagout = 0; iagout < naggout; iagout++)
		{
			//std::cout << "process " << rank << " before reduce of output " << iagout << std::endl;
			//if (rank == rankrec) MPI::COMM_WORLD.Reduce(MPI::IN_PLACE, aggoutvarArray[izone][iagout],outtSteps, MPI::FLOAT, MPI::SUM, rankrec); else 
			MPI::COMM_WORLD.Reduce(aggoutvarArray[izone][iagout], totalAgg, outtSteps, MPI::FLOAT, MPI::SUM, rankrec);
			//std::cout << "process " << rank << " waiting for writing" << std::endl;
			//#_12.28.14 aggregation operation needs defining
			if (rank == rankrec)
			{
				if (strcmp(aggOut[iagout].aggop, "AVE") == 0)
				{
					for (int it = 0; it < outtSteps; it++)
						totalAgg[it] =  totalAgg[it] / totalZonecells;				//aggoutvarArray[izone][iagout][it] / totalZonecells;	//
				}
				/*else
				{
					for (int it = 0; it < outtSteps; it++)
						totalAgg[it] = aggoutvarArray[izone][iagout][it];	// totalAgg[it] / totalZonecells;
				}*/
				//std::cout << "process " << rank << " before write of output " << iagout << " for zone: " << izone << std::endl;
				retvalue = Write_uebaggTS_toNC(aggoutputFile, aggOut[iagout].symbol, aggoutDimord, izone, outtSteps, totalAgg);
				//std::cout << "process: " << rank << " done writing output: " << iagout << " for zone " << izone << std::endl;
			}
		}
	}
	outputWrite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();
	std::cout << "process " << rank << " wrote aggregated outputs" << std::endl;
	//MPI::COMM_WORLD.Barrier();
	//deallocate memory ====#_*_#______Needs revisiting; some of the arrays are not deleted 6.23.13	
	//clear stream 
	/*err = cudaStreamDestroy(oStream);
	cuda_checkERR(err);	
	//free device memory
	err = cudaFree(dev_uebCellArr); // != cudaSuccess)	
	cuda_checkERR(err);
	std::cout << "process " << rank << " device memory freed" << std::endl;*/
	// Free host memory
	delete[] uebCellArr;
	std::cout << "process " << rank << " uebCell arrary freed" << std::endl;
	for (int i = 0; i < dimlen1; i++)
		delete[] wsArray[i];
	delete[] wsArray;
	std::cout<<"proc "<<rank<<" freed watershed arrays"<<std::endl;
	for (int i = 0; i < 32; i++)
	{
		if (strsvArray[i].svType == 1)
		{
			for (int j = 0; j < dimlen1; j++)
				delete[] strsvArray[i].svArrayValues[j];
			delete[] strsvArray[i].svArrayValues;
		}
	}
	delete[] strsvArray;
	std::cout<<"proc "<<rank<<" freed site vars arrays"<<std::endl;
	for (int it = 0; it < 13; it++) {
		if (objCell0.infrContArr[it].infType == 0) {
			delete[] RegArray[it][0];
			delete[] RegArray[it];
			RegArray[it][0] = NULL;
			RegArray[it] = NULL;			
		}
		else if (objCell0.infrContArr[it].infType == 1) {
			delete2DArray_Contiguous(RegArray[it]);
			RegArray[it] = NULL;
		}
	}
	std::cout << "proc " << rank << " freed forcing arrays" << std::endl;
	//delete[] tsvarArray[kx];
	/*for (int it = 0; it < 13; it++)       //10-->12   6.26.14
	{
		delete[] tsvarArray[it];
	}*/
	//delete[] tsvarArray;
	/*if (rank < remLength)
	{
		for (int it = 0; it < 13; it++)
		{
			if (strinpforcArray[it].infType == 1)
				delete[] tcorvar[it];
		}
		//delete[] tcorvar;
	}*/
    /*for (int it = 0; it <13; it++)       //  6.26.14
	{
		if (tsvarArrayTemp[it] != NULL)
		delete3DArrayblock_Contiguous(tsvarArrayTemp[it]);
	}
	std::cout<<"process "<<rank<<"freeed tsvartemp"<<std::endl;
    */
//delete RegArray
	//
	for (int inc = 0; inc < nncout; inc++)
	{
		for (int nindx = 0; nindx < numCells; nindx++)
			delete[] ncoutArray[inc][nindx];
		delete[] ncoutArray[inc];
	}
	delete[] ncoutArray;
	std::cout << "Process " << rank << " freed ncoutArray" << std::endl;
	for (int zk = 0; zk < nZones; zk++)
	{
		for (int ig = 0; ig < naggout; ig++)
			delete[] aggoutvarArray[zk][ig];
		delete[] aggoutvarArray[zk];
	}
	delete[] aggoutvarArray;
	std::cout<<"process "<<rank<<" freed aggregated output array"<<std::endl;
	/*for(int k=0 ;k<numOut; k++)
		delete[] OutVarValues[k];
	delete []OutVarValues; */
	if (remLength > 0)  //if there are remaining compute cells after even distribution 
		delete[] remRanks;
	remRanks = NULL;
	std::cout << "process " << rank << " freed extra rank arrays" << std::endl;
	// data assimilation related
	//if (daAssimlate == 1){	
	for (int ida = 0; ida < objEnKF.nDaStates; ida++) 
	{
		if (objEnKF.daContArr[ida].infType == 0)
		{
			for (int ir = 0; ir < objEnKF.daContArr[ida].nRecs; ir++) {
				delete[] daRegArray[ida][ir];
				daRegArray[ida][ir] = NULL;
			}
			delete[] daRegArray[ida];
			daRegArray[ida] = NULL;  //DA 7.15.16
			delete[] daXcorrArr[ida];
			delete[] daYcorrArr[ida];
			delete[] daTcorrArr[ida];
			daXcorrArr[ida] = NULL;
			daYcorrArr[ida] = NULL;
			daTcorrArr[ida] = NULL;
		}
		else if (objEnKF.daContArr[ida].infType == 1) 
		{
			delete2DArray_Contiguous(daRegArray[ida]);
			daRegArray[ida] = NULL;
		}
	}
	delete[] daRegArray;
	daRegArray = NULL;  //DA 7.15.16
	delete[] daXcorrArr;
	delete[] daYcorrArr;
	delete[] daTcorrArr;
	daXcorrArr = NULL;
	daYcorrArr = NULL;
	daTcorrArr = NULL;

	std::cout << "proc " << rank << " freed da obs arrays" << std::endl;
	for (int is = 0; is < danumStates; is++) {
		for (int nindx = 0; nindx < numCells; nindx++) {
			for (int it = 0; it < objCell0.numTotalTimeSteps; it++)
				delete[] daoutArray[is][nindx][it];
			delete[] daoutArray[is][nindx];
		}
		delete[] daoutArray[is];
		daoutArray[is] = NULL;
	}
	std::cout << "process " << rank << " freed da output array" << std::endl;
	delete2DArray_Contiguous(Xh_stateObsSpace);
	Xh_stateObsSpace = NULL;
	delete2DArray_Contiguous(ensObservationErr);
	ensObservationErr = NULL;
	std::cout << "process " << rank << " da state arrays" << std::endl;
	for (int iforc = 0; iforc < 7; iforc++) {
		delete[] ensForcingMultiplier[iforc];
		ensForcingMultiplier[iforc] = NULL;
	}
	//if (rank = 0) {
	for (int iforc = 0; iforc < 7; iforc++) {
		delete2DArray_Contiguous(multivarNormalDistSamplesForc[iforc]);
		multivarNormalDistSamplesForc[iforc] = NULL;
	}
	delete2DArray_Contiguous(stateOutputUpdateArr);
	stateOutputUpdateArr = NULL;
	delete2DArray_Contiguous(stateOutputArr);
	stateOutputArr = NULL;
	std::cout << "process " << rank << " freed ensemble mean, anomaly, and state aggregation arrays" << std::endl;
	//release mpi datatype contiguous
	//stateEnsembleMatrix.Free();
	//std::cout << "process " << rank << " mpi derived data type freed" << std::endl;
		//}
	//}
	std::cout << "Process " << rank << " finished" << std::endl;
	fflush(stdout);
	//MPI::COMM_WORLD.Barrier();
	if (rank == 0)
	{  
		//endTime = clock();
		TotalTime = MPI::Wtime() - startTimeT;                   //
		std::cout << "Time in seconds" << std::endl;
		std::cout << "Reading param  site state input control:  " << paramSite_Time << std::endl;
		std::cout << "Reading input TS txt arrays:  " << TsReadTime << std::endl;
		std::cout << "Reading input total TS arrays:  " << inputTS_Time << std::endl;
		std::cout << "Model simulation run time:  " << computeRun_Time << std::endl;
		std::cout << "Host<-->Device data copy time: " << dataCopy_Time << std::endl;
		std::cout << "Outptus write time: " << outputWrite_Time << std::endl;
		std::cout << "Total time of including overhead :  " << TotalTime << std::endl;
		std::cout << "Done! return value: " << retvalue << std::endl;
		//fflush(stdout);
	}
	//std::cout << "Done! return value: " << retvalue << std::endl;
	//fflush(stdout);
	MPI::Finalize();
	getchar();
	return 0;
}
