//copied (and modified) from snowdv.f90  
#include "uebdafunctions.h"
#include <ctime>


// read forcing data assimilation control file
void uebEnKF::readDaContr(const char* daconFile)
{
	std::ifstream pinFile(daconFile);
	char headerLine[256];
	pinFile.getline(headerLine, 256);   //skip header
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%d ", &es_enseSize);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f ", &forcEnStdev);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f ", &tempEnStdev);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f ", &forcCorLength);

/*	for (int is = 0; is < danumAssmnStates; is++) {
		pinFile.getline(headerLine, 256, '\n');
		sscanf(headerLine, "%f ", &daStatesStdev[is]);
	}
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f ", &stateCorLength);*/

	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f %f ", &dyC, &dxC);
	
	pinFile.getline(headerLine, 256);   // counter line ---to be updated 7.15.16
	sscanf(headerLine, "%d ", &nDaStates);
	daContArr.resize(nDaStates);
	for (int ida = 0; ida < nDaStates; ida++)
	{
		pinFile.getline(headerLine, 256, ':');
		sscanf(headerLine, "%s ", &daContArr[ida].stateVar);
		pinFile.getline(headerLine, 256, '\n'); //skip the remaining part of the line after :
		pinFile.getline(headerLine, 256, '\n');
		sscanf(headerLine, "%d ", &daContArr[ida].infType);
		switch (daContArr[ida].infType)
		{
		case 0:
			pinFile.getline(headerLine, 256, '\n');
			sscanf(headerLine, "%s ", &daContArr[ida].stateInpFile);
			pinFile.getline(headerLine, 256, '\n');
			sscanf(headerLine, "%f ", &daContArr[ida].obsErrStdev);
			pinFile.getline(headerLine, 256, '\n');
			sscanf(headerLine, "%f ", &daContArr[ida].obsCorLength);
			break;
		case 1:
			pinFile.getline(headerLine, 256, '\n');
			sscanf(headerLine, "%s %s %s %d", &daContArr[ida].stateInpFile, &daContArr[ida].infvarName, &daContArr[ida].inftimeVar, &daContArr[ida].numNcfiles);
			pinFile.getline(headerLine, 256, '\n');
			sscanf(headerLine, "%f ", &daContArr[ida].obsErrStdev);
			pinFile.getline(headerLine, 256, '\n');
			sscanf(headerLine, "%f ", &daContArr[ida].obsCorLength);
			break;
		default:
			std::cout << "Wrong data assimilation type; has to be 0 (time-series text file) or 1 (netcdf)" << std::endl;
			//std::cout << "Using default value..." << std::endl;
			getchar();
			//exit(1);
			break;
		}		
	}
	pinFile.close();
	return;
}

void uebEnKF::initDAMatrices()            //std::vector<std::pair<int, int> > cellCoordinates)  
{
	ng_gridSize = cellCoordinates.size();
	ns_statSize = danumAssmnStates;
	mo_obseSize = 0;
	for (int io = 0; io < nDaStates; io++) {
		//resise matrices to target size
		if (daContArr[io].infType == 0)
			mo_obseSize += daContArr[io].m_numObs;
		else
			mo_obseSize += cellCoordinates.size();          //if gridded (nc) observation then the m = number of active grid cells
	}
	//resise matrices to target size
	
/*	//##To do revise this later
     Q_modelErrCov.resize(ns_statSize, ns_statSize);
	// not used yet vQ_modelErr.resize(nd_xSize * ns_statSize, es_enseSize); 	
	Q_modelErrCov.setZero();  // setOnes();
	//Q_modelErrCov *= daContArr.forcEnStdev*daContArr.forcEnStdev;
*/
	// forcing covariance matrix based on distance between grid cells
	//for perturbation of forcing
	meanF.resize(ng_gridSize);
	//for all forcing except temp, mean of 1 with std.dev from user input
	meanF.setOnes();	
	covarF.resize(ng_gridSize, ng_gridSize);
	double rij = 0;      //distance between grid cells 
	double corij = 1;    //correlation between grid cells i and j
	// covariance matrix based on distance between grid cells
	for (int i = 0; i < ng_gridSize; i++) {
		for (int j = 0; j < ng_gridSize; j++) {
			rij = (cellCoordinates[i].first - cellCoordinates[j].first)*(cellCoordinates[i].first - cellCoordinates[j].first) * dyC * dyC    //(i1-i2)^2 * dy^2 + (j1-j2)^2 *dx^2
				+ (cellCoordinates[i].second - cellCoordinates[j].second)*(cellCoordinates[i].second - cellCoordinates[j].second) * dxC * dxC;
			rij = sqrt(rij);
			corij = exp(-1.0*rij / forcCorLength);
			covarF(i, j) = corij;							//*daContArr.forcEnStdev * daContArr.forcEnStdev;
		}
	}	
	Matrix<double, Dynamic, Dynamic, RowMajor> covarFS;
	covarFS = covarF * forcEnStdev * forcEnStdev;
	// Create a multi variate standard normal distribution 
	const uint64_t seedF = static_cast<uint64_t>(time(0));
	//9.1.16 set rand. generator
	std_norm_dist_Forc_Default.setSeed(seedF);
	std_norm_dist_Forc_Default.setMean(meanF);
	std_norm_dist_Forc_Default.setCovar(covarFS, true);

	//for temperature use additive termwith mean 0
	meanT.resize(ng_gridSize);
	//for all forcing except temp, mean of 1 with std.dev from user input
	meanT.setZero();
	Matrix<double, Dynamic, Dynamic, RowMajor> covarTS;
	covarTS = covarF * tempEnStdev * tempEnStdev;
	// Create a multi variate standard normal distribution 
	const uint64_t seedT = static_cast<uint64_t>(time(0));
	norm_dist_0Mean_Tempr.setSeed(seedT);
	norm_dist_0Mean_Tempr.setMean(meanT);
	norm_dist_0Mean_Tempr.setCovar(covarTS, true);

	if (useHmatrix) {
		Hc_hgVector.resize(mo_obseSize);
		Hc_hsVector.resize(mo_obseSize);
		/*	H_hMaxtirx.resize(mo_obseSize, nds_XstatesSize); //use this if not using function hX()
			//initialize with 0s
			H_hMaxtirx.setZero();
		*/
	}
	Z_obs.resize(mo_obseSize);
	R_obsErrCov.resize(mo_obseSize, mo_obseSize);
	//not used yet vR_obsErr.resize(mo_obseSize, es_enseSize);
//##TODO Revise later
	R_obsErrCov.setZero(); // = covarM.cast<float>() * daContArr.obsErrStdev * daContArr.obsErrStdev;
	int im = 0;
	for (int io = 0; io < nDaStates; io++)
	{
		if (daContArr[io].infType == 0)
		{
			//point observations
			for (int ida = 0; ida < daContArr[io].m_numObs; ++ida) {
				R_obsErrCov(im, im) = daContArr[io].obsErrStdev * daContArr[io].obsErrStdev;
				im++;
			}				
		}
		else
		{
			//gridded observations
			for (int jn = 0; jn < ng_gridSize; ++jn) {
				R_obsErrCov(im, im) = daContArr[io].obsErrStdev * daContArr[io].obsErrStdev;
				im++;
			}
		}
	}
	//for measurement / observatin 
	meanM.resize(mo_obseSize);
	// Create a multi variate normal distribution with mean 0                    8.28.16
	meanM.setZero();  //8.28.16                             //8.28.16 for obs use y' = y + vR, vR ~ N(0,R)	
	const uint64_t seedM = static_cast<uint64_t>(time(0));
	norm_dist_0Mean_Default.setSeed(seedM);
	norm_dist_0Mean_Default.setMean(meanM);
	norm_dist_0Mean_Default.setCovar(R_obsErrCov);

	startIndexDA.resize(nDaStates);
	ncReadStartDA.resize(nDaStates);
	for (int io = 0; io < nDaStates; io++) {
		startIndexDA[io] = 0;
		ncReadStartDA[io] = 0;
	}
	//8.8.16 3th state for snow surface temp
	stateIndex.resize(nDaStates);
	for (int io = 0; io < nDaStates; io++) {
		for (int vindx = 0; vindx < ns_statSize; vindx++) {    //TODO: try to do without looping?  
			if (strcmp(daContArr[io].stateVar, uebdaStates[vindx]) == 0) {
				stateIndex[io] = vindx;
				break;
			}
		}
	}
}
void uebEnKF::setHMatrices(float **daYcorrArr, float **daXcorrArr)
{
	int im = 0;
	for (int io = 0; io < nDaStates; io++)
	{
		if (daContArr[io].infType == 0)
		{
			//point observations
			for (int ida = 0; ida < daContArr[io].m_numObs; ++ida)
				for (int jn = 0; jn < ng_gridSize; ++jn)
				{
					if (abs(daYcorrArr[io][ida] - (y0 + cellCoordinates[jn].first * dyC)) < 0.1 * dyC &&
						abs(daXcorrArr[io][ida] - (x0 + cellCoordinates[jn].second * dxC)) < 0.1 * dxC)
					{
						H_hMaxtirx(im, jn * ns_statSize + stateIndex[io]) = 1;
						im++;
						break;
					}
				}
		}
		else 
		{
			//gridded observations
			for (int jn = 0; jn < ng_gridSize; ++jn)
				H_hMaxtirx(im, jn * ns_statSize + stateIndex[io]) = 1;				
			im++;
		}
	}
	//std::cout << std::endl<<" H Matrix: "<<std::endl << H_hMaxtirx << std::endl;
}
//this finds the indices of the grid cells where there are observations
void uebEnKF::setHcMatrices(float **daYcorrArr, float **daXcorrArr)
{
	int im = 0;
	for (int io = 0; io < nDaStates; io++)
	{
		if (daContArr[io].infType == 0)
		{
			//point observations
			for (int ida = 0; ida < daContArr[io].m_numObs; ++ida)
			{
				for (int jn = 0; jn < ng_gridSize; ++jn)
				{
					if (abs(daYcorrArr[io][ida] - (y0 + cellCoordinates[jn].first * dyC)) < 0.1 * dyC &&
						abs(daXcorrArr[io][ida] - (x0 + cellCoordinates[jn].second * dxC)) < 0.1 * dxC)
					{
						Hc_hgVector(im) = jn;    // *ns_statSize + stateIndex[io];  the index of the grid cell where there is observation
						Hc_hsVector(im) = stateIndex[io];
						im++;
						break;
					}
				}
			}
		}
		else
		{
			//gridded observations
			for (int jn = 0; jn < ng_gridSize; ++jn) 
			{
				Hc_hgVector(im) = jn;         // *ns_statSize + stateIndex[io];
				Hc_hsVector(im) = stateIndex[io];
				im++;
			}
		}
	}
	//std::cout << std::endl<<" H Matrix: "<<std::endl << H_hMaxtirx << std::endl;
}
/*__host__ __device__*/
//using defulat generator--no seeding, setupuse mean of 1s and covariance matrix from input std.dev to sample points from multivar std. normal dist.
void uebEnKF::getMultivarStdNorm_Samples_Forc_Default(float** &stdNormArr)   //for now we use this as am not sure if we can use eigen matrix in mpi broadcast
{	
	//sample ndim * nens points ---ndim vectors of size nens
	Matrix<double, Dynamic, Dynamic, RowMajor> stdnormSamples(ng_gridSize, es_enseSize);
	//Matrix<double, Dynamic, -1> 
	stdnormSamples = std_norm_dist_Forc_Default.samples(es_enseSize);
	//std::cout<< stdnormSamples << std::endl;
	for (int i = 0; i < ng_gridSize; i++)                                             //8.16.16 TBCL     this is an inefficient method
		for (int j = 0; j < es_enseSize; j++)
			stdNormArr[i][j] = stdnormSamples(i, j);

	//if (snowdgt_outflag == 1) {       //if debug mode 
	std::ofstream stdnormalSamplestxt("ensForcingMultipliers.txt");
	stdnormalSamplestxt << stdnormSamples.transpose() << std::endl;
	//}
	return;
}
//using defulat generator--no seeding, setupuse mean of 1s and covariance matrix from input std.dev to sample points from multivar std. normal dist.
void uebEnKF::getMultivarStdNorm_Samples_Forc_Tempr(float** &stdNormArr)   //for now we use this as am not sure if we can use eigen matrix in mpi broadcast
{
	//sample ndim * nens points ---ndim vectors of size nens
	Matrix<double, Dynamic, Dynamic, RowMajor> stdnormSamples(ng_gridSize, es_enseSize);
	//Matrix<double, Dynamic, -1> 
	stdnormSamples = norm_dist_0Mean_Tempr.samples(es_enseSize);
	//std::cout<< stdnormSamples << std::endl;
	for (int i = 0; i < ng_gridSize; i++)                                             //8.16.16 TBCL     this is an inefficient method
		for (int j = 0; j < es_enseSize; j++)
			stdNormArr[i][j] = stdnormSamples(i, j);

	//if (snowdgt_outflag == 1) {       //if debug mode 
	std::ofstream stdnormalSamplestxt("ensTempForcing_err.txt");
	stdnormalSamplestxt << stdnormSamples.transpose() << std::endl;
	//}
	return;
}
/*__host__ __device__*/
//this function gets the obs error samples (not multipliers) using default generator
void uebEnKF::getMultivarNorm_Samples_Obs_Default(float** &multiNormSampleArr)
{
	//sample ndim * nens points ---ndim vectors of size nens
	Matrix<double, Dynamic, Dynamic, RowMajor> stdnormSamples(mo_obseSize, es_enseSize);
	//Matrix<double, Dynamic, -1> 
	stdnormSamples = norm_dist_0Mean_Default.samples(es_enseSize);
	//if (snowdgt_outflag == 1)
	//std::cout << stdnormSamples << std::endl;
	for (int i = 0; i < mo_obseSize; i++)                                             //8.16.16 TBCL     this is an inefficient method
		for (int j = 0; j < es_enseSize; j++)
			multiNormSampleArr[i][j] = stdnormSamples(i, j);
	//if (snowdgt_outflag == 1) {       //if debug mode 
	std::ofstream stdnormalSamplestxt("ensObsErr.txt");
	stdnormalSamplestxt << stdnormSamples.transpose() << std::endl;
	//}
	return;
}
/*__host__ __device__*/ 
void uebEnKF::getEnKFArrays(float** Xh_stateObsSpace, float** &X1_o_ensAnomalyObsSpace, float **&y_obsStateResidual, float **&Pzz_i_ObsStateCov_inv)
{
	VectorXd ensMeanArr(mo_obseSize);
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> Xh_obsState(mo_obseSize, es_enseSize);
	for (int io = 0; io < mo_obseSize; io++) {
		ensMeanArr(io) = 0.0;
		for (int ie = 0; ie < es_enseSize; ie++) {
			ensMeanArr(io) += Xh_stateObsSpace[io][ie];
			Xh_obsState(io, ie) = Xh_stateObsSpace[io][ie];
		}
	}
	ensMeanArr /= es_enseSize;
	//ensemble anomaly
	Matrix<float, Dynamic, Dynamic, RowMajor> ensAnomalyArr(mo_obseSize, es_enseSize);
	for (int id = 0; id < mo_obseSize; id++)
		for (int ie = 0; ie < es_enseSize; ie++) {
			ensAnomalyArr(id, ie) = Xh_obsState(id, ie) - ensMeanArr(id);
			X1_o_ensAnomalyObsSpace[id][ie] = Xh_obsState(id, ie) - ensMeanArr(id);
		}
/*	std::cout << std::endl << " States matrix: " << std::endl;
	std::cout << Xh_obsState << " " << std::endl;
	std::cout << std::endl << " Ensemble mean: " << std::endl;
	std::cout << ensMeanArr << " " << std::endl;
	std::cout << std::endl << " Ensemble anomaly: " << std::endl;
	std::cout << ensAnomalyArr << " " << std::endl;
*/	//observed array 
	std::cout << std::endl << " observed array: " << std::endl;
	std::cout << Z_obs << " " << std::endl;
//	std::cout << std::endl << " model forecasted observation array: " << std::endl;
//	std::cout << Xh_obsState << std::endl;
	//(observed - model forecasted observation), correction / residual/ innovation y = z - HXb (m * Nens) 
	//Yi = Y + vR, vR ~ N(0,R)
	float** ensObservationErr = new float*[mo_obseSize];
	for (int i = 0; i < mo_obseSize; i++)
		ensObservationErr[i] = new float[es_enseSize];
	getMultivarNorm_Samples_Obs_Default(ensObservationErr);
	//if (snowdgt_outflag == 1)        //if debug mode 
/*	std::cout << std::endl << "observation error: " << std::endl;
	for (int id = 0; id < mo_obseSize; id++)
		for (int is = 0; is < es_enseSize; is++)
			std::cout << ensObservationErr[id][is] << " ";
	std::cout << std::endl;
*/	//residual y = z - HXb and R' = function of Yobs
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> y_obsStateRresidual(mo_obseSize, es_enseSize);
	for (int id = 0; id < mo_obseSize; id++)
		for (int is = 0; is < es_enseSize; is++) {
			y_obsStateResidual[id][is] = Z_obs(id) + ensObservationErr[id][is] - Xh_obsState(id, is);
			y_obsStateRresidual(id, is) = Z_obs(id) + ensObservationErr[id][is] - Xh_obsState(id, is);   // *ensObservationMultiplier[id][is] - Xh_obsState(id, is);
		}
/*	std::cout << std::endl << " correction / residulal / Innovation (observed - model forecasted observation) y " << std::endl;
	std::cout << y_obsStateRresidual << std::endl;
	std::cout << std::endl << " R observation Error covariance matrix: " << std::endl;
	std::cout << R_obsErrCov << " ";
*/	//K = Pb * H_T(H*Pb*H_T + R)^-1, Pzz observation uncertainty inovation matrix
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>
		Pzz_obsStateCov = ensAnomalyArr * ensAnomalyArr.transpose();
	Pzz_obsStateCov /= (es_enseSize - 1);
	Pzz_obsStateCov += R_obsErrCov.cast<float>();
//	std::cout << std::endl << " Pzz Innovation covariance matrix(Observation uncertainty) : " << std::endl;
//	std::cout << Pzz_obsStateCov << std::endl;
	//PB_R_1 inverse of Observation uncertainty (innovation?) matrix
	Matrix<float, Dynamic, Dynamic, RowMajor> Pzz_i(mo_obseSize, mo_obseSize);
	float PBR_det = Pzz_obsStateCov.determinant();
	if (PBR_det == 0) {
		std::cout << std::endl << " Warnning: zero determinant of matrix!" << PBR_det << std::endl;
		debugOutputFile.open("debugOutput.txt", std::ios::app);
		debugOutputFile << " Warnning: zero determinant of matrix!" << PBR_det << std::endl;
		debugOutputFile.close();
		std:getchar();
	}
	else
		Pzz_i = Pzz_obsStateCov.inverse();
	//
	for (int io = 0; io < mo_obseSize; io++)
		for (int im = 0; im < mo_obseSize; im++)
			Pzz_i_ObsStateCov_inv[io][im] = Pzz_i(io, im);
	std::cout << std::endl << " Pzzi Inverse of Pzz (Observation uncertainty / innovation covariance) matrix: " << std::endl;
	std::cout << Pzz_i << std::endl;
	
	//free memory	
	for (int i = 0; i < mo_obseSize; i++) {
		delete[] ensObservationErr[i];
		ensObservationErr[i] = NULL;
	}
	delete[] ensObservationErr;
	ensObservationErr = NULL;

	return;
}
//
void uebEnKF::runEnKF(float** X1_o_ensAnomalyObsSpace, float **y_obsStateResidual, float **Pzz_i_ObsStateCov_inv, float** stateOutputArr, float** &stateOutputUpdate)
{
	VectorXd ensMeanArr(ns_statSize);
	Matrix<float, Dynamic, Dynamic, RowMajor> X_stateBackground(ns_statSize, es_enseSize);
	for (int id = 0; id < ns_statSize; id++) {
		ensMeanArr(id) = 0.0;
		for (int ie = 0; ie < es_enseSize; ie++) {
			ensMeanArr(id) += stateOutputArr[id][ie];
			X_stateBackground(id, ie) = stateOutputArr[id][ie];
		}
		ensMeanArr(id) /= es_enseSize;
	}
	//ensemble anomaly
	Matrix<float, Dynamic, Dynamic, RowMajor> ensAnomalyArr(ns_statSize, es_enseSize);
	for (int id = 0; id < ns_statSize; id++)
		for (int ie = 0; ie < es_enseSize; ie++)
			ensAnomalyArr(id, ie) = X_stateBackground(id, ie) - ensMeanArr(id);
	/*	std::cout << std::endl << " States matrix: " << std::endl;
	std::cout << X_stateBackground << " " << std::endl;
	std::cout << std::endl << " Ensemble mean: " << std::endl;
	std::cout << ensMeanArr << " " << std::endl;
	std::cout << std::endl << " Ensemble anomaly: " << std::endl;
	std::cout << ensAnomalyArr << " " << std::endl;
	*/
	//residual y = z - HXb and R' = function of Yobs
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> y_obsStateRresidual(mo_obseSize, es_enseSize);
	for (int id = 0; id < mo_obseSize; id++)
		for (int ie = 0; ie < es_enseSize; ie++)
			y_obsStateRresidual(id, ie) = y_obsStateResidual[id][ie];
	//	std::cout << std::endl << " correction / residulal / Innovation (observed - model forecasted observation) y " << std::endl;
	//	std::cout << y_obsStateRresidual << std::endl;
	// Pxz state obs cross-covariance
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> Pxz_stateObsXCov(ns_statSize, mo_obseSize);
	for (int id = 0; id < ns_statSize; id++)
		for (int io = 0; io < mo_obseSize; io++) {
			Pxz_stateObsXCov(id, io) = 0.0;
			for (int ie = 0; ie < es_enseSize; ie++)
				Pxz_stateObsXCov(id, io) += (ensAnomalyArr(id, ie) * X1_o_ensAnomalyObsSpace[io][ie]);    //ensAnomalyObsSpace[io][ie]) ie changing fast == transpose
		}
	Pxz_stateObsXCov /= (es_enseSize - 1);
	//	std::cout << std::endl << " Pxz state-obs cross-covariance matrix: " << std::endl;
	//	std::cout << Pxz_stateObsXCov << std::endl;
	//Kalman gain K = Pxz * Pzz^-1 
	Matrix<float, Dynamic, Dynamic, RowMajor> K_kalGain(ns_statSize, mo_obseSize);
	for (int id = 0; id < ns_statSize; id++)
		for (int ik = 0; ik < mo_obseSize; ik++) {
			K_kalGain(id, ik) = 0.0;
			for (int io = 0; io < mo_obseSize; io++)
				K_kalGain(id, ik) += (Pxz_stateObsXCov(id, io) * Pzz_i_ObsStateCov_inv[io][ik]);
		}
	//	std::cout << std::endl << " K Kalman gain: " << std::endl;
	//	std::cout << K_kalGain << std::endl;
	//update state Xa = Xb + K*y, y = z-HXb
	Matrix<float, Dynamic, Dynamic, RowMajor>
		X_state = X_stateBackground + (K_kalGain * y_obsStateRresidual);
	//std::cout << std::endl << " Updated state matrix: " << std::endl;
	//std::cout << X_state << std::endl;
	//	std::cout << std::endl << " Updated state matrix mean: " << std::endl;
	float ensupdateMean = 0.0;
	for (int id = 0; id < ns_statSize; id++) {
		ensupdateMean = 0.0;
		for (int ie = 0; ie < es_enseSize; ie++) {
			ensupdateMean += X_state(id, ie);
			stateOutputUpdate[id][ie] = X_state(id, ie);
		}
		stateOutputUpdate[id][es_enseSize] = ensupdateMean / es_enseSize;
		//	    std::cout << stateOutputUpdate[id][es_enseSize] << " ";
	}
	//std::cout << std::endl;

	return;
}
void uebEnKF::runEnKF(float** Xh_stateObsSpace, float** ensObservationErr, float** stateOutputArr, float** &stateOutputUpdate)         // bool NormalDist)   //float* &ensembleUpdateArr,
{
	VectorXf XhensMeanArr(mo_obseSize);
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> Xh_obsState(mo_obseSize, es_enseSize);
	for (int io = 0; io < mo_obseSize; io++) {
		XhensMeanArr(io) = 0.0;
		for (int ie = 0; ie < es_enseSize; ie++) {
			XhensMeanArr(io) += Xh_stateObsSpace[io][ie];
			Xh_obsState(io, ie) = Xh_stateObsSpace[io][ie];
		}
	}
	XhensMeanArr /= es_enseSize;
	//ensemble anomaly
	Matrix<float, Dynamic, Dynamic, RowMajor> Xh_ensAnomalyArr(mo_obseSize, es_enseSize);
	for (int id = 0; id < mo_obseSize; id++)
		for (int ie = 0; ie < es_enseSize; ie++)
			Xh_ensAnomalyArr(id, ie) = Xh_obsState(id, ie) - XhensMeanArr(id);
	if (uebda_debugout2 == 1)
	{
		std::cout << std::endl << " States in Obs space: " << std::endl;
		std::cout << Xh_obsState << " " << std::endl;
		std::cout << std::endl << " Ensemble mean in Obs space: " << std::endl;
		std::cout << XhensMeanArr << " " << std::endl;
		std::cout << std::endl << " Xh' Ensemble anomaly in Obs Space : " << std::endl;
		std::cout << Xh_ensAnomalyArr << " " << std::endl;
	}
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>
		Pzz_obsStateCov = Xh_ensAnomalyArr * Xh_ensAnomalyArr.transpose();
	Pzz_obsStateCov /= (es_enseSize - 1);
	Pzz_obsStateCov += R_obsErrCov.cast<float>();
	if (uebda_debugout2 == 1)
	{
		std::cout << std::endl << " R observation Error covariance matrix: " << std::endl;
		std::cout << R_obsErrCov << " ";
	}
	if (uebda_debugout2 == 1)
	{
		std::cout << std::endl << " Pzz Innovation covariance matrix(Observation uncertainty) : " << std::endl;
		std::cout << Pzz_obsStateCov << std::endl;
	}
	//Pzz_i inverse of Observation uncertainty (innovation?) matrix
	Matrix<float, Dynamic, Dynamic, RowMajor> Pzz_i(mo_obseSize, mo_obseSize);
	float Pzz_det = Pzz_obsStateCov.determinant();
	if (Pzz_det == 0) {
		std::cout << std::endl << " Warnning: zero determinant of matrix! Press 'Enter' to continue with Kalman Gain = 0" << Pzz_det << std::endl;
		debugOutputFile.open("debugOutput.txt", std::ios::app);
		debugOutputFile << " Warnning: zero determinant of matrix!" << Pzz_det << std::endl;
		debugOutputFile.close();
	    std::getchar();
	}
	else
		Pzz_i = Pzz_obsStateCov.inverse();

	if (uebda_debugout2 == 1)
	{
		std::cout << std::endl << " Pzzi Inverse of Pzz (Observation uncertainty / innovation covariance) matrix: " << std::endl;
		std::cout << Pzz_i << std::endl;
	}	
	//states
	VectorXf ensMeanArr(ns_statSize);
	Matrix<float, Dynamic, Dynamic, RowMajor> X_stateBackground(ns_statSize, es_enseSize);
	for (int id = 0; id < ns_statSize; id++) {
		ensMeanArr(id) = 0.0;
		for (int ie = 0; ie < es_enseSize; ie++) {
			ensMeanArr(id) += stateOutputArr[id][ie];
			X_stateBackground(id, ie) = stateOutputArr[id][ie];
		}
	}
	ensMeanArr /= es_enseSize;
	//ensemble anomaly
	Matrix<float, Dynamic, Dynamic, RowMajor> ensAnomalyArr(ns_statSize, es_enseSize);
	for (int id = 0; id < ns_statSize; id++)
		for (int ie = 0; ie < es_enseSize; ie++) 
			ensAnomalyArr(id, ie) = X_stateBackground(id, ie) - ensMeanArr(id);
	if (uebda_debugout2 == 1)
	{
		std::cout << std::endl << " States matrix: " << std::endl;
		std::cout << X_stateBackground << " " << std::endl;
		std::cout << std::endl << " Ensemble mean: " << std::endl;
		std::cout << ensMeanArr << " " << std::endl;
		std::cout << std::endl << " Ensemble anomaly: " << std::endl;
		std::cout << ensAnomalyArr << " " << std::endl;
	}
	// Pxz state obs cross-covariance
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> 
	Pxz_stateObsXCov = ensAnomalyArr * Xh_ensAnomalyArr.transpose();
	Pxz_stateObsXCov /= (es_enseSize - 1);
	if (uebda_debugout2 == 1)
	{
		std::cout << std::endl << " Pxz state-obs cross-covariance matrix: " << std::endl;
		std::cout << Pxz_stateObsXCov << std::endl;
		std::cout << std::endl << " observed array: " << std::endl;
		std::cout << Z_obs << " " << std::endl;
	}
	//(observed - model forecasted observation), correction / residual/ innovation y = z - HXb (m * Nens) 
	//Yi = Y + vR, vR ~ N(0,R)
/*  float** ensObservationErr = new float*[mo_obseSize];
	for (int i = 0; i < mo_obseSize; i++)
		ensObservationErr[i] = new float[es_enseSize];
	getMultivarNorm_Samples_Obs_Default(ensObservationErr);
	*/
	//if (snowdgt_outflag == 1)        //if debug mode 
	if (uebda_debugout2 == 1)
	{
		std::cout << std::endl << "observation error: " << std::endl;
		for (int id = 0; id < mo_obseSize; id++)
			for (int is = 0; is < es_enseSize; is++)
				std::cout << ensObservationErr[id][is] << " ";
		std::cout << std::endl;
	}
	//residual y = z - HXb and R' = function of Yobs
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> y_obsStateRresidual(mo_obseSize, es_enseSize);
	for (int id = 0; id < mo_obseSize; id++) 
		for (int is = 0; is < es_enseSize; is++) 
			y_obsStateRresidual(id, is) = Z_obs(id) + ensObservationErr[id][is] - Xh_obsState(id, is);   //
	if (uebda_debugout2 == 1)
	{
		std::cout << std::endl << " correction / residulal / Innovation (observed - model forecasted observation) y " << std::endl;
		std::cout << y_obsStateRresidual << std::endl;
	}
	//Kalman gain
	Matrix<float, Dynamic, Dynamic, RowMajor> K_kalGain(ns_statSize, mo_obseSize);
	if (Pzz_det == 0) {
		std::cout << std::endl << " Warnning: zero determinant of matrix! Press 'enter' to continue with Kalman Gain = 0" << Pzz_det << std::endl;
		K_kalGain.setZero();                                 //8.28.16 effect is unknown---so don't use update
	    std::getchar();
		goto label1;	    
	}
	//K Kalman gain K = Pxz * Pzz^-1 = Pxz*Pb_R_1
	
	K_kalGain = Pxz_stateObsXCov * Pzz_i;
//8.28.16 ----what if Kalman gain is = 0 for long time?
label1:
	if (uebda_debugout2 == 1)
	{
		std::cout << std::endl << " K Kalman gain: " << std::endl;
		std::cout << K_kalGain << std::endl;
	}
	//update state Xa = Xb + K*y, y = z-HXb
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> X_state(ns_statSize, es_enseSize);
	X_state = X_stateBackground + (K_kalGain * y_obsStateRresidual);
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> Pxx_stateCov(ns_statSize, ns_statSize);
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> Pxx_stateCovUpdate(ns_statSize, ns_statSize);
	//background Pxx = (1 / (es_enseSize - 1)) *ensAnomalyArr * ensAnomalyArr.transpose();
	Pxx_stateCov = ensAnomalyArr * ensAnomalyArr.transpose();
	Pxx_stateCov /= (es_enseSize - 1);
	Pxx_stateCovUpdate = Pxx_stateCov - K_kalGain * Pzz_obsStateCov * K_kalGain.transpose();
	if (uebda_debugout2 == 1)
	{
		std::cout << std::endl << " Updated state matrix: " << std::endl;
		std::cout << X_state << std::endl;
		std::cout << std::endl << " Pxx Input model error covariance: " << std::endl;
		std::cout << Pxx_stateCov << std::endl;
		std::cout << std::endl << " Pxx_u Updated model error covariance: " << std::endl;
		std::cout << Pxx_stateCovUpdate << std::endl;
	}

	//if SWE <= 0: SWE, tausn, dlSage, refD and trefD must be 0 
	for (int ie = 0; ie < es_enseSize; ie++)
	{
		if (X_state(1, ie) <= 0) 
		{
			X_state(1, ie) = 0.0;
			X_state(2, ie) = 0.0;
			X_state(3, ie) = 0.0;
			X_state(4, ie) = 0.0;
		}
	}
	if (uebda_debugout1 == 1)
		std::cout << std::endl << " Updated state matrix mean: " << std::endl;
	float ensupdateMean = 0.0;		
	for (int is = 0; is < ns_statSize; is++)
	{
		ensupdateMean = 0.0;
		for (int ie = 0; ie < es_enseSize; ie++) 
		{
			stateOutputUpdate[is][ie] = X_state(is, ie);          //collecting same variables together for convinience
			ensupdateMean += X_state(is, ie);
		}
		stateOutputUpdate[is][es_enseSize] = ensupdateMean / es_enseSize;
		if (uebda_debugout1 == 1)
			std::cout << stateOutputUpdate[is][es_enseSize] << " ";
	}
	if (uebda_debugout1 == 1)
		std::cout << std::endl;

	//free memory	
/*	for (int i = 0; i < mo_obseSize; i++) {
		delete[] ensObservationErr[i];
		ensObservationErr[i] = NULL;
	}
	delete[] ensObservationErr;
	ensObservationErr = NULL;
*/
	return;
}

void uebEnKF::getDaArr(int ida, int &numNc, float** RegArray, float &tcorVar, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//for (int it = 0; it < 2; it++){
	if (daContArr[ida].infType == 0)
	{
		// for time series from text file read once ---- outside of this function
	}
	else if (daContArr[ida].infType == 1)
	{
		int retvalue = 0;
		//read 3D netcdf (regridded array processed by uebInputs)
		char numtoStr[256];
		sprintf(numtoStr, "%d", numNc);
		char tsInputfile[256];
		strcpy(tsInputfile, daContArr[ida].stateInpFile);
		strcat(tsInputfile, numtoStr);
		strcat(tsInputfile, ".nc");
		//std::cout<<"%s\n",tsInputfile);
		//clear existing memory RegArray[it] before passing to this function // delete[] RegArray[it];
		retvalue = readNC_yxSlub_givenT(tsInputfile, daContArr[ida].infvarName, daContArr[ida].inftimeVar, ncReadStartDA[ida], RegArray, tcorVar, numNc, inpComm, inpInfo);
	}
	//}
}

void uebEnKF::updateDaArr(float*** RegArray, float tcorVar, double** tvarArr)          //, int nDataPoints) 
{
	int Year = YearDA, Month = MonthDA, Day = DayDA;
	double Hour = HourDA;
	//need to call each variable array as each array has to be copied separately to device array in cuda	
	int io = 0;
	for (int jd = 0; jd < nDaStates; jd++)
	{
		if (daContArr[jd].infType == 0)       //
		{
			for (int iot = 0; iot < daContArr[jd].m_numObs; iot++) {
				Z_obs[io] = RegArray[jd][startIndexDA[jd]][iot];
				io++;
			}
			daTime = tvarArr[jd][startIndexDA[jd]];
			startIndexDA[jd] +=1;
		}
		else if (daContArr[jd].infType == 1) 
		{
			for (int ig = 0; ig < ng_gridSize; io++) {            //8.22.16 in this case there is one value at each grid cell
				Z_obs[io] = RegArray[jd][cellCoordinates[ig].first][cellCoordinates[ig].second];
				io++;
			}
			UPDATEtime(Year, Month, Day, Hour, (double)tcorVar);       //8.5.16---for now assume tcorvar is in Hours since start time
			daTime = julian(Year, Month, Day, Hour);
		}
	}
}
// read input text file and record datetime, and list of values; skip no data, get no data value from file
void  uebEnKF::readTStextFile_multiVal(const char* inforcFile, int ida, float* &yCorArr, float* &xCorArr, double* &tvarArr, float** &varvalArr)  // int &nrecords)
{
	FILE* inputFile = fopen(inforcFile, "r");
	if (!inputFile)
	{
		std::cout << "Error opening file: " << inforcFile << std::endl;
		return;
	}
	int nrecords = 0;
	float noDataV = -9999;
	char commentLine[256];                    //string to read header line	
	fscanf(inputFile, "%f %d ", &noDataV, &daContArr[ida].m_numObs); //   get no data value, number of data cols
	fgets(commentLine, 256, inputFile);   //skip remaining line
	yCorArr = new float[daContArr[ida].m_numObs];
	xCorArr = new float[daContArr[ida].m_numObs];
	for (int id = 0; id < daContArr[ida].m_numObs; id++)
		fscanf(inputFile, "%f %f ", &yCorArr[id], &xCorArr[id]); //   coordinates of data points
	fgets(commentLine, 256, inputFile);   //skip remaining contents of line

	fgets(commentLine, 256, inputFile);   //skip header line 
	int Year, Month, Day;
	double Hour, DTimeV;
	float Value;
	while (!feof(inputFile))
	{
		commentLine[0] = ' ';
		fgets(commentLine, 256, inputFile);
		if (commentLine[0] != ' ') 
			++nrecords;
	}//while
	daContArr[ida].nRecs = nrecords;
	//strinpts = new inptimeseries[nrecords];                //assign memory to store data records
	tvarArr = new double[nrecords];
	varvalArr = new float*[nrecords];
	for (int ir = 0; ir < nrecords; ir++)
		varvalArr[ir] = new float[daContArr[ida].m_numObs];
	//
	rewind(inputFile);
	fgets(commentLine, 256, inputFile);   //no data value, number of data cols
	fgets(commentLine, 256, inputFile);   //coordinates of data points
	fgets(commentLine, 256, inputFile);   //skip header line 
	for (int ir = 0; ir<nrecords; ir++)
	{  
		fscanf(inputFile, "%d %d %d %lf %f ", &Year, &Month, &Day, &Hour, &Value);
		//std::cout << " hour " << Hour;
		if (fabs(Value - noDataV) > 0.1) {             //only copy data that is not no-data
			DTimeV = julian(Year, Month, Day, Hour);
			//std::cout << " hour julian " << std::setprecision(15)<< DTimeV;
			tvarArr[ir] = DTimeV;
			varvalArr[ir][0] =  Value;
			for (int id = 1; id <daContArr[ida].m_numObs; id++)
				fscanf(inputFile, "%f ", &varvalArr[ir][id]);
		}	 //
		//fscanf(inputFile, " %*s\n");
	}
	fclose(inputFile);

	return;
}
// read input text file and record datetime, value pair --skip no data, get no data value from file
void  uebEnKF::readTStextFileTimeValPair(const char* inforcFile, std::vector<std::pair<double, float>> &tvar_in, int &nrecords)
{
	std::ifstream inputFile(inforcFile, std::ios::in);
	if (!inputFile)
	{
		std::cout << "Error opening file: " << inforcFile << std::endl;
		return;
	}
	nrecords = 0;
	float noDataV = -9999;
	char commentLine[256];                    //string to read header line
	inputFile.getline(commentLine, 256, '\n');
	sscanf(commentLine, "%f \n", &noDataV); //   get no data value from file

	inputFile.getline(commentLine, 256, '\n');  //skip header line 
	int Year, Month, Day;
	double Hour, DTimeV;
	float Value;
	while (!inputFile.eof())
	{
		commentLine[0] = ' ';
		inputFile.getline(commentLine, 256, '\n');
		if (commentLine[0] != ' ') {  //condition to make sure empty line is not read; 
			sscanf(commentLine, "%d %d %d %lf %f \n", &Year, &Month, &Day, &Hour, &Value); //
																						   //std::cout << " hour " << Hour;
			if (fabs(Value - noDataV) > 0.1) {             //only copy data that is not no-data
				DTimeV = julian(Year, Month, Day, Hour);
				//std::cout << " hour julian " << std::setprecision(15)<< DTimeV;
				tvar_in.push_back(std::make_pair(DTimeV, Value));
				++nrecords;
			}
		}
	}//while
	inputFile.close();

	return;
}

//The following copied from snowxv 8.22.16
/*__host__ __device__*/   
void  uebEnKF::UPDATEtime(int &YEAR, int &MONTH, int &DAY, double &HOUR, double DT)
{
	int DM;				 // 30/03/2004 ITB 
						 // 30/03/2004 ITB  
						 //real  hour, dt  // DGT Dec 10, 2004.  Fixing ITB errors 

	int DMON[12] = { 31,28,31,30,31,30,31,31,30,31,30,31 };
	HOUR = HOUR + DT;
	DM = DMON[MONTH - 1];
	//  check for leap years 
	if (MONTH == 2)
		DM = lyear(YEAR);
	while (HOUR >= 24.0)
	{
		HOUR = HOUR - 24.0;
		DAY++;
	}

	while (DAY > DM)
	{
		DAY = DAY - DM;
		MONTH++;
		if (MONTH>12) {
			MONTH = 1;
			YEAR++;
		}
		//modified from the original by separating the above two lines in the if (month>12)
		//#_6.27.13
		DM = DMON[MONTH - 1];
		if (MONTH == 2)
			DM = lyear(YEAR);
		//}
	}
	return;
}

// ************************** lyear () ***************************
//    function to return number of days in February checking for leap years
/*__host__ __device__*/  
int uebEnKF::lyear(int year)
{
	int lyear; // Leap years are every 4 years 
			   // - except for years that are multiples of centuries (e.g. 1800, 1900)
			   // - except again that when the century is divisible by 4 (e.g. 1600, 2000)
	if ((year % 4  > 0) || ((year % 100 == 0) && (year % 400 != 0)))
		lyear = 28;
	else
		lyear = 29;
	return lyear;
}
//***************************** JULIAN () ****************************
//             To convert the real date to julian date
// YJS The Julian are change to a new version to take the Leap Yean into consideration
//    in the old version, there are 365 days each year.
//     FUNCTION JULIAN(MONTH,DAY)
/*__host__ __device__*/  
int uebEnKF::julian(int yy, int mm, int dd)
{
	int julian;
	int mmstrt[12] = { 0,31,59,90,120,151,181,212,243,273,304,334 };
	int jday = mmstrt[mm - 1] + dd;
	int ileap = yy - ((int)(yy / 4)) * 4;
	if ((ileap == 0) && (mm >= 3))
		jday = jday + 1;
	julian = jday;
	return julian;
}
//The following were copied from functions.f90 # 6.8.13
//THIS SUBROUTINE COMPUTES JULIAN DATE, GIVEN CALENDAR DATE AND time.  INPUT CALENDAR DATE MUST BE GREGORIAN.  INPUT time VALUE
//CAN BE IN ANY UT-LIKE time SCALE (UTC, UT1, TT, ETC.) - OUTPUT. //JULIAN DATE WILL HAVE SAME BASIS.  
//ALGORITHM BY FLIEGEL AND //VAN FLANDERN. //SOURCE: http://aa.usno.navy.mil/software/novas/novas_f/novasf_intro.php
//I = YEAR (IN) //M = MONTH NUMBER (IN) //K = DAY OF MONTH (IN) //H = UT HOURS (IN) //TJD = JULIAN DATE (OUT)
/*__host__ __device__*/ 
double uebEnKF::julian(int I, int M, int K, double H)
{
	double TJD, JD;
	//JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
	JD = K - 32075 + 1461 * (I + 4800 + (M - 14) / 12) / 4 + 367 * (M - 2 - (M - 14) / 12 * 12) / 12 - 3 * ((I + 4900 + (M - 14) / 12) / 100) / 4;
	TJD = JD - 0.5 + H / 24.0;
	//##%^_TBC 6.8.13 //powf(10,0) in place of D0
	return TJD;
}
//THIS SUBROUTINE COMPUTES CALENDAR DATE AND time, GIVEN JULIAN DATE.  INPUT JULIAN DATE CAN BE BASED ON ANY UT-LIKE time SCALE
//(UTC, UT1, TT, ETC.) - OUTPUT time VALUE WILL HAVE SAME BASIS. OUTPUT CALENDAR DATE WILL BE GREGORIAN.  
//ALGORITHM BY FLIEGEL AND VAN FLANDERN. //SOURCE: http://aa.usno.navy.mil/software/novas/novas_f/novasf_intro.php
//TJD = JULIAN DATE (IN) //I = YEAR (OUT) //M = MONTH NUMBER (OUT) //K = DAY OF MONTH (OUT) //H = UT HOURS (OUT)
/*__host__ __device__*/
void uebEnKF::calendardate(double TJD, int &I, int &M, int &K, double &H)
{
	double DJD, JD;
	int L, N;
	DJD = TJD + 0.5;
	JD = DJD;
	H = fmod(DJD, 1.0) * 24;    // 24.D0
								//JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
	L = JD + 68569;
	N = 4 * L / 146097;
	L = L - (146097 * N + 3) / 4;
	//I=YEAR, M=MONTH, K=DAY
	I = 4000 * (L + 1) / 1461001;
	L = L - 1461 * I / 4 + 31;
	M = 80 * L / 2447;
	K = L - 2447 * M / 80;
	L = M / 11;
	M = M + 2 - 12 * L;
	I = 100 * (N - 49) + I + L;
	return;
}
//from ncfunctions 8.22.16
//read multiple slubs (y,x arrays) along the time dim for data with t, y, x config---time as slowely varying array; and pvar_in already allocated
// __host__ __device__ 
int uebEnKF::readNC_yxSlub_givenT(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int &tStart, float** &pvar_in, float &tcorvar, int &numNc, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//float* pvarin_temp = NULL;
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0, ptid; // pxid = 0, pyid = 0, ndims = 0; 	
												  //variable data type
	nc_type varType;
	size_t pdim_sizes;
	//array of dimensions
	int pdimids[3]; //NC_MAX_DIMS]; 3D file only being read here; expected to get error message otherwise
					//dimension names 
	char pdim_Names[80];
	size_t start[3], count[3];
	//Open the netcdf file.  
	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
		ERR(retncval);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid)))
		ERR(retncval);
	// Get the varids of the coordinate variables 
	if ((retncval = nc_inq_varid(ncid, tcor_NAME, &ptid)))
		ERR(retncval);

	//var information, checking the dimension array
	if ((retncval = nc_inq_var(ncid, pvarid, NULL, &varType, NULL, pdimids, NULL)))
		ERR(retncval);
	//check dimension info and set start and count arrays; 
	//int yxDim = 1;
	int tEnd = 0;
	int tIndx = 0;
	for (int i = 0; i < 3; i++) {
		if (retncval = nc_inq_dim(ncid, pdimids[i], pdim_Names, &pdim_sizes))
			ERR(retncval);
		if (strcmp(pdim_Names, tcor_NAME) == 0) {
			tIndx = i;
			start[i] = tStart;
			count[i] = 1;
			tEnd = tStart + 1;
			if (tEnd < pdim_sizes) {
				tStart++;            // count[i];                //new start point
			}
			else {
				numNc++;                               //next time go to the next netcdf file
				tStart = 0;
			}
		}
		else {
			start[i] = 0;
			count[i] = pdim_sizes;
		}
		//start[i] = 0; 		count[i] = pdim_sizes;
		//yxDim *= count[i];		}
	}
	/*if (pvar_in != NULL)
	delete3DArrayblock_Contiguous(pvar_in);
	pvar_in = create3DArrayblock_Contiguous(count[0], count[1], count[2]);*/
	if (retncval = nc_get_vara_float(ncid, pvarid, start, count, &pvar_in[0][0]))
		ERR(retncval);
	//current time value
	if (retncval = nc_get_var1_float(ncid, ptid, &start[tIndx], &tcorvar))
		ERR(retncval);
	//close netcdf file			
	if (retncval = nc_close(ncid))
		ERR(retncval);
	return 0;
}

