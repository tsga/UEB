//copied (and modified) from snowdv.f90  
#include "uebpgdecls.h"
#include "Eigen/Dense"
using namespace Eigen;

//this uses class forcing arrays
/*__host__ __device__ void uebCell::runUEB(int dimLen2)
{
	//dimLen2
	runUEB();
}*/
/*__host__ __device__*/ int uebCell::getforcOffset(int ioffst, int dimLen2)
{
	if (ioffst == 1)
		return uebCellY*dimLen2*numTotalTimeSteps + uebCellX*numTotalTimeSteps;
	else
		return 0;
}
/*__host__ __device__*/ void uebCell::setForcingAndRadiationParamterization()
{
	//int indx = blockIdx.x*blockDim.x + threadIdx.x;
	// UTC to local time conversion
	calendardate(currentModelDateTime, Year, Month, Day, dHour);
	UTCHour = dHour - UTCOffset;
	OHour = UTCHour + lon / 15.0;
	UTCJulDat = julian(Year, Month, Day, OHour);
	calendardate(UTCJulDat, MYear, MMonth, MDay, MHour);
	fHour = (float)MHour;
	fModeldt = (float)modelDT;

	//copy data for observed for da ----7.19.16 till better way found
	int istep = 0;
	
	//      Map from wrapper input variables to UEB variables				
	P = PrecArr ;                               // / 24000;   #12.19.14 --Daymet prcp in mm/day
	V = WindspArr ;
	Ta = TempArr ;
	Tmax = TamaxArr ;
	Tmin = TaminArr ;
	Trange = Tmax - Tmin;
	//std::cout<<Trange<<std::endl;
	if (Trange <= 0)
	{
		if (snowdgtvariteflag == 1)
		{
			std::cout << "Input Diurnal temperature range is less than or equal to 0 which is unrealistic " << std::endl;
			std::cout << "Diurnal temperature range is assumed as 8 degree celsius on " << std::endl;
			std::cout << Year << " " << Month << " " << Day << std::endl;
		}
		Trange = 8.0;
	}
	//		 Flag to control radiation (irad)
	//!     0 is no measurements - radiation estimated from diurnal temperature range
	//!     1 is incoming shortwave radiation read from file (measured), incoming longwave estimated
	//!     2 is incoming shortwave and longwave radiation read from file (measured)
	//!     3 is net radiation read from file (measured)
	switch (irad)
	{
	case 0:
		Qsiobs = infrContArr[8].infdefValue;
		Qli = infrContArr[9].infdefValue;
		Qnetob = infrContArr[10].infdefValue;
		break;
	case 1:
		Qsiobs = SradArr ;                         // *3.6; // Daymet srad in W/m^2
		Qli = infrContArr[9].infdefValue;
		Qnetob = infrContArr[10].infdefValue;
		break;
	case 2:
		Qsiobs = SradArr ;                         // *3.6; // Daymet srad in W/m^2
		Qli = LradArr ;
		Qnetob = infrContArr[10].infdefValue;
		break;
	case 3:
		Qsiobs = infrContArr[8].infdefValue;
		Qli = infrContArr[9].infdefValue;
		Qnetob = NradArr ;
		break;
	default:
		std::cout << " The radiation flag is not the right number; must be between 0 and 3" << std::endl;
		getchar();
		break;
	}
	//atm. pressure from netcdf 10.30.13    //this needs revision 		//####TBC_6.20.13
	if (infrContArr[7].infType == 2)
		sitev[1] = infrContArr[7].infdefValue;
	else
		sitev[1] = ApresArr ;
	//this needs revision 		//####TBC_6.20.13
	if (infrContArr[11].infType == 2)
		Qg = infrContArr[11].infdefValue;
	else
		Qg = QgArr ;
	//!     Flag to control albedo (ireadalb)  				 
	if (infrContArr[12].infType == 2)
		Snowalb = infrContArr[12].infdefValue;
	else
		Snowalb = SnowalbArr ;
	//12.18.14 Vapor pressure of air
	if (infrContArr[6].infType == 2)
		Vp = infrContArr[6].infdefValue;
	else
		Vp = VpArr ;
	//relative humidity computed or read from file
	//#12.18.14 needs revision
	if (infrContArr[5].infType == 2)
	{
		RH = infrContArr[5].infdefValue;
	}
	else if (infrContArr[5].infType == -1)          //RH computed internally 
	{
		float eSat = 611 * exp(17.27*Ta / (Ta + 237.3)); //Pa
		RH = Vp / eSat;
	}
	else
		RH = RhArr ;
	if (RH > 1)
	{
		//std::cout<<"relative humidity >= 1 at time step "<<istep<<std::endl;
		RH = 0.99;
	}

	//  Below is code from point UEB 
	sitev[2] = Qg;
	Inpt[0] = Ta;
	Inpt[1] = P;
	Inpt[2] = V;
	Inpt[3] = RH;
	Inpt[6] = Qnetob;

	//Radiation Input Parameterization  
	hyri(MYear, MMonth, MDay, fHour, fModeldt, slope, azi, lat, HRI, cosZen);
	Inpt[7] = cosZen;
	if (irad <= 2)
	{
		atf(atff, Trange, Month, dtbar, bca, bcc);
		// We found that Model reanalysis and dowscaled data may produce some unreasonably negative solar radiation. this is simply bad data and it is generally better policy to try to give a model good data. 
		// If this is not possible, then the UEB checks will avoid the model doing anything too bad, it handles negative solar radiation in following way:
		// "no data in radiation would be to switch to the temperature method just for time steps when radiation is negative." 

		if (irad == 0 || Qsiobs < 0)     //  For cases where input is strictly negative we calculate QSI from HRI and air temp range.  This covers the case of missing data being flagged with negative number, i.e. -9999.                 
		{
			Inpt[4] = atff* Io *HRI;
			cloud(as, bs, atff, cf);   // For cloudiness fraction
		}
		else   // Here incoming solar is input
		{
			//      Need to call HYRI for horizontal surface to perform horizontal measurement adjustment
			hyri(MYear, MMonth, MDay, fHour, fModeldt, 0.0, azi, lat, HRI0, cosZen);
			//      If HRI0 is 0 the sun should have set so QSIOBS should be 0.  If it is
			//      not it indicates a potential measurement problem. i.e. moonshine
			if (HRI0 > 0)
			{
				//std::cout<<Qsiobs;
				atfimplied = findMin(Qsiobs / (HRI0*Io), 0.9); // To avoid unreasonably large radiation when HRI0 is small
				Inpt[4] = atfimplied * HRI * Io;
			}
			else
			{
				Inpt[4] = Qsiobs;
				if (Qsiobs != 0)
				{
					if (radwarnflag < 3)   //leave this warning only three times--enough to alert to non- -ve night time solar rad
					{
						std::cout << "Warning: you have nonzero nightime incident radiation of " << Qsiobs << std::endl;
						std::cout << "at date " << Year << "   " << Month << "   " << Day << "     " << dHour << std::endl;
						++radwarnflag;
					}
				}
			}
			cloud(as, bs, atff, cf);   // For cloudiness fraction  This is more theoretically correct
		}
		if (irad < 2)
		{
			qlif(Ta, RH, T_k, SB_c, Ema, Eacl, cf, QLif);
			Inpt[5] = QLif;
		}
		else
		{
			Ema = -9999;  //  These values are not evaluated but may need to be written out so are assigned for completeness
			Eacl = -9999;
			Inpt[5] = Qli;
		}
		iradfl = 0;
	}   // Long wave or shortwave either measured and calculated
	else
	{
		iradfl = 1;                    // This case is when given IRAD =3 (From Net Radiation)  
		Inpt[6] = Qnetob;
	}

	//      set control flags
	iflag[0] = iradfl;   // radiation [0=radiation is shortwave in col 5 and longwave in col 6, else = net radiation in column 7]
	//  In the code above radiation inputs were either computed or read from input files
	iflag[1] = 0;        // no 0 [/yes 1] printing
	//iflag[2] = outFile;        // Output unit to which to print
	if (ireadalb == 0)
		iflag[3] = 1;        // Albedo Calculation [a value 1 means albedo is calculated, otherwise statev[3] is albedo
	else
	{
		iflag[3] = 0;
		statev[2] = Snowalb;
	}
	//added 9.16.13
	iflag[4] = 4;
	mtime[0] = Year;
	mtime[1] = Month;
	mtime[2] = Day;
	mtime[3] = dHour;

	return;
}
/*__host__ __device__*/void uebCell::updateSimTime()
{
	//istep++;                 //running one time step at a time 7.23.16
	UPDATEtime(Year, Month, Day, dHour, modelDT);
	currentModelDateTime = julian(Year, Month, Day, dHour);
	//copy next time step
	modelStartDate[0] = Year;
	modelStartDate[1] = Month;
	modelStartDate[2] = Day;
	modelStartHour = dHour;
	return;
}
//run ueb 'ordinarily'--no ensemble/ no data assimilation
/*__host__ __device__*/void uebCell::runUEB()
{
	//
	SNOWUEB2();     
		  
	dStorage = statev[1]-Ws1+ statev[3]-Wc1;
	errMB= cumP-cumMr-cumEs-cumEc -dStorage+cumGm - cumEg; 				
				
	OutVarValues[0] = Year;
	OutVarValues[1] = Month;
	OutVarValues[2] = Day;
	OutVarValues[3] = dHour;
	OutVarValues[4] = atff;
	OutVarValues[5] = HRI;
	OutVarValues[6] = Eacl;
	OutVarValues[7] = Ema;
	OutVarValues[8] = Inpt[7]; //cosZen
	OutVarValues[9] = Inpt[0];
	OutVarValues[10] = Inpt[1];
	OutVarValues[11] = Inpt[2];
	OutVarValues[12] = Inpt[3];
	OutVarValues[13] = Inpt[4];
	OutVarValues[14] = Inpt[5];
	OutVarValues[15] = Inpt[6];	
							
	for (int i=16;i<69;i++)
	{		   
		OutVarValues[i]  = OutArr[i-16];					
	}
	OutVarValues[69]  = errMB;

	if (snowdgt_outflag == 1 )        //if debug mode 
	{
		printf(" time step: %d\n", istep);
		for (int uit = 0; uit<3; uit++)
			printf(" %d   ", (int) OutVarValues[uit] );
		for(int uit = 3; uit< 70; uit++)
			printf(" %16.4f  ", OutVarValues[uit] );
		printf(" \n");				
	}
	
	//to keep track of states not updated
	for (int is = 0; is < 6; is++)
		stateVS0[is] = statev[is];
	//8.8.16 7th state for snow surface temp 
	stateVS0[6] = tsprevday[nstepinaDay - 1];
	stateVS0[7] = taveprevday[nstepinaDay - 1];

	return;

}

//for simulation involving da but not at current time step
/*__host__ __device__*/void uebCell::setNextStepStates(int nEns)
{
	//this ensures consistency when the next time step has da ==> ensembles
	for (int ie = 0; ie < nEns + 1; ie++)
	{
		for (int is = 0; is < 6; is++)
			stateVS[is][ie] = statev[is];
		//8.8.16 7th state for snow surface temp 
		stateVS[6][ie] = tsprevday[nstepinaDay - 1];
		stateVS[7][ie] = taveprevday[nstepinaDay - 1];
	}
	return;
}
//save only state with observation
/*__host__ __device__*/ void uebCell::runUEBEnsembles(int nEns, float** ensForcingMultiplier, int outStateIndex, float* &stateOutput)              //, float* &ensAnomaly, float &ensMean)   //, bool NormalDist)
{
	// save the input forcing
	float Ptemp[7];
	for (int iforc = 0; iforc < 7; iforc++)
		Ptemp[iforc] = Inpt[iforc];               // [inputforcIndex];		    // Inpt[1];  
	//std::cout<< "state variable from previous time step: " << std::endl;
	for (int ie = 0; ie < nEns; ie++)
	{
		//perturb forcing	
		Inpt[0] = ensForcingMultiplier[0][ie] + Ptemp[0];    //Temperature 
		for (int iforc = 1; iforc < 7; iforc++)
			Inpt[iforc] = ensForcingMultiplier[iforc][ie] * Ptemp[iforc];    //Inpt[1] = Ptemp*ep;
																			 //
		for (int is = 0; is < 6; is++)
			statev[is] = stateVS[is][ie];
		//#*$8.8.16 7th state for snow surface temp 
		tsprevday[nstepinaDay - 1] = stateVS[6][ie];
		taveprevday[nstepinaDay - 1] = stateVS[7][ie];

		cumP = 0;
		cumEs = 0;
		cumEc = 0;                // Evaporation from canopy
		cumMr = 0;                 // canopy melt not added
		cumGm = 0;             //  Cumulative glacier melt
		cumEg = 0;
		//std::cout << stateBackgroundInput[ie] << " ";
		/*std::cout << std::endl << " before UEB2, ensemble " << ie << " : ";
		for (int ir = 0; ir < 6; ir++)
		std::cout << statev[ir] << "  ";
		std::cout << std::endl << " forcing: ";
		for (int ir = 0; ir < 8; ir++)
		std::cout<<" "<<Inpt[ir];*/

		//run model
		SNOWUEB2();

		// accumulate for mass balance
		cumPV[ie] = +cumP;
		cumEsV[ie] = +cumEs;
		cumEcV[ie] = +cumEc;                // Evaporation from canopy
		cumMrV[ie] = +cumMr;                 // canopy melt not added
		cumGmV[ie] = +cumGm;             //  Cumulative glacier melt
		cumEgV[ie] = +cumEg;

		dStorageV[ie] = statev[1] - Ws1 + statev[3] - Wc1;
		errMBV[ie] = cumPV[ie] - cumMrV[ie] - cumEsV[ie] - cumEcV[ie] - dStorage + cumGmV[ie] - cumEgV[ie];

		//save ensemble member state for update	
		if (outStateIndex == 5)
			stateOutput[ie] = tsprevday[nstepinaDay - 1];
		/*else if (outStateIndex == 7)
			stateOutput[ie] = taveprevday[nstepinaDay - 1];*/
		else if (outStateIndex > 2 && outStateIndex < 5)
			stateOutput[ie] = statev[outStateIndex + 1];
		else
			stateOutput[ie] = statev[outStateIndex];

		// 10.12.16 we don't need this to save stateVS here --- we are re-running same step next in the update step

		if (snowdgt_outflag == 1)        //if debug mode 
		{
			printf(" time step: %d\n", istep);
			printf(" %d %d %d %8.4f   ", Year, Month, Day, dHour);
			printf(" %16.4f %16.4f %16.4f %16.4f", atff, HRI, Eacl, Ema);
			for (int uit = 0; uit< 8; uit++)
				printf(" %16.4f  ", Inpt[uit]);
			for (int uit = 0; uit< 53; uit++)
				printf(" %16.4f  ", OutArr[uit]);
			printf("ErrMB = %16.4f \n", errMB);
		}
	}
	//the original forc must be maintained
	for (int iforc = 0; iforc < 7; iforc++)
		Inpt[iforc] = Ptemp[iforc];
	//Inpt[inputforcIndex] = Ptemp;     //Inpt[1] = Ptemp;      

	return;
}
//radiation paramterization and other settings done before this---this function runs based on exising inputs held by the object
/*__host__ __device__*/ void uebCell::runUEBEnsembles(int nEns, float** ensForcingMultiplier, float** &stateOutput)              //, float* &ensAnomaly, float &ensMean)   //, bool NormalDist)
{
	//if (snowdgt_outflag == 1) {       //if debug mode 
	/*std::cout <<std::endl<< "forcing multipliers:  "<<std::endl;
	for (int iforc = 0; iforc < 7; iforc++) {
		for (int ie = 0; ie < nEns; ie++)
			std::cout << ensForcingMultiplier[iforc][ie] << "  ";
		std::cout << std::endl;
	}*/
	//}
	// save the input forcing
	float Ptemp[7];
	for (int iforc = 0; iforc < 7; iforc++)
		Ptemp[iforc] = Inpt[iforc];               // [inputforcIndex];		    // Inpt[1];  
	//std::cout<< "state variable from previous time step: " << std::endl;
	for (int ie = 0; ie < nEns; ie++)
	{
		//perturb forcing	
		Inpt[0] = ensForcingMultiplier[0][ie] + Ptemp[0];    //Temperature 
		for (int iforc = 1; iforc < 7; iforc++)
			Inpt[iforc] = ensForcingMultiplier[iforc][ie] * Ptemp[iforc];    //Inpt[1] = Ptemp*ep;
		//
		for (int is = 0; is < 6; is++)
			statev[is] = stateVS[is][ie];
		//#*$8.8.16 7th state for snow surface temp 
		tsprevday[nstepinaDay - 1] = stateVS[6][ie];
		taveprevday[nstepinaDay - 1] = stateVS[7][ie];
		
		cumP = 0;
		cumEs = 0;
		cumEc = 0;                // Evaporation from canopy
		cumMr = 0;                 // canopy melt not added
		cumGm = 0;             //  Cumulative glacier melt
		cumEg = 0;
		//std::cout << stateBackgroundInput[ie] << " ";
		/*std::cout << std::endl << " before UEB2, ensemble " << ie << " : ";
		for (int ir = 0; ir < 6; ir++)
			std::cout << statev[ir] << "  ";
		std::cout << std::endl << " forcing: ";
		for (int ir = 0; ir < 8; ir++)
		std::cout<<" "<<Inpt[ir];*/

		//run model
		SNOWUEB2();

		/*std::cout << std::endl << " after UEB2, ensemble " << ie << " : ";
		for (int ir = 0; ir < 6; ir++)
		std::cout << statev[ir] <<"  ";
		std::cout << std::endl << " forcing: ";
		for (int ir = 0; ir < 8; ir++)
		std::cout<<" "<<Inpt[ir];*/

		/*10.12.16 no need for these states to save --they are getting updated
		for (int is = 0; is < 6; is++)
			stateVS[is][ie] = statev[is];
		//8.8.16 7th state for snow surface temp 
		stateVS[6][ie] = tsprevday[nstepinaDay - 1];
		stateVS[7][ie] = taveprevday[nstepinaDay - 1];*/
		

		// accumulate for mass balance
		cumPV[ie] = +cumP;
		cumEsV[ie] = +cumEs;
		cumEcV[ie] = +cumEc;                // Evaporation from canopy
		cumMrV[ie] = +cumMr;                 // canopy melt not added
		cumGmV[ie] = +cumGm;             //  Cumulative glacier melt
		cumEgV[ie] = +cumEg;

		dStorageV[ie] = statev[1] - Ws1 + statev[3] - Wc1;
		errMBV[ie] = cumPV[ie] - cumMrV[ie] - cumEsV[ie] - cumEcV[ie] - dStorage + cumGmV[ie] - cumEgV[ie];

		//save ensemble member state for update 
		for (int is = 0; is < 3; is++)
			stateOutput[is][ie] = statev[is];
		for (int is = 3; is < 5; is++)
			stateOutput[is][ie] = statev[is + 1];
		stateOutput[5][ie] = tsprevday[nstepinaDay - 1];
		//stateOutput[7][ie] = taveprevday[nstepinaDay - 1];

		if (snowdgt_outflag == 1)        //if debug mode 
		{
			printf(" time step: %d\n", istep);
			printf(" %d %d %d %8.4f   ", Year, Month, Day, dHour);
			printf(" %16.4f %16.4f %16.4f %16.4f", atff, HRI, Eacl, Ema);
			for (int uit = 0; uit< 8; uit++)
				printf(" %16.4f  ", Inpt[uit]);
			for (int uit = 0; uit< 53; uit++)
				printf(" %16.4f  ", OutArr[uit]);
			printf("ErrMB = %16.4f \n", errMB);
		}
	}
	//the original forc must be maintained
	for (int iforc = 0; iforc < 7; iforc++)
		Inpt[iforc] = Ptemp[iforc];
	//Inpt[inputforcIndex] = Ptemp;     //Inpt[1] = Ptemp;     

	/*10.12.16 no need for the states that are getting updated
	//compute and save ensemble mean
	for (int id = 0; id < 8; id++) {
		stateVS[id][nEns] = 0.0;
		for (int ie = 0; ie < nEns; ie++) {
			stateVS[id][nEns] += stateVS[id][ie];
		}
		stateVS[id][nEns] /= nEns;
	}*/
	
	return;
}
//set background states for next time step to filter updated
/*__host__ __device__*/
void uebCell::updateBackgroundStates(int nEns, float** updateStateArr)
{
	//std::cout <<std::endl<< "update state variable: " <<uebStates[is]<< std::endl;	
	
	for (int ie = 0; ie < nEns + 1; ie++) {
		for (int is = 0; is < 3; is++)
			stateVS[is][ie] = updateStateArr[is][ie];
		for (int is = 3; is < 6; is++)
			stateVS[is + 1][ie] = updateStateArr[is][ie];

		stateVS[3][ie] = stateVS0[3];
		stateVS[7][ie] = stateVS0[7];
		//if (snowdgt_outflag == 1)
			//std::cout << updateStateArr[ie] << "  ";
	}

	//this ensures consistency when the next time step has no da--no ensemble
	//#*$8.8.16 7th state for snow surface temp
	for (int is = 0; is < 3; is++)
		statev[is] = updateStateArr[is][nEns];
	for (int is = 3; is < 5; is++)
		statev[is + 1] = updateStateArr[is][nEns];
	tsprevday[nstepinaDay - 1] = updateStateArr[5][nEns];
	//taveprevday[nstepinaDay - 1] = updateStateArr[7][nEns];

	//also update the outputs
	OutVarValues[16] = updateStateArr[0][nEns]; //Us	
	OutVarValues[17] = updateStateArr[1][nEns]; //Ws		
	OutVarValues[18] = updateStateArr[2][nEns]; //tausn			    
	//OutVarValues[56] = updateStateArr[3][nEns]; //Wc     
	OutVarValues[36] = updateStateArr[3][nEns]; //refD
	OutVarValues[37] = updateStateArr[4][nEns]; //totRefD
	OutVarValues[30] = updateStateArr[5][nEns]; //Tsurfs
	//OutVarValues[29] = updateStateArr[7][nEns]; //Tave	
	
	return;
}