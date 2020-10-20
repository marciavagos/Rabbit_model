#include <math.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include "RAMM.h"

using namespace std;

//General Protocol Parameters (both single cell and strand)
bool hideOutput, outputBinaryData;
double dt, ft, startsave, skip;
double* bcls;
int update;
int nrOfBCLs;
string bclsString;

//bool runS1S2;
double *bcl_s2s;
int numS2s, nrOfS2BCLs;	
string bcl2sString;
int cellType;
int levelOfElectroDetail;

Settings *settings;
Constants *dataEpi, *dataEndo, *dataEpi_AF, *dataEndo_AF;

double *RyRPo_times, *RyRPo_levels, *Cao_times, *Cao_levels, *NCX_clamp_times, *NCX_clamp_levels, *RyR_State_times, **RyR_State_levels;
string s_RyRPo_Indices, s_RyRPo_times, s_RyRPo_levels, s_Cao_times, s_Cao_levels, s_NCX_clamp_times, s_NCX_clamp_levels, s_RyR_State_times, s_RyR_State_levels, *s_RyR_State_Indices;

int nrOfRyRPoTimes, RyRPo_index, nrOfCaoTimes, Cao_index, nrOfNCXClampTimes, NCX_clamp_index, nrOfRyRStateTimes, RyR_State_index;

int levelOfSignalingDetail, nrOfCDIBlockTimes;

//Strand Output
bool runStrandSims;
bool outputCai, outputIKs, outputICaL, outputICaT, outputIKr, outputEndoMidEpi;
bool outputAPD;
bool autoStopAtSS;
double ssStopThreshold;
bool loadX0, updateX0;
bool useSignalingPathway;
string x0filename;
int tissueSizeX, tissueSizeY, tissueOutputResolution;

int Cell_ij_To_Index(int x_val, int y_val)
{
	// Zero based system, i_val (= y_val) ranges from 0 ... N_CaRU_y - 1, j_val from 0 ... N_CaRU_x
	return y_val * tissueSizeX + x_val;
}

void setConstantsForAllCellTypes(double* epivar, double* endovar, double* epivar_AF, double* endovar_AF, double val)
{
	*epivar = val;	
	*endovar = val;
	*epivar_AF = val;	
	*endovar_AF = val;
}

void parseKeyValuePair(string key, string value)
{
	if(GenPurpose::sCompI(key,"hideOutput")) hideOutput = GenPurpose::parseBoolVal(value);
	if(GenPurpose::sCompI(key,"BCL")) 
	{
		nrOfBCLs = GenPurpose::countOccurences(value, ",") + 1;
		bcls = GenPurpose::parseArray(value);
		bclsString = value;
	}
	if(GenPurpose::sCompI(key,"FREQ")) settings->freq = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ISTIM")) settings->Istim = atof(value.c_str());
	if(GenPurpose::sCompI(key,"stimdur")) settings->stimdur = atof(value.c_str());
	if(GenPurpose::sCompI(key,"SKIP")) skip = atof(value.c_str());
	if(GenPurpose::sCompI(key,"NUMS2S") || GenPurpose::sCompI(key,"NUM_S2S")) numS2s = (int)atof(value.c_str());
	if(GenPurpose::sCompI(key,"DT_MAX")) settings->PARAM_dTmax = atof(value.c_str());
	if(GenPurpose::sCompI(key,"DT_MIN")) settings->PARAM_dTmin = atof(value.c_str());
	if(GenPurpose::sCompI(key,"FT")) ft = atof(value.c_str());
	if(GenPurpose::sCompI(key,"UPDATE")) update = (int)atof(value.c_str());
	if(GenPurpose::sCompI(key,"StartSave")) startsave = atof(value.c_str());
	if(GenPurpose::sCompI(key,"Params_ICaL")) settings->Params_ICaL = GenPurpose::parseArray(value);
	if(GenPurpose::sCompI(key,"Na_o")) settings->Na_o = atof(value.c_str());
	if(GenPurpose::sCompI(key,"Ca_o_times") || GenPurpose::sCompI(key,"Caotimes"))
	{
		nrOfCaoTimes = GenPurpose::countOccurences(value, ",") + 1;
		Cao_times = GenPurpose::parseArray(value);
		s_Cao_times = value;
	}
	if(GenPurpose::sCompI(key,"Ca_o"))
	{
		Cao_levels = GenPurpose::parseArray(value);
		s_Cao_levels = value;
		settings->Ca_o = Cao_levels[0];
	}
	if(GenPurpose::sCompI(key,"K_o")) settings->K_o = atof(value.c_str());
	if(GenPurpose::sCompI(key,"Cl_o")) settings->Cl_o = atof(value.c_str());
	if(GenPurpose::sCompI(key,"Cl_i")) settings->Cl_i = atof(value.c_str());
	if(GenPurpose::sCompI(key,"Flec") || GenPurpose::sCompI(key,"Flecainide")) settings->Flec = atof(value.c_str());
	if(GenPurpose::sCompI(key,"EGTA")) setConstantsForAllCellTypes(&(dataEpi->EGTA), &(dataEndo->EGTA), &(dataEpi_AF->EGTA), &(dataEndo_AF->EGTA), atof(value.c_str()));
	if(GenPurpose::sCompI(key,"Fluo3"))
	{	
		setConstantsForAllCellTypes(&(dataEpi->EGTA), &(dataEndo->EGTA), &(dataEpi_AF->EGTA), &(dataEndo_AF->EGTA), atof(value.c_str()));
		setConstantsForAllCellTypes(&(dataEpi->kon_EGTA), &(dataEndo->kon_EGTA), &(dataEpi_AF->kon_EGTA), &(dataEndo_AF->kon_EGTA), 100);
		setConstantsForAllCellTypes(&(dataEpi->koff_EGTA), &(dataEndo->koff_EGTA), &(dataEpi_AF->koff_EGTA), &(dataEndo_AF->koff_EGTA), 0.864 * 1E-3 * dataEpi->kon_EGTA);
	}	
	if(GenPurpose::sCompI(key,"fvagal")) settings->fvagal = atof(value.c_str());
	if(GenPurpose::sCompI(key,"Ach")) settings->Ach = atof(value.c_str());
	if(GenPurpose::sCompI(key,"Clamp_K_i") || GenPurpose::sCompI(key,"Clamp_Ki")) settings->clamp_K_i = atof(value.c_str());
	if(GenPurpose::sCompI(key,"Clamp_Na_i") || GenPurpose::sCompI(key,"Clamp_Nai")) settings->clamp_Na_i = atof(value.c_str());
	if(GenPurpose::sCompI(key,"Clamp_Ca_i") || GenPurpose::sCompI(key,"Clamp_Cai")) settings->clamp_Ca_i = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ICaLB")) settings->ICaLB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ICaTB")) settings->ICaTB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ICaPB")) settings->ICaPB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ISKB")) settings->ISKB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IK1B")) settings->IK1B = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IKsB")) settings->IKsB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IKrB")) settings->IKrB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IupB")) settings->IupB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"INaB")) settings->INaB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"INaKB")) settings->INaKB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"INaCaB")) settings->INaCaB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ISKB")) settings->ISKB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ITo1B")) settings->ITo1B = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IleakB")) settings->IleakB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IrelB")) settings->IrelB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IpCaB")) settings->IpCaB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ICabB")) settings->ICabB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"INabB")) settings->INabB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IKpB")) settings->IKpB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IClbB")) settings->IClbB = atof(value.c_str());
        if(GenPurpose::sCompI(key,"IClCaB")) settings->IClCaB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IsusB")) settings->IsusB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IKAchB")) settings->IKAchB = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IbarNCX")) dataEpi->IbarNCX = atof(value.c_str());
	if(GenPurpose::sCompI(key,"INa_junc_scl")) settings->INa_junc_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"INa_sl_scl")) settings->INa_sl_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ITo_scl")) settings->ITo_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IK1_scl")) settings->IK1_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IKr_scl")) settings->IKr_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IKs_scl")) settings->IKs_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"Isus_scl")) settings->Isus_scl = atof(value.c_str());
        if(GenPurpose::sCompI(key,"ISK_scl")) settings->ISK_scl = atof(value.c_str());
        if(GenPurpose::sCompI(key,"IKAch_scl")) settings->IKAch_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ICaL_junc_scl")) settings->ICaL_junc_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ICaL_sl_scl")) settings->ICaL_sl_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ICaT_junc_scl")) settings->ICaT_junc_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ICaT_sl_scl")) settings->ICaT_sl_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"INaK_junc_scl")) settings->INaK_junc_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"INaK_sl_scl")) settings->INaK_sl_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"INaCa_junc_scl")) settings->INaCa_junc_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"INaCa_sl_scl")) settings->INaCa_sl_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ICaP_sl_scl")) settings->ICaP_sl_scl = atof(value.c_str());
        if(GenPurpose::sCompI(key,"ICaP_junc_scl")) settings->ICaP_junc_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ICab_junc_scl")) settings->ICab_junc_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"ICab_sl_scl")) settings->ICab_sl_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"INab_junc_scl")) settings->INab_junc_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"INab_sl_scl")) settings->INab_sl_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"IClb_scl")) settings->IClb_scl = atof(value.c_str());
        if(GenPurpose::sCompI(key,"IClCa_scl")) settings->IClCa_scl = atof(value.c_str());
	if(GenPurpose::sCompI(key,"APDRepLevel")) settings->APDRepLevel = atof(value.c_str());
	if(GenPurpose::sCompI(key,"Ito_model")) settings->Ito_model = atof(value.c_str());
        if(GenPurpose::sCompI(key,"NCX_model")) settings->NCX_model = atof(value.c_str());
	if(GenPurpose::sCompI(key,"Ca_Buff_Factor")) settings->Ca_Buff_Factor = atof(value.c_str());

	if(GenPurpose::sCompI(key,"simIKAch") || GenPurpose::sCompI(key,"sim_IKAch")) settings->simIKAch = GenPurpose::parseBoolVal(value);
	if(GenPurpose::sCompI(key,"clampIK1NaDependence") || GenPurpose::sCompI(key,"clamp_IK1_Na_Dependence")) settings->clampIK1NaDependence = atof(value.c_str());
	if(GenPurpose::sCompI(key,"clampIKAChNaDependence") || GenPurpose::sCompI(key,"clamp_IKACh_Na_Dependence")) settings->clampIKAChNaDependence = atof(value.c_str());
	if(GenPurpose::sCompI(key,"simIKAChNaDependence") || GenPurpose::sCompI(key,"simIKAChNaDep")) settings->simIKAChNaDep = (int)atof(value.c_str());
	if(GenPurpose::sCompI(key,"BCL_S2"))
	{
		nrOfS2BCLs = GenPurpose::countOccurences(value, ",") + 1;
		bcl_s2s = GenPurpose::parseArray(value);
		bcl2sString = value;
	}
	if(GenPurpose::sCompI(key,"cellType"))
	{
		cellType = -1;
		if(GenPurpose::sCompI(value,"Epi")) cellType = 3;
		if(GenPurpose::sCompI(value,"Mid")) cellType = 2;
		if(GenPurpose::sCompI(value,"Endo")) cellType = 1;
		if(cellType < 0) cellType = (int)atof(value.c_str());
	}
	if(GenPurpose::sCompI(key,"Output_Binary_Data") || GenPurpose::sCompI(key,"OutputBinaryData") || GenPurpose::sCompI(key,"OutputDataBinary")) outputBinaryData = GenPurpose::parseBoolVal(value);
	if(GenPurpose::sCompI(key,"Output_Cai")) outputCai = GenPurpose::parseBoolVal(value);
	if(GenPurpose::sCompI(key,"Output_IKs")) outputIKs = GenPurpose::parseBoolVal(value);
	if(GenPurpose::sCompI(key,"Output_ICaL")) outputICaL = GenPurpose::parseBoolVal(value);
	if(GenPurpose::sCompI(key,"Output_ICaT")) outputICaT = GenPurpose::parseBoolVal(value);
	if(GenPurpose::sCompI(key,"Output_IKr")) outputIKr = GenPurpose::parseBoolVal(value);
	if(GenPurpose::sCompI(key,"Output_EndoMidEpi")) outputEndoMidEpi = GenPurpose::parseBoolVal(value);
	if(GenPurpose::sCompI(key,"Output_APD")) outputAPD = GenPurpose::parseBoolVal(value);
	if(GenPurpose::sCompI(key,"levelOfElectroDetail")) levelOfElectroDetail = (int)atof(value.c_str());
	if(GenPurpose::sCompI(key,"X0_File")) x0filename = value;
	if(GenPurpose::sCompI(key,"Load_X0")) loadX0 = GenPurpose::parseBoolVal(value);
	if(GenPurpose::sCompI(key,"Update_X0")) updateX0 = GenPurpose::parseBoolVal(value);
	if(GenPurpose::sCompI(key,"DivisionEndoMid")) settings->strandDivisionEndoMid = (int)atof(value.c_str());
	if(GenPurpose::sCompI(key,"DivisionMidEpi" )) settings->strandDivisionMidEpi = (int)atof(value.c_str());
	if(GenPurpose::sCompI(key,"LowConductionStart")) settings->strandLowConductionStart = (int)atof(value.c_str());
	if(GenPurpose::sCompI(key,"LowConductionEnd")) settings->strandLowConductionEnd = (int)atof(value.c_str());
	if(GenPurpose::sCompI(key,"LowConductionFactor")) settings->strandLowConductionFactor = atof(value.c_str());
	if(GenPurpose::sCompI(key,"GapJunction")) settings->Rgap = atof(value.c_str());
	if(GenPurpose::sCompI(key,"PacingLocation")) settings->PacingLocation = atof(value.c_str());
	if(GenPurpose::sCompI(key,"StrandLength")) 
	{
		settings->strandLength = (int)atof(value.c_str());
		tissueSizeX = settings->strandLength;
		tissueSizeY = 1;
		settings->ECGEnd = settings->strandLength - 15;
		delete[] settings->stimlocation;
		delete[] settings->stimlocationAddStim;
		settings->stimlocation = new bool[settings->strandLength];
		settings->stimlocationAddStim = new bool[settings->strandLength];
		for(int ab=0; ab<settings->strandLength; ab++)
		{
			settings->stimlocation[ab] = false;
			settings->stimlocationAddStim[ab] = false;
		}
		//settings->stimlocation[0] = true;
		//settings->stimlocationAddStim[0] = true;
                settings->stimlocation[settings->PacingLocation] = true;
                settings->stimlocationAddStim[settings->PacingLocation] = true;
	}
	if(GenPurpose::sCompI(key,"RyRPo_Times") || GenPurpose::sCompI(key,"RyRPoTimes"))
	{
		nrOfRyRPoTimes = GenPurpose::countOccurences(value, ",") + 1;
		RyRPo_times = GenPurpose::parseArray(value);
		s_RyRPo_times = value;
	}
	if(GenPurpose::sCompI(key,"RyRPo_levels") || GenPurpose::sCompI(key,"RyRPolevels"))
	{
		RyRPo_levels = GenPurpose::parseArray(value);
		s_RyRPo_levels = value;
	}
	if(GenPurpose::sCompI(key,"RyR_State_Times") || GenPurpose::sCompI(key,"RyRState_Times") || GenPurpose::sCompI(key,"RyRStateTimes"))
	{
		nrOfRyRStateTimes = GenPurpose::countOccurences(value, ",") + 1;
		RyR_State_times = GenPurpose::parseArray(value);
		s_RyR_State_times = value;
	}
	if(GenPurpose::sCompI(key,"RyR_State_Levels") || GenPurpose::sCompI(key,"RyRState_levels") || GenPurpose::sCompI(key,"RyRStatelevels"))
	{
		RyR_State_levels = GenPurpose::parseMatrix(value);
		s_RyR_State_levels = value;
	}
	if(GenPurpose::sCompI(key,"RyR_State_Indices") || GenPurpose::sCompI(key,"RyRState_Indices") || GenPurpose::sCompI(key,"RyRStateIndices"))
	{
		s_RyR_State_Indices = GenPurpose::parseStringArray(value, ";");
	}
	if(GenPurpose::sCompI(key,"RyRPo_Indices") || GenPurpose::sCompI(key,"RyRPoIndices"))
	{
		settings->numRyR_Po_Indices = GenPurpose::countOccurences(value, ",") + 1;
		settings->RyR_Po_Indices = GenPurpose::parseIntArray(value);
		s_RyRPo_Indices = value;
	}
	if(GenPurpose::sCompI(key,"NCX_Clamp_Times") || GenPurpose::sCompI(key,"NCXClampTimes"))
	{
		nrOfNCXClampTimes = GenPurpose::countOccurences(value, ",") + 1;
		NCX_clamp_times = GenPurpose::parseArray(value);
		s_NCX_clamp_times = value;
	}
	if(GenPurpose::sCompI(key,"NCX_Clamp_levels") || GenPurpose::sCompI(key,"NCXClampLevels"))
	{
		NCX_clamp_levels = GenPurpose::parseArray(value);
		s_NCX_clamp_levels = value;
	}
	if(GenPurpose::sCompI(key,"Adaptive_dt_Alpha")) settings->PARAM_Alpha = atof(value.c_str());	
	if(GenPurpose::sCompI(key,"applyVoltageClamp")) settings->applyVoltageClamp = GenPurpose::parseBoolVal(value);
	if(GenPurpose::sCompI(key,"VC_Values")) settings->VC_Values = GenPurpose::parseArray(value);
	if(GenPurpose::sCompI(key,"VC_Times")) settings->VC_Times = GenPurpose::parseArray(value);	
	if(GenPurpose::sCompI(key,"RyR_NP_Params")) { delete[] settings->RyRParams_NP; settings->RyRParams_NP = GenPurpose::parseArray(value); }
	if(GenPurpose::sCompI(key,"RyRParams_NP_5")) settings->RyRParams_NP[5] = atof(value.c_str());
        if(GenPurpose::sCompI(key,"RyRParams_NP_11")) settings->RyRParams_NP[11] = atof(value.c_str());

	if(GenPurpose::sCompI(key,"numRyRs") || GenPurpose::sCompI(key,"numIRelChannels")) settings->numRyRs = (int)atof(value.c_str());
	if(GenPurpose::sCompI(key,"RyRModelType") || GenPurpose::sCompI(key,"IRelModelType"))
	{			
		//if(value == "1" || value == "HH") settings->RyRModelType = 1;
		if(value == "2" || value == "Det" || value == "Markov" || GenPurpose::sCompI(value,"Markov_Det")) settings->RyRModelType = 2;				
		if(value == "3" || value == "Stoch" || GenPurpose::sCompI(value,"Markov_Stoch")) settings->RyRModelType = 3;
	}
	if(GenPurpose::sCompI(key,"N_Segments") || GenPurpose::sCompI(key,"NSegments"))
	{
		settings->N_Segments = atof(value.c_str());
	}
	if(GenPurpose::sCompI(key,"N_CaDomains") || GenPurpose::sCompI(key,"NCaDomains"))
	{
		settings->N_CaDomains = atof(value.c_str());
	}
	if(GenPurpose::sCompI(key,"segOutputResolution") || GenPurpose::sCompI(key,"segOutputRes")) settings->segOutputRes = atof(value.c_str());
	if(GenPurpose::sCompI(key,"domOutputResolution") || GenPurpose::sCompI(key,"domOutputRes")) settings->domOutputRes = atof(value.c_str());
	if(GenPurpose::sCompI(key,"tau_diff_Cai_Segment") || GenPurpose::sCompI(key,"taudiffCaiSegment")) setConstantsForAllCellTypes(&(dataEpi->tau_diff_Cai_Segment), &(dataEndo->tau_diff_Cai_Segment), &(dataEpi_AF->tau_diff_Cai_Segment), &(dataEndo_AF->tau_diff_Cai_Segment), atof(value.c_str()));
	if(GenPurpose::sCompI(key,"tau_diff_Cai_Domain") || GenPurpose::sCompI(key,"taudiffCaiDomain")) setConstantsForAllCellTypes(&(dataEpi->tau_diff_Cai_Domain), &(dataEndo->tau_diff_Cai_Domain), &(dataEpi_AF->tau_diff_Cai_Domain), &(dataEndo_AF->tau_diff_Cai_Domain), atof(value.c_str()));
	if(GenPurpose::sCompI(key,"tau_diff_Casr_Segment") || GenPurpose::sCompI(key,"taudiffCasrSegment")) setConstantsForAllCellTypes(&(dataEpi->tau_diff_Casr_Segment), &(dataEndo->tau_diff_Casr_Segment), &(dataEpi_AF->tau_diff_Casr_Segment), &(dataEndo_AF->tau_diff_Casr_Segment), atof(value.c_str()));
	if(GenPurpose::sCompI(key,"tau_diff_Casr_Domain") || GenPurpose::sCompI(key,"taudiffCasrDomain")) setConstantsForAllCellTypes(&(dataEpi->tau_diff_Casr_Domain), &(dataEndo->tau_diff_Casr_Domain), &(dataEpi_AF->tau_diff_Casr_Domain), &(dataEndo_AF->tau_diff_Casr_Domain), atof(value.c_str()));
	if(GenPurpose::sCompI(key,"tau_diff_Casrs_Segment") || GenPurpose::sCompI(key,"taudiffCasrsSegment")) setConstantsForAllCellTypes(&(dataEpi->tau_diff_Casrs_Segment), &(dataEndo->tau_diff_Casrs_Segment), &(dataEpi_AF->tau_diff_Casrs_Segment), &(dataEndo_AF->tau_diff_Casrs_Segment), atof(value.c_str()));
	if(GenPurpose::sCompI(key,"tau_diff_Casrs_Domain") || GenPurpose::sCompI(key,"taudiffCasrsDomain")) setConstantsForAllCellTypes(&(dataEpi->tau_diff_Casrs_Domain), &(dataEndo->tau_diff_Casrs_Domain), &(dataEpi_AF->tau_diff_Casrs_Domain), &(dataEndo_AF->tau_diff_Casrs_Domain), atof(value.c_str()));
	if(GenPurpose::sCompI(key,"tau_diff_Casl_Segment") || GenPurpose::sCompI(key,"taudiffCaslSegment")) setConstantsForAllCellTypes(&(dataEpi->tau_diff_Casl_Segment), &(dataEndo->tau_diff_Casl_Segment), &(dataEpi_AF->tau_diff_Casl_Segment), &(dataEndo_AF->tau_diff_Casl_Segment), atof(value.c_str()));
	if(GenPurpose::sCompI(key,"tau_diff_srs_i") || GenPurpose::sCompI(key,"taudiffsrsi")) setConstantsForAllCellTypes(&(dataEpi->tau_diff_SRS_i), &(dataEndo->tau_diff_SRS_i), &(dataEpi_AF->tau_diff_SRS_i), &(dataEndo_AF->tau_diff_SRS_i), atof(value.c_str()));
	if(GenPurpose::sCompI(key,"fact_Bmax_SLlowj") || GenPurpose::sCompI(key,"factBmaxSLlowj")) setConstantsForAllCellTypes(&(dataEpi->Bmax_SLlowj), &(dataEndo->Bmax_SLlowj), &(dataEpi_AF->Bmax_SLlowj), &(dataEndo_AF->Bmax_SLlowj), (4.6e-3*dataEpi->vmyo/dataEpi->vjunc*0.1) * atof(value.c_str()));
	if(GenPurpose::sCompI(key,"fact_Bmax_SLhighj") || GenPurpose::sCompI(key,"factBmaxSLhighj")) setConstantsForAllCellTypes(&(dataEpi->Bmax_SLhighj), &(dataEndo->Bmax_SLhighj), &(dataEpi_AF->Bmax_SLhighj), &(dataEndo_AF->Bmax_SLhighj), (1.65e-3*dataEpi->vmyo/dataEpi->vjunc*0.1) * atof(value.c_str()));
	if(GenPurpose::sCompI(key,"SERCA_Kmf_NP") || GenPurpose::sCompI(key,"SERCAKmfNP")) setConstantsForAllCellTypes(&(dataEpi->Kmf_NP), &(dataEndo->Kmf_NP), &(dataEpi_AF->Kmf_NP), &(dataEndo_AF->Kmf_NP), atof(value.c_str()));
	if(GenPurpose::sCompI(key,"SERCA_Kmr") || GenPurpose::sCompI(key,"SERCAKmr")) setConstantsForAllCellTypes(&(dataEpi->Kmr), &(dataEndo->Kmr), &(dataEpi_AF->Kmr), &(dataEndo_AF->Kmr), atof(value.c_str()));
	if(GenPurpose::sCompI(key,"Vmax_SRCaP") || GenPurpose::sCompI(key,"VmaxSRCaP")) setConstantsForAllCellTypes(&(dataEpi->Vmax_SRCaP), &(dataEndo->Vmax_SRCaP), &(dataEpi_AF->Vmax_SRCaP), &(dataEndo_AF->Vmax_SRCaP), atof(value.c_str()));
	if(GenPurpose::sCompI(key,"TissueSize") || GenPurpose::sCompI(key,"Tissue_Size"))
	{
		double* vals = GenPurpose::parseArray(value);
		tissueSizeX = vals[0]; tissueSizeY = vals[1];
		settings->strandLength = tissueSizeX * tissueSizeY;
		delete[] vals;
		settings->ECGEnd = settings->strandLength - 15;
		delete[] settings->stimlocation;
		delete[] settings->stimlocationAddStim;
		settings->stimlocation = new bool[settings->strandLength];
		settings->stimlocationAddStim = new bool[settings->strandLength];
		for(int ab=0; ab<settings->strandLength; ab++)
		{
			settings->stimlocation[ab] = false;
			settings->stimlocationAddStim[ab] = false;
		}
		for(int ab=0; ab<tissueSizeY; ab++)
		{
			settings->stimlocation[Cell_ij_To_Index(0,ab)] = true;
			settings->stimlocationAddStim[Cell_ij_To_Index(0,ab)] = true;
		}
	}
	if(GenPurpose::sCompI(key,"tissue_Output_Resolution") || GenPurpose::sCompI(key,"tissueOutputResolution")) tissueOutputResolution = (int)atof(value.c_str());
}

void parseSettingsFile(int argc, char *argv[], string fname)
{
	string::size_type posBeginIdx, posEndIdx;    	
	string            sLine, sValue;    
	string            sKeyWord;    
	const string      sDelim( "=" );    
	int numcommands = 0;

	ifstream myInputFile(fname.c_str(), ios::in);    
	
	if( !myInputFile )    
	{       
		cout << "File " << fname << " could not be opened\n";       
		return; 
	}
	else
	{
		cout << "File " << fname << " successfully opened\n";       
	}    
	
	while( getline(myInputFile,sLine) )    
	{       
		if( sLine.empty() );   // Ignore empty lines       
		else       
		{	
			//cout << "Line = [" << sLine << "]" << endl;
			string scomment = sLine.substr(0, 2);
			GenPurpose::trim(scomment);

			if(scomment != "//")
			{
				if(sLine.substr(0,8) == "#include")
				{
					parseSettingsFile(argc, argv, sLine.substr(9, sLine.length() - 9));
				}
				else
				{
					posEndIdx = sLine.find_first_of( sDelim );          
					sKeyWord  = sLine.substr( 0, posEndIdx ); // Extract word
					GenPurpose::trim(sKeyWord);
					posBeginIdx = posEndIdx + 1;  // Beginning of next word (after ':')       
					sValue    = sLine.substr( posBeginIdx, sLine.length() - posBeginIdx );						
					GenPurpose::trim(sValue);			

					if(sValue == "#ASK")
					{
						cout << "User input value for key [" << sKeyWord << "] = ";
						cin >> sValue;
					}
					else if(sValue.substr(0,1) == "#")
					{				
						string temp = sValue.substr(1,sValue.length()-1);				
						int nr = (int)atof(temp.c_str());				
						sValue = argv[nr+3];
					}

					parseKeyValuePair(sKeyWord, sValue);
					numcommands++;
				}
			}
			else
			{
				if(!hideOutput) cout << "Ignoring: [" << sLine << "]" << endl;
			}
		}    
	}

	if(!hideOutput) cout << "Processed " << numcommands << " commands" << endl;    	
}

double checkForTimeDependentSettingChanges(double t)
{
	bool resetDt = false;
	
	if(RyRPo_index < nrOfRyRPoTimes && t >= RyRPo_times[RyRPo_index])
	{				
		settings->clamp_RyR_Po = RyRPo_levels[RyRPo_index];
		cout << "RyR Po clamp now equals" << settings->clamp_RyR_Po << endl;
		RyRPo_index++;
		resetDt = true;
	}
	
	if(NCX_clamp_index < nrOfNCXClampTimes && t >= NCX_clamp_times[NCX_clamp_index])
	{				
		settings->clamp_NCX = NCX_clamp_levels[NCX_clamp_index];
		cout << "NCX clamp now equals" << settings->clamp_NCX << " pA/pF" << endl;
		NCX_clamp_index++;
		resetDt = true;
	}
	
	if(Cao_index < nrOfCaoTimes && t >= Cao_times[Cao_index])
	{				
		settings->Ca_o = Cao_levels[Cao_index];
		cout << "Ca_o now equals" << settings->Ca_o << endl;
		Cao_index++;
		resetDt = true;
	}
	
	if(RyR_State_index < nrOfRyRStateTimes && t >= RyR_State_times[RyR_State_index])
	{
		if(settings->set_RyR_State == NULL)
		{
			cout << "Creating matrix for RyR State (" << settings->N_Segments * settings->N_CaDomains << ",4)" << endl;
			settings->set_RyR_State = new double*[settings->N_Segments * settings->N_CaDomains];
			
			for(int abc=0; abc<(settings->N_Segments * settings->N_CaDomains); abc++)
			{
				settings->set_RyR_State[abc] = new double[4];
				settings->set_RyR_State[abc][0] = -1;
				settings->set_RyR_State[abc][1] = -1;
				settings->set_RyR_State[abc][2] = -1;
				settings->set_RyR_State[abc][3] = -1;
			}
		}
		
		cout << "Setting RyR State for " << s_RyR_State_Indices[RyR_State_index] << " at t = " << t << " to [" << RyR_State_levels[RyR_State_index][0] << "," << RyR_State_levels[RyR_State_index][1] << "," << RyR_State_levels[RyR_State_index][2] << "," << RyR_State_levels[RyR_State_index][3] << "]" << endl;
		string s_expanded = GenPurpose::expandArrayRanges(s_RyR_State_Indices[RyR_State_index]);
		//cout << "Expanded indices = [" << s_expanded << "]" << endl;
		int* numIndices = GenPurpose::parseIntArray(s_expanded);
		int countOfnumIndices = GenPurpose::countOccurences(s_expanded, ",") + 1;
		//cout << "# Expanded indices = [" << countOfnumIndices << "]" << endl;
		
		for(int abc=0; abc<countOfnumIndices; abc++)
		{
			int curIndex = numIndices[abc];
			//cout << " curIndex = " << curIndex << endl;
			settings->set_RyR_State[curIndex][0] = RyR_State_levels[RyR_State_index][0];
			settings->set_RyR_State[curIndex][1] = RyR_State_levels[RyR_State_index][1];
			settings->set_RyR_State[curIndex][2] = RyR_State_levels[RyR_State_index][2];
			settings->set_RyR_State[curIndex][3] = RyR_State_levels[RyR_State_index][3];
			//cout << "    settings->set_RyR_State[" << curIndex << "] = [" << settings->set_RyR_State[numIndices[abc]][0] << "," << settings->set_RyR_State[numIndices[abc]][1] << "," << settings->set_RyR_State[numIndices[abc]][2] << "," << settings->set_RyR_State[numIndices[abc]][3] << "]" << endl;
		}
		
		delete[] numIndices;
		
		RyR_State_index++;
		resetDt = true;
	}

	return resetDt;
}

int singlecellsims(int argc, char *argv[])
{
	
	//Default initial settings for single cell simulations
	tissueSizeX = 1; tissueSizeY = 1; tissueOutputResolution = 1;
	bool S1S2flag = false;	
	hideOutput = false; outputBinaryData = false;
	levelOfSignalingDetail = 1; levelOfElectroDetail = 2;
	numS2s= 0;	cellType = 2;
	bcl2sString = "";
	update = 7500000;
	double t = 0;
	skip = 0.5; dt = 0.005;	ft = 750000; startsave = 740000;	
	
	nrOfRyRPoTimes = 0; nrOfCaoTimes = 0; nrOfNCXClampTimes = 0;
	outputAPD = false; updateX0 = false; autoStopAtSS = false;
	ssStopThreshold = 0.0005;

	settings = new Settings();
		settings->bcl = 1000;
		settings->freq = 15;
	
	dataEndo = new Constants(true, false, true);
	dataEpi = new Constants(false, false, true);
	dataEndo_AF = new Constants(true, true, true);
	dataEpi_AF = new Constants(false, true, true);
	
	//dataSignaling = new BARConstants();	

	parseSettingsFile(argc, argv, argv[1]);
	
		dt = settings->PARAM_dTmin;
	//dataSignaling->init();

	if(bcl2sString == "")
	{
		numS2s = 0;
		bcl_s2s = new double[nrOfBCLs];
		for(int i=0; i<nrOfBCLs; i++) bcl_s2s[i] = 0;
	}

	RAMM *cells = new RAMM[settings->strandLength];		
	//BARsignaling *pathways = new BARsignaling[settings->strandLength];

	//BARsignaling *signaling = &pathways[0];	
	RAMM *cell = &cells[0];
	long i=0;

	double* listOfAPDs = NULL, *listOfCaTs = NULL, *listOfVmins = NULL, *listOfAPAs = NULL;
	if(autoStopAtSS) 
	{
		listOfAPDs = new double[30];	
		listOfCaTs = new double[30];
		listOfVmins = new double[30];
		listOfAPAs = new double[30];
	}
	bool autoStopInitiated = false;

	if(!hideOutput)
	{
	cout << "==================================================" << endl;
	cout << "==          Rabbit Atrial Myocyte Model         ==" << endl;
	cout << "=================================================="<< endl;
	cout << "==        Version 28.01.2019 - Marcia Vagos         ==" << endl;
	cout << "=================================================="<< endl;
	cout << "== General Parameters:                         ==" << endl;
	cout << "==   dt (adap. a / b)      = " << dt << " (" << settings->PARAM_Alpha << ", " << settings->PARAM_Beta << ")" << endl;
	cout << "==   skip                  = " << skip << endl;
	cout << "==   ft                    = " << ft << endl;
	cout << "==   startsave             = " << startsave << endl;
	if(settings->useFixedDiastolePacing)  cout << "==   diastolic interval    = " << bclsString << endl;
	if(!settings->useFixedDiastolePacing) cout << "==   bcl                   = " << bclsString << endl;
	cout << "==   S2 cycle length       = " << bcl2sString << endl;
	cout << "==   Number of S2s         = " << numS2s << endl;
	cout << "==   K_o                   = " << settings->K_o << endl;
	cout << "==   Na_o                  = " << settings->Na_o << endl;
	cout << "==   Ca_o                  = " << settings->Ca_o << endl;
	cout << "==   Cl_o                  = " << settings->Cl_o << endl;
	if(settings->clamp_Na_i >= 0) cout << "==   Na_i                  = " << settings->clamp_Na_i << endl;
	if(settings->clamp_K_i >= 0) cout << "==   K_i                   = " << settings->clamp_K_i << endl;
	if(settings->clamp_Ca_i >= 0) cout << "==   Ca_i                  = " << settings->clamp_Ca_i << endl;
	if(settings->Cl_i != 15) cout << "==   Cl_i                  = " << settings->Cl_i << endl;
	cout << "==   Output Electr / Sign  = " << levelOfElectroDetail << " / " << levelOfSignalingDetail << endl;
	cout << "==   Istim (duration)      = " << settings->Istim << " pA/pF (" << settings->stimdur << " ms)" << endl;
	cout << "==   Voltage Clamp         = " << settings->applyVoltageClamp << endl;
	cout << "==   Load / Update X0      = " << loadX0 << " / " << updateX0 << endl;
	
	if(runStrandSims)
	{
			cout << "==================================================" << endl;
			cout << "== Strand Parameters:                           ==" << endl;
			cout << "==   Tissue Size           = " << tissueSizeX << " x " << tissueSizeY << endl;
			cout << "==   Tissue Output         = " << "res: " << tissueOutputResolution << ", bin: " << outputBinaryData << endl;
			cout << "==   Res. Gap Junction     = " << settings->Rgap << endl;
			cout << "==   Division Endo Mid     = " << settings->strandDivisionEndoMid << endl;
			cout << "==   Division Mid Epi      = " << settings->strandDivisionMidEpi << endl;
			cout << "==   Low Conduction Start  = " << settings->strandLowConductionStart << endl;
			cout << "==   Low Conduction End    = " << settings->strandLowConductionEnd << endl;
			cout << "==   Low Conduction Factor = " << settings->strandLowConductionFactor << endl;
			cout << "==   ECG Start             = " << settings->ECGStart << endl;
			cout << "==   ECG End               = " << settings->ECGEnd << endl;
			cout << "==   ECG Distance          = " << settings->ECGDistance << endl;
			cout << "==   Pacing location       = " << "[";
				for(int ab=0; ab<settings->strandLength; ab++) cout << settings->stimlocation[ab] << ",";
				cout << "]" << endl;
			cout << "==   2nd Pacing location   = " << "[";
				for(int ab=0; ab<settings->strandLength; ab++) cout << settings->stimlocationAddStim[ab] << ",";
				cout << "]" << endl;			
		}
		
		cout << "==================================================" << endl;
		cout << "== Cell Parameters:                             ==" << endl;
		cout << "==   Cell Type             = " << cellType << endl;	
		if(dataEpi->kon_EGTA == 100)
		{
			cout << "==   [Fluo3]                = " << dataEpi->EGTA << endl;	
		}
		else
		{
			cout << "==   [EGTA]                = " << dataEpi->EGTA << endl;	
		}
		cout << "==   Vjunc                 = " << 1E6 * dataEpi->vjunc << endl;	
		cout << "==   Sim. IKAch (Na Dep.)  = " << settings->simIKAch << " (" << settings->simIKAChNaDep << ", clamp = " << settings->clampIKAChNaDependence << ")" << endl;
		cout << "==   Clamp IK1 Na Dep.     = " << settings->clampIK1NaDependence << endl;
		cout << "==   Use Courtemanche INa  = " << settings->useCourtemancheINa << endl;
		cout << "==   RyR Model Type (#)    = " << settings->RyRModelType << " (" << settings->numRyRs << ")" << endl;
		cout << "==   RyR Parameters        = " << "[" << settings->RyRParams_NP[0] << "," << settings->RyRParams_NP[1] << "," << settings->RyRParams_NP[2] << ", ... ," << settings->RyRParams_NP[8] << "," << settings->RyRParams_NP[9] << "," << settings->RyRParams_NP[10] << "," << settings->RyRParams_NP[11] << "]" << endl;
		cout << "==   Tau Ca Diffusion (i)  = " << dataEpi->tau_diff_Cai_Segment << " , " << dataEpi->tau_diff_Cai_Domain << endl;
		cout << "==   Tau Ca Diffusion (sr) = " << dataEpi->tau_diff_Casr_Segment << " , " << dataEpi->tau_diff_Casr_Domain << endl;
		cout << "==   Tau Ca Diffusion (sl) = " << dataEpi->tau_diff_Casl_Segment << " , " << "N/A" << endl;
		cout << "==   Tau Ca Diffusion (srs)= " << dataEpi->tau_diff_Casrs_Segment << " , " << dataEpi->tau_diff_Casrs_Domain << endl;
		cout << "==   Tau Ca Diff (srs<->i) = " << dataEpi->tau_diff_SRS_i << endl;
		cout << "==   Fact Ca Diff          = " << dataEpi->J_ca_juncsl * 1.2134e12 << " (junc<->sl) / " << dataEpi->J_ca_slmyo * 2.68510e11 << " (sl<->myo)" << endl;
		cout << "==   Cell architecture     = " << settings->N_Segments << " x " << settings->N_CaDomains << endl;
		cout << "==   Cell Output           = " << settings->segOutputRes << " x " << settings->domOutputRes << endl;
		if(nrOfRyRPoTimes > 0) cout << "==   RyR Po Times          = " << s_RyRPo_times << endl;
		if(nrOfRyRPoTimes > 0) cout << "==   RyR Po Levels         = " << s_RyRPo_levels << endl;
		if(settings->numRyR_Po_Indices > 0) cout << "==   RyR Po Indices        = " << s_RyRPo_Indices << endl;
		if(nrOfNCXClampTimes > 0) cout << "==   NCX Clamp Times       = " << s_NCX_clamp_times << endl;
		if(nrOfNCXClampTimes > 0) cout << "==   NCX Clamp Levels      = " << s_NCX_clamp_levels << endl;
		if(settings->Ito_model == 2){cout << "==   Ito_model: 	Aslanidi " << endl;}
		if(settings->Ito_model == 1){ cout << "==   Ito_model:  Lindblad " << endl;}
		if(settings->NCX_model == 1){cout << "==   NCX_model:   Lindblad " << endl;}
		else{ cout << "==   NCX_model:  Voigt&Heijman " << endl;};
		cout << "==================================================" << endl;
		cout << "== Current block: " << endl;
		if(settings->Flec != 0) 	cout << "==   Flec   :  " << settings->Flec << endl;
		if(settings->CaMKIIB != 0) 	cout << "==   CaMKII :  " << settings->CaMKIIB << endl;
		if(settings->ICaLB != 0) 	cout << "==   ICaL   :  " << settings->ICaLB << endl;
		if(settings->ICaTB != 0)        cout << "==   ICaT   :  " << settings->ICaTB << endl;
		if(settings->ICaPB != 0)        cout << "==   ICaP   :  " << settings->ICaPB << endl;
		if(settings->ICabB != 0) 	cout << "==   ICab   :  " << settings->ICabB << endl;
		if(settings->IpCaB != 0) 	cout << "==   IpCa   :  " << settings->IpCaB << endl;
		if(settings->clamp_ICaL_CDI >= 0 || settings->clamp_ICaL_VDI >= 0) 	cout << "==   ICaL   :  " << settings->clamp_ICaL_CDI << " (CDI), " << settings->clamp_ICaL_VDI << " (VDI)" << endl;
		if(settings->IKsB != 0)	cout << "==   IKs    :  " << settings->IKsB << endl;
		if(settings->IKrB != 0) 	cout << "==   IKr    :  " << settings->IKrB << endl;
		if(settings->IK1B != 0)		cout << "==   IK1    :  " << settings->IK1B << endl;
		if(settings->INaKB != 0) 	cout << "==   INaK   :  " << settings->INaKB << endl;
		if(settings->IsusB != 0)     cout << "==   Isus   :  " << settings->IsusB << endl;
		if(settings->IKpB != 0)     cout << "==   IKp    :  " << settings->IKpB << endl;
		if(settings->IupB != 0)     cout << "==   Iup    :  " << settings->IupB << endl;
		if(settings->INaB != 0)     cout << "==   INa    :  " << settings->INaB << endl;
		if(settings->ISKB != 0)    cout << "==   ISK   :  " << settings->ISKB << endl;
		if(settings->IKAchB != 0)    cout << "==   IKAch   :  " << settings->IKAchB << endl;
		if(settings->INaCaB != 0) 	cout << "==   INaCa  :  " << settings->INaCaB << endl;
		if(settings->ITo1B != 0) cout << "==   ITo1   :  " << settings->ITo1B << endl;
		if(settings->INabB != 0) cout << "==   Bckgrnd:  " << settings->INabB << " (INabB), " << settings->ICabB << " (ICabB)" << endl;
		if(settings->IClbB != 0) cout << "==   IClb   :  " << settings->IClbB << endl;
		if(settings->IClCaB != 0) cout << "==   IClCa   :  " << settings->IClCaB << endl;
		if(settings->IleakB != 0) cout << "==   Ileak  :  " << settings->IleakB << endl;
		if(settings->IrelB != 0) cout << "==   Irel  :  " << settings->IrelB << endl;

		cout << "==================================================" << endl;
                cout << "== Current scalings: " << endl;
		if(dataEpi->IbarNCX != 0) cout            << "==   IbarNCX        :  " << dataEpi->IbarNCX << endl;
		if(settings->INa_junc_scl != 0) cout   << "==   INa_junc_scl   :  " << settings->INa_junc_scl << endl;
		if(settings->INa_sl_scl != 0) cout     << "==   INa_sl_scl     :  " << settings->INa_sl_scl << endl;
		if(settings->ICaL_junc_scl != 0) cout  << "==   ICaL_junc_scl  :  " << settings->ICaL_junc_scl << endl;
		if(settings->ICaL_sl_scl != 0) cout    << "==   ICaL_sl_scl    :  " << settings->ICaL_sl_scl << endl;
		if(settings->ICaT_junc_scl != 0) cout  << "==   ICaT_junc_scl  :  " << settings->ICaT_junc_scl << endl;
		if(settings->ICaT_sl_scl != 0) cout    << "==   ICaT_sl_scl    :  " << settings->ICaT_sl_scl << endl;
		if(settings->IK1_scl != 0) cout        << "==   IK1_scl        :  " << settings->IK1_scl << endl;
		if(settings->IKr_scl != 0) cout        << "==   IKr_scl        :  " << settings->IKr_scl << endl;
		if(settings->IKs_scl != 0) cout        << "==   IKs_scl        :  " << settings->IKs_scl << endl;
		if(settings->Isus_scl != 0) cout       << "==   Isus_scl       :  " << settings->Isus_scl << endl;
		if(settings->INaK_junc_scl != 0) cout  << "==   INaK_junc_scl  :  " << settings->INaK_junc_scl << endl;
		if(settings->INaK_sl_scl != 0) cout    << "==   INaK_sl_scl    :  " << settings->INaK_sl_scl << endl;
		if(settings->INaCa_junc_scl != 0) cout << "==   INaCa_junc_scl :  " << settings->INaCa_junc_scl << endl;
		if(settings->INaCa_sl_scl != 0) cout   << "==   INaCa_sl_scl   :  " << settings->INaCa_sl_scl << endl;
		if(settings->ITo_scl != 0) cout        << "==   ITo_scl        :  " << settings->ITo_scl << endl;
                if(settings->ISK_scl != 0) cout        << "==   ISK_scl        :  " << settings->ISK_scl << endl;
                if(settings->IKAch_scl != 0) cout      << "==   IKAch_scl      :  " << settings->IKAch_scl << endl;
		if(settings->ICab_junc_scl != 0) cout  << "==   ICab_junc_scl  :  " << settings->ICab_junc_scl << endl;
		if(settings->ICab_sl_scl != 0) cout    << "==   ICab_sl_scl    :  " << settings->ICab_sl_scl << endl;
                if(settings->ICaP_junc_scl != 0) cout  << "==   ICaP_junc_scl  :  " << settings->ICaP_junc_scl << endl;
                if(settings->ICaP_sl_scl != 0) cout    << "==   ICaP_sl_scl    :  " << settings->ICaP_sl_scl << endl;
		if(settings->INab_junc_scl != 0) cout  << "==   INab_junc_scl  :  " << settings->INab_junc_scl << endl;
		if(settings->INab_sl_scl != 0) cout    << "==   INab_sl_scl    :  " << settings->INab_sl_scl << endl;
		if(settings->IClb_scl != 0) cout       << "==   IClb_scl       :  " << settings->IClb_scl << endl;
                if(settings->IClCa_scl != 0) cout      << "==   IClCa_scl      :  " << settings->IClCa_scl << endl;
		cout << "==================================================" << endl;

                cout << "== Buffering parameters: " << endl;
                if(settings->Ca_Buff_Factor != 0) cout    << "==   Ca_Buff_Factor    :  " << settings->Ca_Buff_Factor << endl;
                cout << "==================================================" << endl;

                cout << "== RyR parameters: " << endl;
                if(settings->RyRParams_NP[5] != 0) cout    << "==   RyRParams_NP5    :  " << settings->RyRParams_NP[5] << endl;
		if(settings->RyRParams_NP[11] != 0) cout    << "==   RyRParams_NP11    :  " << settings->RyRParams_NP[11] << endl;
		if(settings->numRyRs != 0) cout    << "==   numRyRs    :  " << settings->numRyRs << endl;
                cout << "==================================================" << endl;

		cout << "== SERCA parameters: " << endl;
		if(dataEpi->Kmf_NP != 0) cout    << "==   SERCA_Kmf_NP    :  " << dataEpi->Kmf_NP << endl;
		if(dataEpi->Kmr != 0) cout    << "==   SERCA_Kmr    :  " << dataEpi->Kmr << endl;
		cout << "==================================================" << endl;
	}

	char newx0filename[200];
	
	clock_t t_start = clock() / CLOCKS_PER_SEC;
	
	FILE *x0file;
	if(updateX0)
	{	
		if(loadX0)
		{
			sprintf(newx0filename, "%s_new", x0filename.c_str());
			x0file = fopen(newx0filename, "w");
		}
		else { x0file = fopen(x0filename.c_str(), "w"); }
	}	
	
	double original_ft = ft;
	for(int bclindex=0; bclindex<nrOfBCLs; bclindex++)
	{
		double* VmList = new double[settings->strandLength];
		
		if(runStrandSims)
		{			
			for(int cc=0; cc<settings->strandLength; cc++)
			{
				if(cc < settings->strandDivisionEndoMid)
				{				
					if(cellType <= 2) cells[cc].init(dataEndo, settings);
					if(cellType > 2) cells[cc].init(dataEndo_AF, settings);
				}
				else
				{
					if(cellType <= 2) cells[cc].init(dataEpi, settings);
					if(cellType > 2) cells[cc].init(dataEpi_AF, settings);
				}
				VmList[cc] = cells[cc].Vm;				
			}
		}
		else
		{
			if(cellType == 4) cell->init(dataEpi_AF, settings);
			if(cellType == 3) cell->init(dataEndo_AF, settings);
			if(cellType == 2) cell->init(dataEpi, settings);
			if(cellType == 1) cell->init(dataEndo, settings);			
		}
		
		if(loadX0)
		{
			string sLine;
			ifstream myInputFile(x0filename.c_str(), ios::in);    
		
			if(!myInputFile )    
			{       
				cout << "File " << x0filename << " could not be opened\n";       				
			}    
			else
			{
				bool stopreading = false;
				while(!stopreading)
				{
					stopreading = !(getline(myInputFile,sLine));
					int nrOfVars = GenPurpose::countOccurences(sLine, "\t");
					double* linedata = GenPurpose::parseArray(sLine, "\t");
					//cout << "sLine: " << sLine << endl;
					if(!hideOutput) cout << " X0 CL Data = " << linedata[0] << ", tot vars = " << nrOfVars << endl; 

					if(fabs(linedata[0] - (int)bcls[bclindex]) < 1E-6)
					{					
						stopreading = true;						
						int varspercell = (nrOfVars-1)/settings->strandLength;
						if(varspercell >= cells[0].lengthMinimumStateVector())
						{
							double* x0datapercell = new double[varspercell];
							for(int cc=0; cc<settings->strandLength; cc++)
							{
								for(int ii=0; ii<varspercell; ii++)	{ x0datapercell[ii] = linedata[cc*varspercell+ii+1]; }								
								cells[cc].loadX0Data(x0datapercell, varspercell);
								VmList[cc] = cells[cc].Vm;
							}
							delete[] x0datapercell;
						}
						else
						{
							if(!hideOutput) cout << "Not enough x0 variables for all cells. Using first cell to initialize strand." << endl;
							double* x0datapercell = new double[nrOfVars-1];
							for(int ii=0; ii<nrOfVars-1; ii++)	{ x0datapercell[ii] = linedata[ii+1]; }
							
							for(int cc=0; cc<settings->strandLength; cc++)
							{								
								cells[cc].loadX0Data(x0datapercell, nrOfVars-1);
								VmList[cc] = cells[cc].Vm;
							}
							delete[] x0datapercell;
						}
						
						
						if(!runStrandSims)
						{
							if(!hideOutput) cout << "Start conditions successfully loaded from file (Na_i = " << cell->segments[0].Na_i << ")" << endl;
						}
						else
						{
							if(!hideOutput) cout << "Start conditions successfully loaded from file" << endl;
							if(!hideOutput) cout << "  Endo: Vm = " << cells[settings->ECGStart].Vm << ", Cai = " << cells[settings->ECGStart].segments[0].caUnits[0].Ca_i << endl;
							if(!hideOutput) cout << "  Epi : Vm = " << cells[settings->ECGEnd].Vm << ", Cai = " << cells[settings->ECGEnd].segments[0].caUnits[0].Ca_i << endl;
						}
					}

					delete[] linedata;
				}
				myInputFile.close();											
			}
		}
				
		
		ft = original_ft;
		int beat = 0;
		RyRPo_index = 0; NCX_clamp_index = 0; Cao_index = 0;
		i = 0; t = 0;	S1S2flag = false;	
		double lastsavetime = -1, timeOfNextPacing = 0, lastPaceTime = 0;
		settings->bcl = bcls[bclindex];
		autoStopInitiated = false;
		double reducval = floor(startsave / bcls[bclindex]) * bcls[bclindex];

		if(!settings->useFixedDiastolePacing && !hideOutput) cout << "Current cycle length = " << settings->bcl << endl;
		if(settings->useFixedDiastolePacing && !hideOutput) cout << "Current diastolic interval = " << settings->bcl << ", ft = " << ft << endl;
		
		char vmfilename[200], ecgfilename[200], caifilename[200], iksfilename[200], icalfilename[200], icatfilename[200], ikrfilename[200], endofilename[200], midfilename[200], epifilename[200];
		char outputfilename[200], APDfilename[200];
		sprintf(outputfilename,"%s_%d_ms.txt", argv[2], (int)settings->bcl);
		sprintf(APDfilename,"%s_%d_ms_APD.txt", argv[2], (int)settings->bcl);
		sprintf(ecgfilename,"%s_%d_ECG.txt", argv[2], (int)(bcls[bclindex]));
		if(!outputBinaryData)
		{
			sprintf(vmfilename,"%s_%d_Vm.txt", argv[2], (int)(bcls[bclindex]));
			sprintf(caifilename, "%s_%d_Cai.txt", argv[2], (int)(bcls[bclindex]));
			sprintf(iksfilename, "%s_%d_IKs.txt", argv[2], (int)(bcls[bclindex]));
			sprintf(ikrfilename, "%s_%d_IKr.txt", argv[2], (int)(bcls[bclindex]));
			sprintf(icalfilename, "%s_%d_ICaL.txt", argv[2], (int)(bcls[bclindex]));
			sprintf(icatfilename, "%s_%d_ICaT.txt", argv[2], (int)(bcls[bclindex]));
		}
		else
		{
			sprintf(vmfilename,"%s_%d_Vm.dat", argv[2], (int)(bcls[bclindex]));
			sprintf(caifilename, "%s_%d_Cai.dat", argv[2], (int)(bcls[bclindex]));
			sprintf(iksfilename, "%s_%d_IKs.dat", argv[2], (int)(bcls[bclindex]));
			sprintf(ikrfilename, "%s_%d_IKr.dat", argv[2], (int)(bcls[bclindex]));
			sprintf(icalfilename, "%s_%d_ICaL.dat", argv[2], (int)(bcls[bclindex]));
			sprintf(icatfilename, "%s_%d_ICaT.txt", argv[2], (int)(bcls[bclindex]));
		}
		sprintf(endofilename, "%s_%d_Endo.txt", argv[2], (int)(bcls[bclindex]));
		sprintf(midfilename, "%s_%d_Mid.txt", argv[2], (int)(bcls[bclindex]));
		sprintf(epifilename, "%s_%d_Epi.txt", argv[2], (int)(bcls[bclindex]));

		FILE *vmfile; FILE *ecgfile;
		FILE* caifile; FILE* iksfile; FILE* icalfile; FILE* icatfile; FILE* ikrfile; 
		FILE* endofile; FILE* midfile; FILE* epifile;
		FILE *apdfile; FILE *outputfile;
		if(runStrandSims)
		{
			if(!outputBinaryData)
			{
				vmfile = fopen(vmfilename,"w");
				ecgfile = fopen(ecgfilename,"w");

				if(outputCai) caifile = fopen(caifilename, "w");
				if(outputIKs) iksfile = fopen(iksfilename, "w");
				if(outputIKr) ikrfile = fopen(ikrfilename, "w");
				if(outputICaL) icalfile = fopen(icalfilename, "w");
				if(outputICaT) icatfile = fopen(icatfilename, "w");
			}
			else
			{
				vmfile = fopen(vmfilename,"wb");
				ecgfile = fopen(ecgfilename,"wb");

				if(outputCai) caifile = fopen(caifilename, "wb");
				if(outputIKs) iksfile = fopen(iksfilename, "wb");
				if(outputIKr) ikrfile = fopen(ikrfilename, "wb");
				if(outputICaL) icalfile = fopen(icalfilename, "wb");
				if(outputICaT) icatfile = fopen(icatfilename, "wb");
			}
			
			if(outputEndoMidEpi)
			{
				endofile = fopen(endofilename, "w");
				midfile = fopen(midfilename, "w");
				epifile = fopen(epifilename, "w");
			}
		}
		else
		{
			if(!outputBinaryData)
			{
				outputfile = fopen(outputfilename,"w");
			}
			else
			{
				outputfile = fopen(outputfilename,"wb");
			}
		}
		
		if(outputAPD || autoStopAtSS) 
		{			
			apdfile = fopen(APDfilename, "w");
			if(!apdfile) cout << "APD file " << APDfilename << " could not be opened..." << endl;
		}

		double CaMKII_Int = 0;
		while(t<ft + (nrOfS2BCLs <= nrOfBCLs) * (numS2s >= 0) * numS2s * bcl_s2s[bclindex])			
		{		
			bool resetDt = checkForTimeDependentSettingChanges(t);
			if(resetDt)
			{
				dt = settings->PARAM_dTmin;
			}
			
			if(t >= ft && !S1S2flag)
			{
				cout << "Starting S2 cycle length: " << bcl_s2s[bclindex] << endl;
				S1S2flag = true;
				beat = 0;
				settings->bcl = bcl_s2s[bclindex];
			}
						
			if((!settings->useFixedDiastolePacing && t >= timeOfNextPacing) || (settings->useFixedDiastolePacing && cell->APD > 0 && t >= lastPaceTime + cell->APD + settings->bcl))
			{	
				lastPaceTime = t;
				beat++;
				settings->VC_index = 0;
				dt = settings->PARAM_dTmin;
				if(outputAPD || autoStopAtSS)
				{
					if(!runStrandSims)
					{
						double CaMKII_Avg = CaMKII_Int / settings->bcl;
						if(outputAPD && !hideOutput) cout << "APD: " << cell->APD << ", dVdt: " << cell->valmaxdvdt << ", CaT: " << 1E3 * (cell->CaT_Max - cell->CaT_Min) << ", SR: " << cell->Ca_sr << ", Nai = " << cell->Na_i << ", Cai = " << 1E3 * cell->Ca_i << ", Vmin: " << cell->Vmin << ", APA: " << cell->APA << endl;
							
						fprintf(apdfile, "%-4e\t%-4e\t%-4e\t%-6e\t%-6e\t%-6e\t%-6e\t%-4e\t%-4e\t%d\n", t, cell->APD, cell->valmaxdvdt, 1E3 * (cell->CaT_Max - cell->CaT_Min), cell->Ca_sr, cell->Na_i, 1E3 * cell->Ca_i, cell->Vmin, cell->APA, beat);
					}
					else
					{
						//double CV = 1000 * ((settings->ECGEnd - settings->ECGStart) * 0.01) / (cells[settings->ECGEnd].tmaxdvdt - cells[settings->ECGStart].tmaxdvdt);
						double CV = 1000 * ((1.000 * settings->ECGEnd - settings->ECGStart) / settings->strandLength * 1.29) / (cells[settings->ECGEnd].tmaxdvdt - cells[settings->ECGStart].tmaxdvdt);
						double CT = cells[settings->strandLength-1].tmaxdvdt - cells[0].tmaxdvdt;
						int midindex = (int)(0.5 * (settings->strandDivisionEndoMid + settings->strandDivisionMidEpi));
						if(midindex < 0 || midindex >= settings->strandLength) midindex = (int)(0.5 * settings->strandLength);
						
						fprintf(apdfile, "%d\t", beat);
						double minRT = 1E8, maxRT = -1E8;
						for(int x=0; x<settings->strandLength; x++)
						{
							fprintf(apdfile, "%-3e\t", cells[x].APD);
							double reptime = cells[x].tmaxdvdt + cells[x].APD;
							if(reptime < minRT) minRT = reptime;
							if(reptime > maxRT) maxRT = reptime;
						}
						fprintf(apdfile, "%-3e\t%-3e\t%-3e\n", CV, CT, maxRT - minRT);
						if(outputAPD) cout << "APDs & dVdt: ENDO = " << cells[settings->ECGStart].APD << " & " << cells[settings->ECGStart].valmaxdvdt <<  ", MID = " << cells[midindex].APD << " & " << cells[midindex].valmaxdvdt << ", EPI = " << cells[settings->ECGEnd].APD << " & " << cells[settings->ECGEnd].valmaxdvdt << ", CV = " << CV << ", CT = " << CT << ", TDR = " << (maxRT - minRT) << endl;
					}
				}								
				
				if(autoStopAtSS && !autoStopInitiated)
				{
					int a = (beat - 1) % 30;
					listOfAPDs[a] = cell->APD;
					listOfCaTs[a] = 1E3 * (cell->CaT_Max - cell->CaT_Min);
					listOfVmins[a] = cell->Vmin;
					listOfAPAs[a] = cell->APA;

					if(beat > 30)
					{
						double APD_m30 = listOfAPDs[a+1];
						double CaT_m30 = listOfCaTs[a+1];
						double Vmin_m30 = listOfVmins[a+1];
						double APA_m30 = listOfAPAs[a+1];

						if(a == 29)
						{
							APD_m30 = listOfAPDs[0];
							CaT_m30 = listOfCaTs[0];
							Vmin_m30 = listOfVmins[0];
							APA_m30 = listOfAPAs[0];
						}
						if(outputAPD) cout << "APD Change = " << fabs((APD_m30 - listOfAPDs[a]) / listOfAPDs[a]) << ", CaT Change = " << fabs((CaT_m30 - listOfCaTs[a]) / listOfCaTs[a]) << ", Vmin Change = " << fabs((Vmin_m30 - listOfVmins[a]) / listOfVmins[a]) << ", APA Change = " << fabs((APA_m30 - listOfAPAs[a]) / listOfAPAs[a]) << endl;
						if(fabs((APD_m30 - listOfAPDs[a]) / listOfAPDs[a]) < ssStopThreshold && fabs((CaT_m30 - listOfCaTs[a]) / listOfCaTs[a]) < ssStopThreshold && fabs((Vmin_m30 - listOfVmins[a]) / listOfVmins[a]) < ssStopThreshold && fabs((APA_m30 - listOfAPAs[a]) / listOfAPAs[a]) < ssStopThreshold)
						{
							t= startsave;
							beat = floor(t / settings->bcl) + 1;
							autoStopInitiated = true;
						}
						
					}
				}
				cell->APD = 0;
				cell->APA = 0;
				CaMKII_Int = 0;
			}
						
			//Update Pacing instant
			if(settings->useFixedDiastolePacing)
			{
				// When pacing with diastolic interval fixed, update timeOfNextPacing when APD is known
				if(t >= timeOfNextPacing && cell->APD > 0)
				{
					timeOfNextPacing = lastPaceTime + cell->APD + settings->bcl;
					if(timeOfNextPacing > ft) 
					{
						cout << "Setting ft to " << timeOfNextPacing << endl;
						ft = timeOfNextPacing;
					}
				}
			}
			else
			{
				if(t >= timeOfNextPacing) timeOfNextPacing += settings->bcl;
			}
					
			double teffective = t - lastPaceTime;
			
			double dt_adj = dt;
			CaMKII_Int += 0; //cell->Ca_CaMKII() * dt_adj;
						
			//Output time
			if(runStrandSims && t - lastsavetime >= skip && t >= startsave)
			{
				double toutput = t-reducval;
				fprintf(ecgfile,"%-5e\t",toutput);
				if(!outputBinaryData)
				{
					fprintf(vmfile,"%-5e\t",toutput);
					if(outputCai) fprintf(caifile,"%-5e\t",toutput);
					if(outputIKs) fprintf(iksfile,"%-5e\t",toutput);
					if(outputIKr) fprintf(ikrfile,"%-5e\t",toutput);
					if(outputICaL) fprintf(icalfile,"%-5e\t",toutput);
					if(outputICaT) fprintf(icatfile,"%-5e\t",toutput);
				}
				else
				{
					fwrite(&toutput, sizeof(double), 1, vmfile);
					if(outputCai) fwrite(&toutput, sizeof(double), 1, caifile);
					if(outputIKs) fwrite(&toutput, sizeof(double), 1, iksfile);
					if(outputIKr) fwrite(&toutput, sizeof(double), 1, ikrfile);
					if(outputICaL) fwrite(&toutput, sizeof(double), 1, icalfile);
					if(outputICaT) fwrite(&toutput, sizeof(double), 1, icatfile);
				}
				
				if(outputEndoMidEpi)
				{
					fprintf(endofile, "%-5e\t",t-reducval);
					fprintf(midfile, "%-5e\t",t-reducval);
					fprintf(epifile, "%-5e\t",t-reducval);
				}
			}
			
			//Loop through individual cells
			double ECG_Single_Far = 0, ECG_Single_Near = 0, ECG_Double_Far = 0, ECG_Double_Near = 0, ECG_Double_Alt = 0;
			double Max_State_Change = 0;			
			for(int x=0; x<settings->strandLength; x++)
			{
				VmList[x] = cells[x].Vm; // Synchronization
			}

			//Output time
			if(t - lastsavetime >= skip && t >= startsave)
			{
				if(runStrandSims)
				{
					for(int xj=0; xj<tissueSizeY; xj+=tissueOutputResolution)
					{
						for(int xi=0; xi<tissueSizeX; xi+=tissueOutputResolution)
						{
							int x = xj * tissueSizeX + xi;
							if(!outputBinaryData)
							{
								fprintf(vmfile,"%-5e\t",cells[x].Vm);
								if(outputCai) fprintf(caifile,"%-8e\t",cells[x].Ca_i);
								if(outputIKs) fprintf(iksfile,"%-8e\t",cells[x].IKs);
								if(outputIKr) fprintf(ikrfile,"%-8e\t",cells[x].IKr);
								if(outputICaL) fprintf(icalfile,"%-8e\t",(cells[x].ICaL_junc));
								if(outputICaT) fprintf(icatfile,"%-8e\t",(cells[x].ICaT_junc));
							}
							else
							{
								fwrite(&(cells[x].Vm), sizeof(double), 1, vmfile);
								if(outputCai) fwrite(&(cells[x].Ca_i), sizeof(double), 1, caifile);
								if(outputIKs) fwrite(&(cells[x].IKs), sizeof(double), 1, iksfile);
								if(outputIKr) fwrite(&(cells[x].IKr), sizeof(double), 1, ikrfile);
								if(outputICaL) fwrite(&(cells[x].ICaL_junc), sizeof(double), 1, icalfile);
								if(outputICaT) fwrite(&(cells[x].ICaT_junc), sizeof(double), 1, icatfile);
							}
						}
					}
				}
				else
				{
					if(!outputBinaryData)
					{
						fprintf(outputfile, "%-12e\t", (t - reducval));
						cell->outputData(outputfile, levelOfElectroDetail, levelOfSignalingDetail, outputBinaryData);
						fprintf(outputfile, "\n");
					}
					else
					{
						double toutput = t-reducval;
						fwrite(&(toutput), sizeof(double), 1, outputfile);
						cell->outputData(outputfile, levelOfElectroDetail, levelOfSignalingDetail, outputBinaryData);
					}
				}
			}
						
			//#pragma omp parallel for
			for(int x=0; x<settings->strandLength; x++)
			{
				int y_index_2d = (int)(x / tissueSizeX);
				int x_index_2d = x - y_index_2d * tissueSizeX;
				
				double vleft = cells[x].Vm;
				double vright = cells[x].Vm;
				double vup = cells[x].Vm;
				double vdown = cells[x].Vm;

				if(x_index_2d>0) vleft = VmList[Cell_ij_To_Index(x_index_2d-1,y_index_2d)];
				if(x_index_2d<tissueSizeX-1) vright = VmList[Cell_ij_To_Index(x_index_2d+1,y_index_2d)];
				if(y_index_2d>0) vup = VmList[Cell_ij_To_Index(x_index_2d,y_index_2d-1)];
				if(y_index_2d<tissueSizeY-1) vdown = VmList[Cell_ij_To_Index(x_index_2d,y_index_2d+1)];

				if(x >= settings->ECGStart && x <= settings->ECGEnd)
				{
					//ECG_Single_Far += (vleft-vright)/pow(2.0+(settings->strandLength-x)/100.0,2.0);
					ECG_Single_Far += (vleft-vright)/pow(2.0+(settings->strandLength-1.0 * x)/settings->strandLength * 1.29,2.0);					
					ECG_Single_Near += (vleft-vright)/pow(settings->ECGDistance+(settings->strandLength-x)/100.0,2.0);				
					ECG_Double_Far += (vleft-vright)* ( 1.0 / pow(2.0+0.01*(settings->strandLength-x),2.0) - 1.0 / pow(2.0+0.01*x,2.0));
					ECG_Double_Near += (vleft-vright)* ( (1.0 / pow(settings->ECGDistance+0.01*(settings->strandLength-x),2.0)) - (1.0 / pow(settings->ECGDistance+0.01*x,2.0)));
					ECG_Double_Alt += (vleft-vright)* ( (1.0 / pow(settings->ECGDistance+0.01*(settings->strandLength-x),2.0)) + (1.0 / pow(settings->ECGDistance+0.01*x,2.0)));
				}
								
				cells[x].update(x, teffective, dt_adj, vleft, vright, vup, vdown);
				if(cells[x].Max_State_Change > Max_State_Change) Max_State_Change = cells[x].Max_State_Change;
				
				// Something went wrong / output state variables to debug!
				if (isnan(cells[x].Vm) || cells[x].Vm < -150 || cells[x].Vm>120)
				{
					cout << "t = " << teffective << " x = " << x << "(" << x_index_2d << "," << y_index_2d << "), dt_adj = " << dt_adj << ", vleft = " << vleft << ", vright = " << vright << ", vup = " << vup << ", vdown = " << vdown << endl;
					cout << "caiont = " << cells[x].caiont << ", kiont = " << cells[x].kiont << ", naiont = " << cells[x].naiont << ", clont = " << cells[x].clont << endl;
					
					if(isnan(cells[x].caiont) || fabs(cells[x].caiont) > 300)
					{
						cout << "Vm = " << cells[x].Vm << "(was: " << VmList[x] << "), Ca_i = " << cells[x].Ca_i << ", Ca_sl = " << cells[x].Ca_sl << ", Ca_srs = " << cells[x].Ca_srs << ", ICaL_junc = " << cells[x].ICaL_junc << ", ICaL_sl = " << cells[x].ICaL_sl << ", ICaT_junc = " << cells[x].ICaT_junc << ", ICaT_sl = " << cells[x].ICaT_sl << ", INaCa_junc = " << cells[x].INaCa_junc << ", INaCa_sl = " << cells[x].INaCa_sl << ", ICab_junc = " << cells[x].ICab_junc << ", ICab_sl = " << cells[x].ICab_sl << ", ICaP_junc = " << cells[x].ICaP_junc << ", ICaP_sl = " << cells[x].ICaP_sl << endl;
						//cout << "ICaLState[0] = " << cells[x].segments[0].membrane[0].ICaLState_NP[0] << ", ICaLState[1] = " << cells[x].segments[0].membrane[0].ICaLState_NP[1] << ", ICaLState[2] = " << cells[x].segments[0].membrane[0].ICaLState_NP[2] << ", ICaLState[3] = " << cells[x].segments[0].membrane[0].ICaLState_NP[3] << endl;
						//cout << "ICaLState[4] = " << cells[x].segments[0].membrane[0].ICaLState_NP[4] << ", ICaLState[5] = " << cells[x].segments[0].membrane[0].ICaLState_NP[5] << ", ICaLState[6] = " << cells[x].segments[0].membrane[0].ICaLState_NP[6] << ", ICaLState[7] = " << cells[x].segments[0].membrane[0].ICaLState_NP[7] << endl;
					}
					
					if(isnan(cells[x].kiont) || fabs(cells[x].kiont) > 300)
					{
						cout << "Vm = " << cells[x].Vm << "(was: " << VmList[x] << "), IK1 = " << cells[x].IK1 << ", IKr = " << cells[x].IKr << ", IKs = " << cells[x].IKs << endl;
					}

					if(isnan(cells[x].naiont) || fabs(cells[x].naiont) > 300)
					{
						//cout << "Vm = " << cells[x].Vm << "(was: " << VmList[x] << "), INa = " << cells[x].INa << ", m_NP = " << cells[x].m_NP << ", H_NP = " << cells[x].H_NP << ", J_NP = " << cells[x].J_NP << endl;
					}
					
					if(isnan(cells[x].Vm))	perror("V is nan");
                    if(cells[x].Vm < -150)	{
                        cout << "INaCa_junc=" << cells[x].INaCa_junc << ", INaCa_sl=" << cells[x].INaCa_sl << endl;
                        perror("V is below -150 mV");
                    }
					GenPurpose::wait(1.0);
					
					if(runStrandSims)
					{
						//fclose(vmfile); 
						fclose(vmfile); fclose(ecgfile);
						if(outputCai) fclose(caifile);
						if(outputIKs) fclose(iksfile);
						if(outputIKr) fclose(ikrfile);
						if(outputICaL) fclose(icalfile);
						if(outputICaT) fclose(icatfile);
						if(outputEndoMidEpi)
						{
							fclose(endofile);
							fclose(midfile);
							fclose(epifile);
						}
					}
					else
					{
						fclose(outputfile);
					}
					if(outputAPD || autoStopAtSS) fclose(apdfile);
					return 1;
				}				
			} // End of Looping through individual cells
			
			if (t - lastsavetime >= skip && t >= startsave)
			{
				if(runStrandSims)
				{
					fprintf(ecgfile, "%-8e\t%-8e\t%-8e\t%-8e\t%8e\n", ECG_Single_Far, ECG_Single_Near, ECG_Double_Far, ECG_Double_Near, ECG_Double_Alt);
					if(!outputBinaryData)
					{
						fprintf(vmfile, "\n");
						if(outputCai) fprintf(caifile,"\n");
						if(outputIKs) fprintf(iksfile,"\n");
						if(outputIKr) fprintf(ikrfile,"\n");
						if(outputICaL) fprintf(icalfile,"\n");
						if(outputICaT) fprintf(icatfile,"\n");
					}
					if(outputEndoMidEpi)
					{
						int cc = settings->ECGStart;
						cells[cc].outputData(endofile, levelOfElectroDetail, 0, outputBinaryData); fprintf(endofile, "\n");						

						cc = (int)(0.5 * (settings->strandDivisionEndoMid + settings->strandDivisionMidEpi));
						if(cc < 0 || cc >= settings->strandLength) cc = (int)(0.5 * settings->strandLength);
						cells[cc].outputData(midfile, levelOfElectroDetail, 0, outputBinaryData); fprintf(midfile, "\n");

						cc = settings->ECGEnd;
						cells[cc].outputData(epifile, levelOfElectroDetail, 0, outputBinaryData); fprintf(epifile, "\n");
					}
				}
				lastsavetime = t;
			}

			t+=dt_adj;
			i++;
		
			double dt_bar = settings->PARAM_dTmax - (settings->PARAM_dTmax - settings->PARAM_dTmin) / (1 + pow(settings->PARAM_Km / Max_State_Change, settings->PARAM_h));
			double g_dtalt = settings->PARAM_Alpha + (settings->PARAM_Beta - settings->PARAM_Alpha) / (1 + exp((dt_bar - dt) / 0.0001));			
			dt = g_dtalt * dt + (1 - g_dtalt) * dt_bar;
			if((!settings->useFixedDiastolePacing && t+dt >= timeOfNextPacing) || (settings->useFixedDiastolePacing && cell->APD > 0 && t+dt >= lastPaceTime + cell->APD + settings->bcl))
			{
				dt = timeOfNextPacing - t;
				//cout << " beat = " << beat << ", t = " << t << ", dt = " << dt << endl; 
			}
			//if(S1S2flag && t+dt - ft >= beat * settings->bcl) dt = beat * settings->bcl - t + ft;
			
			if (i%update==0)
			{
				if(autoStopAtSS)
				{
					int a = (beat - 1) % 30;
					double APD_m30 = listOfAPDs[a+1];
					double CaT_m30 = listOfCaTs[a+1];
					double Vmin_m30 = listOfVmins[a+1];
					double APA_m30 = listOfAPAs[a+1];

					if(a == 29)
					{
						APD_m30 = listOfAPDs[0];
						CaT_m30 = listOfCaTs[0];
						Vmin_m30 = listOfVmins[0];
						APA_m30 = listOfAPAs[0];
					}
					cout << 100*t / ft << "% complete; (steady state measures: " << fabs((APD_m30 - listOfAPDs[a]) / listOfAPDs[a]) << " (APD), " << fabs((CaT_m30 - listOfCaTs[a]) / listOfCaTs[a]) << " (CaT) )"  << fabs((Vmin_m30 - listOfVmins[a]) / listOfVmins[a]) << " (Vmin) )" << fabs((APA_m30 - listOfAPAs[a]) / listOfAPAs[a]) << " (APA) )" << endl;
				}
				else
				{
					cout << 100*t / ft << "% complete" << endl;
				}
			}
		} //End While t < ft
		
		dt = settings->PARAM_dTmin;
		delete[] VmList;
		delete[] listOfAPDs;
		delete[] listOfCaTs;
		delete[] listOfVmins;
		delete[] listOfAPAs;

		if(updateX0)
		{			
			fprintf(x0file, "%d\t", (int)settings->bcl);
			for(int x=0; x<settings->strandLength; x++)
			{
				cells[x].outputData(x0file, 10, 10, false);
			}	
			fprintf(x0file, "\n");			
		}		

		if(runStrandSims)
		{
			//if(updateX0) fclose(x0file);
			fclose(vmfile); 
			fclose(ecgfile);
			if(outputCai) fclose(caifile);
			if(outputIKs) fclose(iksfile);
			if(outputIKr) fclose(ikrfile);
			if(outputICaL) fclose(icalfile);
			if(outputICaT) fclose(icatfile);
			if(outputEndoMidEpi)
			{
				fclose(endofile);
				fclose(midfile);
				fclose(epifile);
			}
			fclose(apdfile);
		}

		// Finalize APD and determine BVR (single cell only)
		if(!runStrandSims && (outputAPD || autoStopAtSS))
		{
			if(outputAPD)
			{
				//double CaMKII_Avg = CaMKII_Int / settings->bcl;
				cout << "Almost done..." << endl;
				cout << "APD: " << cell->APD << ", dVdt: " << cell->valmaxdvdt << ", CaT: " << 1E3 * (cell->CaT_Max - cell->CaT_Min) << ", SR: " << cell->Ca_sr << ", Nai = " << cell->Na_i << ", Cai = " << cell->Ca_i << ", Vmin = " << cell->Vmin << ", APA = " << cell->APA << endl;
				fprintf(apdfile, "%-4e\t%-4e\t%-4e\t%-6e\t%-6e\t%-6e\t%-6e\t%-4e\t%-4e\t%d\n", t, cell->APD, cell->valmaxdvdt, 1E3 * (cell->CaT_Max - cell->CaT_Min), cell->Ca_sr, cell->Na_i, cell->Ca_i, cell->Vmin, cell->APA, beat);
			}
			fclose(apdfile);
			
			GenPurpose::wait(0.1);

			ifstream myInputFile(APDfilename, ios::in);	
			if(myInputFile && beat > 30)    
			{
				string sLine;  
				int countBeats = 0;
				double prevAPD = 0, curAPD = 0;
				double* dAPD = new double[beat-4];
				double totAPD = 0;
				
				// Skip first 5 beats
				for(int abc = 0; abc < 4; abc++) getline(myInputFile,sLine);

				while( getline(myInputFile,sLine) )    
				{
					//cout << "  counter = " << countBeats << ", line = [" << sLine << "] ";
					//Determine curAPD
					double* vals = GenPurpose::parseArray(sLine, "\t");
					curAPD = vals[1];
					//cout << "curAPD = " << curAPD << endl;
					totAPD += curAPD;					
					
					if(countBeats >= 1 && countBeats <= beat) dAPD[countBeats-1] = abs(curAPD - prevAPD);
					prevAPD = curAPD;
					delete[] vals;
					countBeats++;
				}
				myInputFile.close();
				
				int offset = 30;
				double* STV = new double[countBeats - offset];
				for(int abc = 0; abc < countBeats-offset; abc++) STV[abc] = GenPurpose::sum(&(dAPD[abc]), 29) / (30 * sqrt(2));
				
				cout << "APD = " << totAPD / countBeats <<" STV = " << GenPurpose::average(STV, countBeats-offset) << "+/-" << GenPurpose::std(STV, countBeats-offset) << " (ms)" << endl;
				delete[] STV;
				delete[] dAPD;
				
			}			
			
			fclose(outputfile);
		}	

	} //End for int bclindex < nrOfBCLs
	
	if(updateX0) fclose(x0file);
	clock_t t_end = clock() / CLOCKS_PER_SEC;

	if(!hideOutput) cout << "Computations took " << (int)(t_end - t_start) << " seconds" << endl;

	cout << endl;
	
	return 0;
}

int main(int argc,char *argv[])
{
	// 3 Inputs required:
	// - SettingsFile				= argv[1]
	// - Outputname					= argv[2]
	// - Single cell vs Strand		= argv[3]
	// Note: argv[0] = process name

	if(argc < 4)
	{
		cout << "Error: 3 input arguments required\nUsage: " << argv[0] << " Settingsfile Outputname cell/strand\n";
	}

	string sType = argv[3];
	if(sType == "CELL" || sType == "cell" || sType == "1")
	{
		cout << "Starting single cell simulations...\n";
		runStrandSims = false;
		return singlecellsims(argc, argv);
	}

	if(sType == "STRAND" || sType == "strand" || sType == "2")
	{
		cout << "Starting strand simulations...\n";
		runStrandSims = true;
		return singlecellsims(argc, argv);
		//return strandsims(argc, argv);
	}

	return(0);
}
