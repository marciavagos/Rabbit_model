class Settings
{
	public:
		int N_Segments, N_CaDomains;
		int segOutputRes, domOutputRes;
		double PARAM_dTmin, PARAM_dTmax, PARAM_h, PARAM_Km, PARAM_Alpha, PARAM_Beta;
		
		double bcl, LastBCL;
		int freq, CT, PM;
		//int isAF, isRA;
		
		double K_o, Ca_o, Na_o, Cl_o;
		double Mg_i, Cl_i;
		double clamp_K_i, clamp_Na_i, clamp_Ca_i;
		double clamp_ICaL_CDI, clamp_ICaL_VDI;
		double clamp_RyR_Po, clamp_NCX;
		double Ca_Buff_Factor;
		double Bmax_TnClow, kon_tncl, koff_tncl, kon_tnchca, koff_tnchca, Bmax_TnChigh, kon_tnchmg, koff_tnchmg;
		double kon_cam, koff_cam, Bmax_CaM, kon_myoca, koff_myoca, Bmax_myosin, kon_myomg, koff_myomg, kon_sr, koff_sr, Bmax_SR;
		double Buff_Ca_TnCL, Buff_Ca_CaM, Buff_Ca_Myosin_ca;
		//double Buff_Ca_TnCL, Buff_Ca_TnCHc, Buff_Ca_TnCHm, Buff_Ca_CaM, Buff_Ca_Myosin_ca, Buff_Ca_Myosin_mg, Buff_Ca_SRB;
		int* RyR_Po_Indices;
		double **set_RyR_State;
		int numRyR_Po_Indices;
		
		double Istim, stimdur;
		double IaddStim, addStimEnd, addStimStart;
		double* addStimParams;
		bool* stimlocation; bool* stimlocationAddStim;
		int PacingLocation;
		bool useFixedDiastolePacing; // When used, bcl defines the fixed diastolic period

		double* RyRParams_NP;
		double* Params_ICaL;
		
		bool simIKAch;
		double clampIK1NaDependence, clampIKAChNaDependence;
		int simIKAChNaDep;
		bool useCourtemancheINa;
		
		double IKsB, IKrB, INaKB, INaCaB, IKpB, IsusB, IK1B, ICaLB, ICaTB, ICaPB, INabB, ITo1B, INaB, IClbB, IpCaB,  ICabB, IrelB, IupB, IKAchB, IClCaB;
		double ISKB, CaMKIIB, CTKClB, CTNaClB, IleakB, ItrB, IdiffB;

		double INa_junc_scl, INa_sl_scl, ICaL_junc_scl, ICaL_sl_scl, ICaT_junc_scl, ICaT_sl_scl, ICaP_junc_scl, ICaP_sl_scl;
		double IK1_scl, Isus_scl, ITo_scl, IKs_scl, IKr_scl, ISK_scl, IKAch_scl, IClb_scl, IClCa_scl;
		double INaCa_junc_scl, INaCa_sl_scl, INaK_junc_scl, INaK_sl_scl, INab_junc_scl, INab_sl_scl, ICab_junc_scl, ICab_sl_scl;

		double Flec;

		double fvagal, Ach;
		
		int strandLength, varsPerCell, strandDivisionEndoMid, strandDivisionMidEpi, strandLowConductionStart, strandLowConductionEnd, ECGStart, ECGEnd;		
		double strandLowConductionFactor, Rgap, CVAnisotropyFactor;
		double ECGDistance;

		double APDRepLevel;

		bool applyVoltageClamp;
		double* VC_Values;
		double* VC_Times;
		int VC_index;
		
		int numRyRs, RyRModelType;

		int Ito_model, NCX_model;

		Settings()
		{
			N_Segments = 50; N_CaDomains = 18;
			segOutputRes = 2; domOutputRes = 1;
			PARAM_dTmin = 0.005; PARAM_dTmax = 0.2; PARAM_h = 4; PARAM_Km = 0.05; PARAM_Alpha = 1.0; PARAM_Beta = 0.1;
			APDRepLevel = 0.90;

			PM = 0; CT = 0;
			
			bcl = 1000;	LastBCL = bcl;
			freq = 10;
	
			K_o = 5.4; Ca_o = 1.8; Na_o = 140; Cl_o = 132;
			Mg_i = 1.0; Cl_i = 30.0;
			clamp_K_i = -1;
			clamp_Na_i = -1;
			clamp_Ca_i = -1;
			clamp_ICaL_CDI = -1;
			clamp_ICaL_VDI = -1;
			clamp_RyR_Po = -1; clamp_NCX = -10;
			numRyR_Po_Indices = -1;

			Ca_Buff_Factor = 1.0;
                        Bmax_TnClow = 70e-3;    // [mM]                      % TnC low affinity
                        koff_tncl = 5*19.6E-3; //(1+0.5*ISO)*19.6e-3;    % [1/ms]
                        kon_tncl = 5*32.7;
                        Bmax_TnChigh = 140e-3;  // [mM]                      % TnC high affinity
                        koff_tnchca = 0.032e-3;
                        kon_tnchca = 2.37;
                        koff_tnchmg = 3.33e-3;
                        kon_tnchmg = 3e-3;
                        Bmax_CaM = 24e-3;       // [mM] **? about setting to 0 in c-code**   % CaM buffering
                        koff_cam = 238e-3;
                        kon_cam = 34;
                        Bmax_myosin = 140e-3;   // [mM]                      % Myosin buffering
                        koff_myoca = 0.46e-3;
                        kon_myoca = 13.8;
                        koff_myomg = 0.057e-3;
                        kon_myomg = 0.0157;
                        Bmax_SR = 19*.9e-3;     // [mM] (Bers text says 47e-3) 19e-3
                        koff_sr = 60e-3;
                        kon_sr = 100;
			Buff_Ca_TnCL = 1.8065e-02;
			//Buff_Ca_TnCHc = 1.2774e-1;
			//Buff_Ca_TnCHm = 5.7527e-3;
			Buff_Ca_CaM = 6.9090e-4;
			Buff_Ca_Myosin_ca = 3.9052e-3;
			//Buff_Ca_Myosin_mg = 1.3558e-1;
			//Buff_Ca_SRB = 4.3974e-3;

			// IKAch
			fvagal = 0; //[per second]
			Ach = 0.0001; //[millimolar]

			set_RyR_State = NULL; //[0] = -1; set_RyR_State[1] = -1; set_RyR_State[2] = -1; set_RyR_State[3] = -1;
			
			Istim = -12.5; IaddStim = 0; stimdur = 5.0;
			addStimStart = 0; addStimEnd = 0;
			addStimParams = 0; useFixedDiastolePacing = false;
	
			strandLength = 1; varsPerCell = 77;
			strandDivisionEndoMid = 35;
			strandDivisionMidEpi = 100;
			strandLowConductionStart = 95;
			strandLowConductionEnd = 105;
			strandLowConductionFactor = 5.0;
			Rgap = 2.6;
			CVAnisotropyFactor = 5.0;
			ECGStart = 16; ECGEnd = 145;
			ECGDistance = 1.5;
			
			stimlocation = new bool[1];
				stimlocation[0] = true;
		
			stimlocationAddStim = new bool[1];
				stimlocationAddStim[0] = true;

			PacingLocation = 0;

			applyVoltageClamp = false;
			//VC_Values = NULL; VC_Times = NULL;
			VC_index = 0;
			
			ICaLB = 0.0; ICaTB = 0.0; IKsB = 0.0; IKrB = 0.0; INaKB = 0.0; INaCaB = 0.0; //-0.35; 
				IKpB = 0.0; IK1B = 0.0; IsusB = 0.0;
			INabB= 0.0; ITo1B = 0.0; INaB = 0.0; IpCaB = 0.0; ICabB = 0.0; //-0.85; 
				IrelB = 0.0; IupB = 0.0; //0.5;
			CTKClB = 0.0; CTNaClB = 0.0; ItrB = 0.0; IdiffB = 0.0; CaMKIIB = 0.0; 
			IKAchB = 0.0; IleakB = 0.0; IClbB = 0.0; ISKB= 0.0; IClCaB = 0.0;

			INa_junc_scl = 1.0; INa_sl_scl = 1.0; ICaL_junc_scl = 1.0; ICaL_sl_scl = 1.0; ICaT_junc_scl = 1.0; ICaT_sl_scl = 1.0; 
			IK1_scl = 1.0; Isus_scl = 1.0; ITo_scl = 1.0; IKs_scl = 1.0; IKr_scl = 1.0; INaCa_junc_scl = 1.0; INaCa_sl_scl = 1.0; 
			INaK_junc_scl = 1.0; INaK_sl_scl = 1.0; INab_junc_scl = 1.0; INab_sl_scl = 1.0; ICab_junc_scl = 1.0; ICab_sl_scl = 1.0;
			IClb_scl = 1.0; ISK_scl = 1.0; IKAch_scl = 1.0; IClCa_scl = 1.0; 

			RyRParams_NP = new double[25];
				RyRParams_NP[0] = 0.2; RyRParams_NP[1] = 0.22; RyRParams_NP[2] = 12; RyRParams_NP[3] = 0.63;
				RyRParams_NP[4] = 0.001; RyRParams_NP[5] = 8E-005; RyRParams_NP[6] = 0.0035; RyRParams_NP[7] = 6.3;
				RyRParams_NP[8] = 0.0; RyRParams_NP[9] = 1; RyRParams_NP[10] = 1; RyRParams_NP[11] = 0.007;
				RyRParams_NP[12] = 20; RyRParams_NP[13] = 0.9; RyRParams_NP[14] = 0.0008; RyRParams_NP[15] = 0.00015;
				RyRParams_NP[16] = 0.05; RyRParams_NP[17] = 2; RyRParams_NP[18] = 0.00065; RyRParams_NP[19] = 1;
				RyRParams_NP[20] = 1.75; RyRParams_NP[21] = 1; RyRParams_NP[22] = 50; RyRParams_NP[23] = 1.5; RyRParams_NP[24] = 0.5;
			
			simIKAch = false; clampIK1NaDependence = -1; simIKAChNaDep = -1;
			clampIKAChNaDependence = -1;
			
			Params_ICaL = new double[25];
				Params_ICaL[0] = 0.59; Params_ICaL[1] = 0.8; Params_ICaL[2] = 0.052;
				Params_ICaL[3] = 13; Params_ICaL[4] = 0.132; Params_ICaL[5] = 10;
				Params_ICaL[6] = 9.45; Params_ICaL[7] = 70; Params_ICaL[8] = 65;
				Params_ICaL[9] = 6; Params_ICaL[10] = 12; Params_ICaL[11] = 0.213;
				Params_ICaL[12] = 10.807; Params_ICaL[13] = 27.5; Params_ICaL[14] = 5.0;
				Params_ICaL[15] = 0.2474; Params_ICaL[16] = 43.825; Params_ICaL[17] = 45;
				Params_ICaL[18] = 0.001; Params_ICaL[19] = 5; Params_ICaL[20] = 1E-06;
				Params_ICaL[21] = 0.00017; Params_ICaL[22] = 0.0035; Params_ICaL[23] = 0.1; Params_ICaL[24] = 0.5;
			
			useCourtemancheINa = true;
						
			numRyRs = 198000;
			RyRModelType = 2;
			
			Flec = 0;
			Ito_model = 1; //Lindblad 
			//Ito_model = 2; //Asladini
			NCX_model = 1; //Lindblad [else: Aslanidi]
		}
		
		~Settings()
		{
		}
};
