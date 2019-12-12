#include "RAMM.h"

void RAMM::init(Constants *idata, Settings *isettings)
{		
	showlog = false;
	data = idata;
	settings = isettings;	

	tmaxdvdt = 0;
	valmaxdvdt = 0;
	vrest = Vm;
	vmax = 0;
	APD = 0;
	CaT_Min = 1;
	CaT_Max = -1;
	Vmin = 1;
	APA = 0;

	Ca_junc = 7.3e-05; //2.0738e-04; //3.2231e-4;
	if(settings->clamp_Ca_i >= 0) Ca_junc = settings->clamp_Ca_i;
	
	double Km_SLL = data->koff_sll / data->kon_sll;
	double Km_SLH = data->koff_slh / data->kon_slh;
	
	double Buff_Ca_SLLj = data->Bmax_SLlowj * Ca_junc / (Ca_junc + Km_SLL);
	double Buff_Ca_SLHj = data->Bmax_SLhighj * Ca_junc / (Ca_junc + Km_SLH);
	double Buff_Ca_SLLsl = data->Bmax_SLlowsl * Ca_junc / (Ca_junc + Km_SLL);
	double Buff_Ca_SLHsl = data->Bmax_SLhighsl * Ca_junc / (Ca_junc + Km_SLH);

	Vm = -73.496; 
	//Na_i = 9.1366; Na_junc = 9.1361; Na_sl = 9.1364; 
	Na_i = settings->Nai_initial; Na_junc = settings->Nai_initial; Na_sl = settings->Nai_initial;
	//Na_i = 6.86; Na_junc = 6.86; Na_sl = 6.86;
	K_i = 120;
	Ca_i = 2.0738e-04; Ca_sl = 2.2918e-4;
	//Ca_i = 7.3e-05; Ca_sl = 7.3e-5;
	Cl_i = 30.0; Cl_o = 132.0;
	if(settings->clamp_Ca_i >= 0) Ca_i = settings->clamp_Ca_i;
	if(settings->clamp_Ca_i >= 0) Ca_sl = settings->clamp_Ca_i;
	Ca_sr = 0.25;
	Ca_srs = 7.3e-05;
	
	E_Cl = log(Cl_i/Cl_o)/data->frt;
        EB_Cl = E_Cl - 0.49 * (E_Cl + 30.59);

	Copy_Ca_sl = new double*[settings->N_Segments];
	Copy_Na_sl = new double*[settings->N_Segments];
	Copy_Ca_sr = new double*[settings->N_Segments];
	Copy_Ca_srs = new double*[settings->N_Segments];
	Copy_Ca_i = new double*[settings->N_Segments];
	Copy_Na_i = new double[settings->N_Segments];
	
	double Ca_cis_observed = Ca_srs;
	double exp_Ca_cis = 1.0 / (1 + exp(-(Ca_cis_observed - settings->RyRParams_NP[4])/settings->RyRParams_NP[5]));
	double Act_ss =  (settings->RyRParams_NP[6] * Ca_cis_observed / 1E-4 + settings->RyRParams_NP[3] * exp_Ca_cis);		
	double Inact_ss = settings->RyRParams_NP[0] + (settings->RyRParams_NP[9] - settings->RyRParams_NP[0]) / (1 + pow(Ca_sr / settings->RyRParams_NP[1], settings->RyRParams_NP[2]));
					
	segments = new Segment[settings->N_Segments];
		for(int ii=0; ii<settings->N_Segments; ii++)
		{
			Copy_Ca_sl[ii] = new double[2];
			Copy_Na_sl[ii] = new double[2];
			Copy_Ca_i[ii] = new double[settings->N_CaDomains];
			Copy_Ca_sr[ii] = new double[settings->N_CaDomains];
			Copy_Ca_srs[ii] = new double[settings->N_CaDomains];
			
			// Loop through the two membrane sides
			for(int jj=0; jj<2; jj++)
			{
				//Initializing state variables
				segments[ii].membrane[jj].m_Na = 0.01309; segments[ii].membrane[jj].h1_Na = 0.706; segments[ii].membrane[jj].h2_Na = 6.1493e-01; 
				segments[ii].membrane[jj].to_r = 0.00006; segments[ii].membrane[jj].to_s1 = 0.5753; segments[ii].membrane[jj].to_s2 = 0.39871; segments[ii].membrane[jj].to_s3 = 0.57363;
				segments[ii].membrane[jj].z_Ks = 0.02032; segments[ii].membrane[jj].pa_Kr = 0.00016; segments[ii].membrane[jj].pi_Kr = 0.76898;
				segments[ii].membrane[jj].d_L = 0.00003; segments[ii].membrane[jj].f_L = 0.99981;
                                segments[ii].membrane[jj].d_T = 0.00046; segments[ii].membrane[jj].f_T = 0.30752;
                                segments[ii].membrane[jj].O_KCa_junc = 1.0; segments[ii].membrane[jj].O_KCa_sl = 1.0;
				segments[ii].membrane[jj].a_KAch = 1.0;

				//segments[ii].membrane[jj].sus_a = 3.9041e-01; segments[ii].membrane[jj].sus_i = 9.5862e-01;
				//segments[ii].membrane[jj].mL_Na = 1.0092e-02; segments[ii].membrane[jj].hL_Na = 4.053e-02; segments[ii].membrane[jj].xr = 1.490e-3;
				//Markov formulation of ICaL only
				//segments[ii].membrane[jj].ICaLState_NP[0] = 1.0; segments[ii].membrane[jj].ICaLState_NP[1] = 0.0; segments[ii].membrane[jj].ICaLState_NP[2] = 0.0; segments[ii].membrane[jj].ICaLState_NP[3] = 0.0;
				//segments[ii].membrane[jj].ICaLState_NP[4] = 0.0; segments[ii].membrane[jj].ICaLState_NP[5] = 0.0; segments[ii].membrane[jj].ICaLState_NP[6] = 0.0; segments[ii].membrane[jj].ICaLState_NP[7] = 0.0;
				//segments[ii].membrane[jj].f_CaL_P = 9.9841e-01; segments[ii].membrane[jj].f_CaL_Ca_junc = 4.4350e-2; segments[ii].membrane[jj].f_CaL_Ca_sl = 3.2145e-2;
		
				segments[ii].membrane[jj].Buff_Na_junc = 3.6109; segments[ii].membrane[jj].Buff_Na_sl = 7.8794e-1;
				segments[ii].membrane[jj].Buff_Ca_SLLj = Buff_Ca_SLLj; segments[ii].membrane[jj].Buff_Ca_SLLsl = Buff_Ca_SLLsl; //2.1059e-02;
				segments[ii].membrane[jj].Buff_Ca_SLHj = Buff_Ca_SLHj; segments[ii].membrane[jj].Buff_Ca_SLHsl = Buff_Ca_SLHsl; //1.8888e-01;
				segments[ii].membrane[jj].Buff_Ca_EGTA_sl = 0;
				
				segments[ii].membrane[jj].Ca_junc = Ca_junc;
				segments[ii].membrane[jj].Ca_sl = Ca_sl;
				segments[ii].membrane[jj].Na_junc = Na_junc;
				segments[ii].membrane[jj].Na_sl = Na_sl;
				
				//segments[ii].membrane[jj].SK_State[0] = 1.0; segments[ii].membrane[jj].SK_State[1] = 0.0;
				//segments[ii].membrane[jj].SK_State[2] = 0.0; segments[ii].membrane[jj].SK_State[3] = 0.0;
			}
			
			// Loop through the cytosolic Ca2+ domains
			segments[ii].caUnits = new CaDomain[settings->N_CaDomains];
			for(int jj=0; jj<settings->N_CaDomains; jj++)
			{
				if(settings->RyRModelType == 3)
				{
					double fact = settings->numRyRs * 1.0 / (settings->N_Segments * settings->N_CaDomains);
					segments[ii].caUnits[jj].RyRState_NP[0] = 0; segments[ii].caUnits[jj].RyRState_NP[1] = 0; segments[ii].caUnits[jj].RyRState_NP[2] = 0; 
					segments[ii].caUnits[jj].RyRState_NP[3] = (int)(fact); segments[ii].caUnits[jj].RyRState_NP[4] = 0; segments[ii].caUnits[jj].RyRState_NP[5] = 0;
					segments[ii].caUnits[jj].RyRState_NP[6] = 0; segments[ii].caUnits[jj].RyRState_NP[7] = 0;
				}
				else
				{		
					segments[ii].caUnits[jj].RyRState_NP[0] = 0; segments[ii].caUnits[jj].RyRState_NP[1] = 0; segments[ii].caUnits[jj].RyRState_NP[2] = 0; 
					segments[ii].caUnits[jj].RyRState_NP[3] = 1.0; segments[ii].caUnits[jj].RyRState_NP[4] = 0; segments[ii].caUnits[jj].RyRState_NP[5] = 0;
					segments[ii].caUnits[jj].RyRState_NP[6] = 0; segments[ii].caUnits[jj].RyRState_NP[7] = 0;
				}
				segments[ii].caUnits[jj].RyR_Integr = 0.0;
				segments[ii].caUnits[jj].RyR_SEP = 1.0;
				
				segments[ii].caUnits[jj].Buff_Ca_TnCL = 1.8065e-02;
				if(segments[ii].caUnits[jj].Buff_Ca_TnCL > data->Bmax_TnClow) segments[ii].caUnits[jj].Buff_Ca_TnCL = data->Bmax_TnClow;
				segments[ii].caUnits[jj].Buff_Ca_TnCHc = 1.2774e-1;
				segments[ii].caUnits[jj].Buff_Ca_TnCHm = 5.7527e-3;
				segments[ii].caUnits[jj].Buff_Ca_CaM = 6.9090e-4;
				segments[ii].caUnits[jj].Buff_Ca_Myosin_ca = 3.9052e-3;
				segments[ii].caUnits[jj].Buff_Ca_Myosin_mg = 1.3558e-1;
				if(segments[ii].caUnits[jj].Buff_Ca_Myosin_ca + segments[ii].caUnits[jj].Buff_Ca_Myosin_mg > data->Bmax_myosin) 
				{
					segments[ii].caUnits[jj].Buff_Ca_Myosin_ca = 0;
					segments[ii].caUnits[jj].Buff_Ca_Myosin_mg = 0;
				}
				segments[ii].caUnits[jj].Buff_Ca_SRB = 4.3974e-3;
				segments[ii].caUnits[jj].Buff_Ca_CSQN = data->Bmax_Csqn * Ca_sr / (Ca_sr + (data->koff_csqn / data->kon_csqn)); //1.1352;
				segments[ii].caUnits[jj].Buff_Ca_EGTA_cyt = 0;
				
				segments[ii].caUnits[jj].Buff_Ca_SLLsrs = Buff_Ca_SLLj; 
				segments[ii].caUnits[jj].Buff_Ca_SLHsrs = Buff_Ca_SLHj;
		
				segments[ii].caUnits[jj].Ca_i = Ca_i;
				segments[ii].caUnits[jj].Ca_sr = Ca_sr;
				segments[ii].caUnits[jj].Ca_srs = Ca_srs;
			}
			
			// Other stuff
			segments[ii].Na_i = Na_i;
			segments[ii].K_i = K_i;
			
		}
}

double RAMM::check_State_Change(double oldval, double dval, double dt)
{
	if(fabs(dval) > Max_State_Change) 
	{
		Max_State_Change = fabs(dval);
	}
	return oldval + dval * dt;
}

void RAMM::update(double t, double dt)
{
	update(0, t, dt, Vm, Vm, Vm, Vm);
} 

void RAMM::update(int cellnumber, double t, double dt, double vleft, double vright)
{
	update(cellnumber, t, dt, vleft, vright, Vm, Vm);
}

void RAMM::update(int cellnumber, double t, double dt, double vleft, double vright, double vup, double vdown)
{
	
	if(showlog) cout << " t = " << t << ", Vm = " << Vm << ", Cai = " << Ca_i << ", dt = " << dt << endl;
	//cout << "t = " << t << "  "; 
	
	Max_State_Change = 0.0000;

	if(t <= dt)	
	{
		CaT_Min = 1; CaT_Max = -1;
		tmaxdvdt = 0; valmaxdvdt = 0;
		vrest = Vm; vmax = Vm;
		APD = 0; APA = 0;
	}

	// Allow clamping of intracellular ion concentrations when values are set in settings
	// ===================================================================
	if(settings->clamp_K_i >= 0.0) 
	{
		for(int ii=0; ii<settings->N_Segments; ii++)
		{
			segments[ii].K_i = settings->clamp_K_i;
		}
	}
	if(settings->clamp_Na_i >= 0.0)
	{
		for(int ii=0; ii<settings->N_Segments; ii++)
		{
			segments[ii].Na_i = settings->clamp_Na_i;
			segments[ii].membrane[0].Na_sl = settings->clamp_Na_i;
			segments[ii].membrane[1].Na_sl = settings->clamp_Na_i;
			segments[ii].membrane[0].Na_junc = settings->clamp_Na_i;
			segments[ii].membrane[1].Na_junc = settings->clamp_Na_i;
		}
	}
	if(settings->clamp_Ca_i >= 0.0)
	{
		for(int ii=0; ii<settings->N_Segments; ii++)
		{
			for(int jj=0; jj<settings->N_CaDomains; jj++)
			{
				segments[ii].caUnits[jj].Ca_i = settings->clamp_Ca_i; 			
				segments[ii].caUnits[jj].Ca_srs = settings->clamp_Ca_i; 			
			}
			
			segments[ii].membrane[0].Ca_sl = settings->clamp_Ca_i;
			segments[ii].membrane[1].Ca_sl = settings->clamp_Ca_i;
		}
	}

    // Determine stimulus current
	if(settings->applyVoltageClamp)
	{
		// Switch to next step in VC protocol
		if(t >= settings->VC_Times[settings->VC_index]) 
		{
			settings->VC_index++;
			Max_State_Change = 100;
		}
		Istim = settings->Istim * (t > 0 && t <= 0+settings->stimdur && settings->stimlocation[cellnumber]);
		if(settings->VC_index > 0) Istim = -50 * (0.005 / dt) * (settings->VC_Values[settings->VC_index-1] - Vm);
	}
	else
	{
		// Current clamp
		Istim = settings->Istim * (t > 0 && t <= 0+settings->stimdur && settings->stimlocation[cellnumber]);
		if(settings->addStimParams == 0)
		{
			Istim += settings->IaddStim * (t > settings->addStimStart && t <= settings->addStimEnd && settings->stimlocationAddStim[cellnumber]);
		}
	}
    
	// Make a copy of concentrations to be able to calculate diffusion
	// ===================================================================		
	for(int ii=0; ii<settings->N_Segments; ii++)
	{
		for(int jj=0; jj<settings->N_CaDomains; jj++)
		{
			Copy_Ca_i[ii][jj] = segments[ii].caUnits[jj].Ca_i;
			Copy_Ca_sr[ii][jj] = segments[ii].caUnits[jj].Ca_sr;
			Copy_Ca_srs[ii][jj] = segments[ii].caUnits[jj].Ca_srs;
		}
		
		Copy_Ca_sl[ii][0] = segments[ii].membrane[0].Ca_sl;
		Copy_Ca_sl[ii][1] = segments[ii].membrane[1].Ca_sl;
		Copy_Na_sl[ii][0] = segments[ii].membrane[0].Na_sl;
		Copy_Na_sl[ii][1] = segments[ii].membrane[1].Na_sl;
		Copy_Na_i[ii] = segments[ii].Na_i;
	}
	

	Ca_i = 0; Ca_sl = 0; Ca_junc = 0; Ca_sr = 0; Ca_srs = 0; Na_i = 0; 
	INaCa_junc = 0; INaCa_sl = 0; ICab_junc = 0; ICab_sl = 0; IpCa_junc = 0; IpCa_sl = 0; 
	ICaL_junc = 0; ICaL_sl = 0; ICaT_junc = 0; ICaT_sl = 0; ICaP_junc = 0; ICaP_sl = 0;
	Jrel = 0; Jleak = 0; Jup = 0; Buff_Ca_EGTA_cyt = 0;
	
	INa_junc = 0; INa_sl = 0; INaK_junc = 0; INaK_sl = 0; INab_junc = 0; INab_sl = 0; IKr = 0; IKs = 0; IK1 = 0; Isus = 0; ITo = 0; IClb = 0; ISK_junc = 0; ISK_sl = 0; IKAch = 0; IClCa_junc = 0; IClCa_sl = 0;
	
	// == Main update loop of all segments
	// ===================================================================	
	for(int ii=0; ii<settings->N_Segments; ii++)
	{
		for(int jj=0; jj<2; jj++)
		{
			Ion_HH_Na(ii, jj, dt); 
			//Ion_HH_NaL(ii, jj, dt);
			Ion_HH_NaK(ii, jj, dt);
			Ion_HH_Nab(ii, jj, dt);
			Ion_HH_Clb(ii, jj, dt);
			Ion_HH_K1(ii, jj, dt);
			Ion_HH_sus(ii, jj, dt);
			Ion_HH_Kr(ii, jj, dt);
			Ion_HH_To(ii, jj, dt);
			Ion_HH_Ks(ii, jj, dt);
			Ion_HH_NaCa(ii, jj, dt);	
			Ion_HH_CaL(ii, jj, dt);
			Ion_HH_CaT(ii, jj, dt);
                        Ion_HH_CaP(ii, jj, dt);
			Ion_HH_SK(ii, jj, dt);
			Ion_HH_KAch(ii, jj, dt);
			Ion_HH_ClCa(ii, jj, dt);
			Ion_Ca_Handling(ii, jj, dt);     // time-independent calcium currents
		}
	
		// == Calculate all local Ca domains
		// ===================================================================		
		for(int jj=0; jj<settings->N_CaDomains; jj++)
		{			
			double Ca_i_L_Domain = segments[ii].caUnits[jj].Ca_i;
			double Ca_sr_L_Domain = segments[ii].caUnits[jj].Ca_sr;
			double Ca_srs_L_Domain = segments[ii].caUnits[jj].Ca_srs;
			if(jj > 0) { Ca_i_L_Domain = Copy_Ca_i[ii][jj-1]; Ca_sr_L_Domain = Copy_Ca_sr[ii][jj-1]; Ca_srs_L_Domain = Copy_Ca_srs[ii][jj-1]; }
			
			double Ca_i_R_Domain = segments[ii].caUnits[jj].Ca_i;
			double Ca_sr_R_Domain = segments[ii].caUnits[jj].Ca_sr;
			double Ca_srs_R_Domain = segments[ii].caUnits[jj].Ca_srs;
			if(jj < settings->N_CaDomains - 1) { Ca_i_R_Domain = Copy_Ca_i[ii][jj+1]; Ca_sr_R_Domain = Copy_Ca_sr[ii][jj+1]; Ca_srs_R_Domain = Copy_Ca_srs[ii][jj+1];}
			
			double Ca_i_U_Segment = segments[ii].caUnits[jj].Ca_i;
			double Ca_sr_U_Segment = segments[ii].caUnits[jj].Ca_sr;
			double Ca_srs_U_Segment = segments[ii].caUnits[jj].Ca_srs;
			if(ii > 0) { Ca_i_U_Segment = Copy_Ca_i[ii-1][jj]; Ca_sr_U_Segment = Copy_Ca_sr[ii-1][jj]; Ca_srs_U_Segment = Copy_Ca_srs[ii-1][jj]; }
			
			double Ca_i_D_Segment = segments[ii].caUnits[jj].Ca_i;
			double Ca_sr_D_Segment = segments[ii].caUnits[jj].Ca_sr;
			double Ca_srs_D_Segment = segments[ii].caUnits[jj].Ca_srs;
			if(ii < settings->N_Segments - 1) { Ca_i_D_Segment = Copy_Ca_i[ii+1][jj]; Ca_sr_D_Segment = Copy_Ca_sr[ii+1][jj]; Ca_srs_D_Segment = Copy_Ca_srs[ii+1][jj]; }

			updateCaDomain(dt, t, ii, jj, Ca_i_L_Domain, Ca_i_R_Domain, Ca_i_U_Segment, Ca_i_D_Segment, Ca_sr_L_Domain, Ca_sr_R_Domain, Ca_sr_U_Segment, Ca_sr_D_Segment, Ca_srs_L_Domain, Ca_srs_R_Domain, Ca_srs_U_Segment, Ca_srs_D_Segment);
			
			Ca_i += segments[ii].caUnits[jj].Ca_i;
			Ca_sr += segments[ii].caUnits[jj].Ca_sr;
			Ca_srs += segments[ii].caUnits[jj].Ca_srs;
			Jrel += segments[ii].caUnits[jj].Jrel;
			Jleak += segments[ii].caUnits[jj].Jleak;
			Jup += segments[ii].caUnits[jj].Jup;
			Buff_Ca_EGTA_cyt += segments[ii].caUnits[jj].Buff_Ca_EGTA_cyt;
		}
		
		// == Calculate all membrane domains
		// ===================================================================
			double Na_i_U_Segment = segments[ii].Na_i;
			if(ii > 0) Na_i_U_Segment = Copy_Na_i[ii-1];
			double Na_i_D_Segment = segments[ii].Na_i;
			if(ii < settings->N_Segments - 1) Na_i_D_Segment = Copy_Na_i[ii+1];
				
			double dNa_i = data->J_na_slmyo/data->vmyo * (segments[ii].membrane[0].Na_sl - segments[ii].Na_i + segments[ii].membrane[1].Na_sl - segments[ii].Na_i) + (Na_i_U_Segment - segments[ii].Na_i + Na_i_D_Segment - segments[ii].Na_i) / data->tau_diff_Cai_Segment;
			
			for(int jj=0; jj<2; jj++)
			{				
				double Na_sl_U_Segment = segments[ii].membrane[jj].Na_sl;
				double Ca_sl_U_Segment = segments[ii].membrane[jj].Ca_sl;
				if(ii > 0) { Na_sl_U_Segment = Copy_Na_sl[ii-1][jj]; Ca_sl_U_Segment = Copy_Ca_sl[ii-1][jj]; }
				double Na_sl_D_Segment = segments[ii].membrane[jj].Na_sl;
				double Ca_sl_D_Segment = segments[ii].membrane[jj].Ca_sl;
				if(ii < settings->N_Segments - 1) { Na_sl_D_Segment = Copy_Na_sl[ii+1][jj]; Ca_sl_D_Segment = Copy_Ca_sl[ii+1][jj]; }
				double Ca_i_neighbor = Copy_Ca_i[ii][0];			
				if(jj == 1) Ca_i_neighbor = Copy_Ca_i[ii][settings->N_CaDomains-1];
				double Ca_srs_neighbor = Copy_Ca_srs[ii][0];
				if(jj == 1) Ca_srs_neighbor = Copy_Ca_srs[ii][settings->N_CaDomains-1];
				
				UpdateConcentrations(ii, jj, dt, t, Na_sl_U_Segment, Na_sl_D_Segment, Ca_sl_U_Segment, Ca_sl_D_Segment, Ca_i_neighbor, Ca_srs_neighbor);
				
				// Update whole cell currents
				INa_junc += segments[ii].membrane[jj].INa_junc; INa_sl += segments[ii].membrane[jj].INa_sl; 
				INaK_junc += segments[ii].membrane[jj].INaK_junc; INaK_sl += segments[ii].membrane[jj].INaK_sl; 
				INab_junc += segments[ii].membrane[jj].INab_junc; INab_sl += segments[ii].membrane[jj].INab_sl; 
				IKr += segments[ii].membrane[jj].IKr; 
				IKs += segments[ii].membrane[jj].IKs; 
				IK1 += segments[ii].membrane[jj].IK1; 
				//IKp += segments[ii].membrane[jj].IKp; 
				Isus += segments[ii].membrane[jj].Isus; 
				ITo += segments[ii].membrane[jj].ITo; 
				IClb += segments[ii].membrane[jj].IClb; 
				IClCa_junc += segments[ii].membrane[jj].IClCa_junc; IClCa_sl += segments[ii].membrane[jj].IClCa_sl;
				IKAch += segments[ii].membrane[jj].IKAch;
				ISK_junc += segments[ii].membrane[jj].ISK_junc;	ISK_sl += segments[ii].membrane[jj].ISK_sl;
				Ca_sl += segments[ii].membrane[jj].Ca_sl; Ca_junc += segments[ii].membrane[jj].Ca_junc;
				//ICaL_K += segments[ii].membrane[jj].ICaL_K;
				//ICaL_Na_junc += segments[ii].membrane[jj].ICaL_Na_junc; ICaL_Na_sl += segments[ii].membrane[jj].ICaL_Na_sl;
				ICaL_junc += segments[ii].membrane[jj].ICaL_junc; ICaL_sl += segments[ii].membrane[jj].ICaL_sl;
                                ICaT_junc += segments[ii].membrane[jj].ICaT_junc; ICaT_sl += segments[ii].membrane[jj].ICaT_sl;
				ICaP_junc += segments[ii].membrane[jj].ICaP_junc; ICaP_sl += segments[ii].membrane[jj].ICaP_sl;
				INaCa_junc += segments[ii].membrane[jj].INaCa_junc; INaCa_sl += segments[ii].membrane[jj].INaCa_sl;
				ICab_junc += segments[ii].membrane[jj].ICab_junc; ICab_sl += segments[ii].membrane[jj].ICab_sl;
				IpCa_junc += segments[ii].membrane[jj].IpCa_junc; IpCa_sl += segments[ii].membrane[jj].IpCa_sl;
			}
			
			segments[ii].Na_i = check_State_Change(segments[ii].Na_i, dNa_i, dt);
			Na_i += segments[ii].Na_i;
			double IK_tot = segments[ii].membrane[0].ITo + segments[ii].membrane[1].ITo + segments[ii].membrane[0].IKr + segments[ii].membrane[1].IKr + segments[ii].membrane[0].IKs + segments[ii].membrane[1].IKs + segments[ii].membrane[0].IK1 + segments[ii].membrane[1].IK1; 
			IK_tot = IK_tot - 2*(segments[ii].membrane[0].INaK_junc + segments[ii].membrane[1].INaK_junc + segments[ii].membrane[0].INaK_sl + segments[ii].membrane[1].INaK_sl) + segments[ii].membrane[0].Isus + segments[ii].membrane[1].Isus + 2*Istim;
			double dK_i = -0.5*IK_tot*data->CmemF/data->vmyo; // 0.0;
			segments[ii].K_i = check_State_Change(segments[ii].K_i, dK_i, dt);
	}		
	
	//Determine whole-cell currents (Normalization step)
	INa_junc = INa_junc / (2 * settings->N_Segments); INa_sl = INa_sl / (2 * settings->N_Segments);
	INaK_junc = INaK_junc / (2 * settings->N_Segments); INaK_sl = INaK_sl / (2 * settings->N_Segments);
	INab_junc = INab_junc / (2 * settings->N_Segments); INab_sl = INab_sl / (2 * settings->N_Segments);

	IKr = IKr / (2 * settings->N_Segments);
	IKs = IKs / (2 * settings->N_Segments);
	IK1 = IK1 / (2 * settings->N_Segments);
	//IKp = IKp / (2 * settings->N_Segments);
	Isus = Isus / (2 * settings->N_Segments);
	ITo = ITo / (2 * settings->N_Segments);
	IClb = IClb / (2 * settings->N_Segments);
	IClCa_junc = IClCa_junc / (2 * settings->N_Segments); IClCa_sl = IClCa_sl / (2 * settings->N_Segments);
	IKAch = IKAch / (2 * settings->N_Segments);
	ISK_junc = ISK_junc / (2 * settings->N_Segments); ISK_sl = ISK_sl / (2 * settings->N_Segments);
	
	Ca_junc = Ca_junc / (2 * settings->N_Segments); Ca_sl = Ca_sl / (2 * settings->N_Segments);	
	Na_junc = Na_junc / (2 * settings->N_Segments); Na_sl = Na_sl / (2 * settings->N_Segments);
	Na_i = Na_i / settings->N_Segments;
	
	//ICaL_K = ICaL_K / (2 * settings->N_Segments);
	//ICaL_Na_junc = ICaL_Na_junc / (2 * settings->N_Segments); ICaL_Na_sl = ICaL_Na_sl / (2 * settings->N_Segments);
	INaCa_junc = INaCa_junc / (2 * settings->N_Segments); INaCa_sl = INaCa_sl / (2 * settings->N_Segments);
	ICab_junc = ICab_junc / (2 * settings->N_Segments); ICab_sl = ICab_sl / (2 * settings->N_Segments);
	IpCa_junc = IpCa_junc / (2 * settings->N_Segments); IpCa_sl = IpCa_sl / (2 * settings->N_Segments);
	ICaL_junc = ICaL_junc / (2 * settings->N_Segments); ICaL_sl = ICaL_sl / (2 * settings->N_Segments);
	ICaT_junc = ICaT_junc / (2 * settings->N_Segments); ICaT_sl = ICaT_sl / (2 * settings->N_Segments);
	ICaP_junc = ICaP_junc / (2 * settings->N_Segments); ICaP_sl = ICaP_sl / (2 * settings->N_Segments);

	Ca_i = Ca_i / (settings->N_Segments * settings->N_CaDomains);
	Ca_sr = Ca_sr / (settings->N_Segments * settings->N_CaDomains);
	Ca_srs = Ca_srs / (settings->N_Segments * settings->N_CaDomains);
	Jrel = Jrel / (settings->N_Segments * settings->N_CaDomains);
	Jleak = Jleak / (settings->N_Segments * settings->N_CaDomains);
	Jup = Jup / (settings->N_Segments * settings->N_CaDomains);
	Buff_Ca_EGTA_cyt = Buff_Ca_EGTA_cyt / (settings->N_Segments * settings->N_CaDomains);
	
	// == Update Vm and global ion concentrations
	// ===================================================================	
		caiont = (ICaL_junc+ICaL_sl)+(ICaT_junc+ICaT_sl)+(ICab_junc+ICab_sl)+(ICaP_junc+ICaP_sl)-2*(INaCa_junc+INaCa_sl);
		//caiont = (ICaL_junc+ICaL_sl)+(ICaT_junc+ICaT_sl)+(ICab_junc+ICab_sl)+(ICaP_junc+ICaP_sl)-2*(INaCa_junc+INaCa_sl)+Jup-Jrel;
		naiont = (INa_junc+INa_sl)+3*(INaCa_junc+INaCa_sl)+3*(INaK_junc+INaK_sl)+(INab_junc+INab_sl);
		kiont = IKr+IKs+IK1+Isus-2*(INaK_junc+INaK_sl)+ITo+(ISK_junc+ISK_sl)+IKAch+Istim;
		clont = IClb+(IClCa_junc+IClCa_sl);
		
		double Rgap = settings->Rgap;
		if(cellnumber >= settings->strandLowConductionStart && cellnumber <= settings->strandLowConductionEnd) Rgap = Rgap * settings->strandLowConductionFactor;
		double RgapupDown = Rgap * settings->CVAnisotropyFactor;
		double loc_l = data->l / 1.0E4; // conversion to cm
		double loc_a = data->a / 1.0E4;		
		
		// New 2D, more optimized, correct units (note sqrt(PI) = 1.772453, RCG = 2, dx = loc_l, dy = 2 * loc_a)
		double dvdt = 1000 * ((vleft - 2*Vm + vright) / ((data->Rmyo + Rgap / loc_l) * 4 * loc_l * loc_l / loc_a) + (vup - 2*Vm + vdown) / ((data->Rmyo + RgapupDown / (2 * loc_a)) * 8 * loc_a * 1.77245)) - (naiont+kiont+caiont+clont);		
		Vm = check_State_Change(Vm, dvdt, dt);

		if(Ca_i < CaT_Min) CaT_Min = Ca_i;
		if(Ca_i > CaT_Max) CaT_Max = Ca_i;
		if (Vm < Vmin) Vmin = Vm;
		
		if(dvdt > valmaxdvdt)
		{
			valmaxdvdt = dvdt;
			tmaxdvdt = t;
		}
		else
		{
			if(t < tmaxdvdt + 4 && Vm > vmax) vmax = Vm;
		}
		if(vmax > -50 && t > tmaxdvdt + 4 && Vm < (vrest + (1 - settings->APDRepLevel) * (vmax - vrest)) && APD == 0) 
		{
			APD = t - tmaxdvdt;
			APA = (vmax - vrest);
		}

}

void RAMM::Ion_HH_Na(int seg_i, int mem_i, double dt)
{
	
	if(fabs(settings->INaB - 1.0) < 1E-10)
	{
		segments[seg_i].membrane[mem_i].INa_junc = 0;
		segments[seg_i].membrane[mem_i].INa_sl = 0;
		return;
	}
	
	//double ENa = log(settings->Na_o/segments[seg_i].Na_i)/data->frt;
	double ENa_junc = log(settings->Na_o/segments[seg_i].membrane[mem_i].Na_junc)/data->frt;
	double ENa_sl = log(settings->Na_o/segments[seg_i].membrane[mem_i].Na_sl)/data->frt;

	double Po_Na1 = segments[seg_i].membrane[mem_i].m_Na * segments[seg_i].membrane[mem_i].m_Na * segments[seg_i].membrane[mem_i].m_Na;
	double Po_Na2 = 0.635*segments[seg_i].membrane[mem_i].h1_Na + 0.365*segments[seg_i].membrane[mem_i].h2_Na;
	double INa_junc_NP = data->Fjunc * data->GNa * Po_Na1 * Po_Na2 * settings->Na_o * Vm * data->F * data-> frt * (exp((Vm-ENa_junc)*data->frt)-1)/(exp(Vm * data->frt) - 1);
	double INa_sl_NP = (1-data->Fjunc) * data->GNa * Po_Na1 * Po_Na2 * settings->Na_o * Vm * data->F * data-> frt * (exp((Vm-ENa_sl)*data->frt)-1)/(exp(Vm * data->frt) - 1);

        segments[seg_i].membrane[mem_i].INa_junc = settings->INa_junc_scl * (1 - settings->INaB) * INa_junc_NP;
        segments[seg_i].membrane[mem_i].INa_sl = settings->INa_sl_scl * (1 - settings->INaB) * INa_sl_NP;

	if(seg_i == 0 && mem_i == 0) //only depends on Vm ?
	{
		// Sodium m gate
		am = -0.460 * (Vm + 44.4)/(exp(-(Vm + 44.4)/12.673) - 1);
		bm = 18.400 * exp(-(Vm + 44.4)/12.673);

		ah = 0.0449 * exp(-(Vm + 66.9)/5.57);
        	bh = 1.491/(1 + 323.3 * exp(-(Vm + 94.6)/12.9));
        	h_infinity = ah / (ah + bh);

		// Sodium h1 gate
		tauh1 = 30/(1 + exp((Vm + 40)/6)) + 0.15;

		// Sodium  h2 gate
		tauh2 = 120/(1 + exp((Vm + 60)/2)) + 0.45;
	}

	segments[seg_i].membrane[mem_i].m_Na = check_State_Change(segments[seg_i].membrane[mem_i].m_Na, am * (1 - segments[seg_i].membrane[mem_i].m_Na) - bm * segments[seg_i].membrane[mem_i].m_Na, dt);
		if(segments[seg_i].membrane[mem_i].m_Na < 0.0) segments[seg_i].membrane[mem_i].m_Na = 0.0;
		if(segments[seg_i].membrane[mem_i].m_Na > 1.0) segments[seg_i].membrane[mem_i].m_Na = 1.0;
		
	segments[seg_i].membrane[mem_i].h1_Na = check_State_Change(segments[seg_i].membrane[mem_i].h1_Na, (h_infinity - segments[seg_i].membrane[mem_i].h1_Na) / tauh1, dt);
        segments[seg_i].membrane[mem_i].h2_Na = check_State_Change(segments[seg_i].membrane[mem_i].h2_Na, (h_infinity - segments[seg_i].membrane[mem_i].h2_Na) / tauh2, dt);
	
}

void RAMM::Ion_HH_NaCa(int seg_i, int mem_i, double dt)
{
        
	if(fabs(settings->INaCaB - 1.0) < 1E-10)
        {
                segments[seg_i].membrane[mem_i].INaCa_junc = 0;
                segments[seg_i].membrane[mem_i].INaCa_sl = 0;
                return;
        }

        if(settings->clamp_NCX > -9.9)
        {
                segments[seg_i].membrane[mem_i].INaCa_junc = settings->clamp_NCX;
                segments[seg_i].membrane[mem_i].INaCa_sl = settings->clamp_NCX;
        }
        else
        {
                double Ca_junc_observed;
                if(mem_i == 0) Ca_junc_observed = segments[seg_i].caUnits[0].Ca_srs;
                if(mem_i == 1) Ca_junc_observed = segments[seg_i].caUnits[settings->N_CaDomains-1].Ca_srs;

		double Vm_s1 = exp(data->nu*Vm*data->frt) * settings->Ca_o;
                double Vm_s2 = exp((data->nu-1)*Vm*data->frt);

		if (settings->NCX_model==1)
		{
			double s1_junc = pow(segments[seg_i].membrane[mem_i].Na_junc,data->n_NaCa)*settings->Ca_o;
			double s1_sl = pow(segments[seg_i].membrane[mem_i].Na_sl,data->n_NaCa)*settings->Ca_o;
			double s2_junc = pow(settings->Na_o,data->n_NaCa) * Ca_junc_observed;
			double s2_sl = pow(settings->Na_o,data->n_NaCa) * segments[seg_i].membrane[mem_i].Ca_sl;
			double s3_junc = 1+data->dNaCa*(pow(settings->Na_o,data->n_NaCa))*Ca_junc_observed + pow(segments[seg_i].membrane[mem_i].Na_junc,data->n_NaCa)*settings->Ca_o;
			double s3_sl = 1+data->dNaCa*(pow(settings->Na_o,data->n_NaCa))*segments[seg_i].membrane[mem_i].Ca_sl + pow(segments[seg_i].membrane[mem_i].Na_sl,data->n_NaCa)* settings->Ca_o;

			segments[seg_i].membrane[mem_i].INaCa_junc = settings->INaCa_junc_scl * (1 - settings->INaCaB) * data->Fjunc * data->IbarNCX * (s1_junc*Vm_s1-s2_junc*Vm_s2)/s3_junc;
			segments[seg_i].membrane[mem_i].INaCa_sl = settings->INaCa_sl_scl * (1 - settings->INaCaB) * (1 - data->Fjunc) * data->IbarNCX * (s1_sl*Vm_s1-s2_sl*Vm_s2)/s3_sl;
		}
		else
		{                
               		double Na_o_p3 = settings->Na_o * settings->Na_o * settings->Na_o;
                	double Na_sl_p3 =  segments[seg_i].membrane[mem_i].Na_sl * segments[seg_i].membrane[mem_i].Na_sl * segments[seg_i].membrane[mem_i].Na_sl;
                	double Na_junc_p3 =  segments[seg_i].membrane[mem_i].Na_junc * segments[seg_i].membrane[mem_i].Na_junc * segments[seg_i].membrane[mem_i].Na_junc;
                
                	double Ka_junc = 1 / (1+pow(data->Kdact/Ca_junc_observed, 2));
                	double Ka_sl = 1 / (1+pow(data->Kdact/segments[seg_i].membrane[mem_i].Ca_sl, 2));

                	if(settings->Ca_o > 0)
                	{
                        	double Vm_s1 = exp(data->nu*Vm*data->frt) * settings->Ca_o;
                        	double s1_junc = Vm_s1 * Na_junc_p3;
                        	double s1_sl = Vm_s1 * Na_sl_p3;
                        	double Vm_s2 = exp((data->nu-1)*Vm*data->frt);
                        	double s2_junc = Vm_s2 * Na_o_p3 * Ca_junc_observed;
                        	double s3_junc = data->KmCai*Na_o_p3*(1+pow(segments[seg_i].membrane[mem_i].Na_junc/data->KmNai,3)) + data->KmNao_p3*Ca_junc_observed*(1+Ca_junc_observed/data->KmCai) + data->KmCao*Na_junc_p3 + Na_junc_p3*settings->Ca_o + Na_o_p3*Ca_junc_observed;
                        	double s2_sl = Vm_s2 * Na_o_p3 * segments[seg_i].membrane[mem_i].Ca_sl;
                        	double s3_sl = data->KmCai*Na_o_p3*(1+pow(segments[seg_i].membrane[mem_i].Na_sl/data->KmNai,3)) + data->KmNao_p3*segments[seg_i].membrane[mem_i].Ca_sl*(1+segments[seg_i].membrane[mem_i].Ca_sl/data->KmCai) + data->KmCao*Na_sl_p3 + Na_sl_p3*settings->Ca_o + Na_o_p3*segments[seg_i].membrane[mem_i].Ca_sl;

 	                        segments[seg_i].membrane[mem_i].INaCa_junc = settings->INaCa_junc_scl * (1 - settings->INaCaB) * data->Fjunc * data->IbarNCX * Ka_junc * (s1_junc-s2_junc)/s3_junc/(1+data->ksat*Vm_s2);
                        	segments[seg_i].membrane[mem_i].INaCa_sl = settings->INaCa_sl_scl * (1 - settings->INaCaB) * (1 - data->Fjunc) * data->IbarNCX * Ka_sl * (s1_sl-s2_sl)/s3_sl/(1+data->ksat*Vm_s2);
                	}
                	else
                	{
                        	segments[seg_i].membrane[mem_i].INaCa_junc = 0;
                        	segments[seg_i].membrane[mem_i].INaCa_sl = 0;
                	}
		}
        }
        
}

void RAMM::Ion_HH_NaK(int seg_i, int mem_i, double dt)
{ 
	
	if(fabs(settings->INaKB - 1.0) < 1E-10)
        {
                segments[seg_i].membrane[mem_i].INaK_junc = 0;
		segments[seg_i].membrane[mem_i].INaK_sl = 0;
                return;
        }

	if(seg_i == 0 && mem_i == 0)
	{
		// Only need to calculate these once every time step
		double INaK_Po_junc = (pow(segments[seg_i].membrane[mem_i].Na_junc,data->kNaK_Na_n)/(pow(segments[seg_i].membrane[mem_i].Na_junc,data->kNaK_Na_n)+pow(data->kNaK_Na,data->kNaK_Na_n)));
		double INaK_Po_sl = (pow(segments[seg_i].membrane[mem_i].Na_sl,data->kNaK_Na_n)/(pow(segments[seg_i].membrane[mem_i].Na_sl,data->kNaK_Na_n)+pow(data->kNaK_Na,data->kNaK_Na_n)));

		// Aslanidi
		inak_junc = data->iNaKmax * INaK_Po_junc * (settings->K_o/(settings->K_o+data->kNaK_K)) * 1.6/(1.5+exp(-(Vm+60)/40));
		inak_sl = data->iNaKmax * INaK_Po_sl * (settings->K_o/(settings->K_o+data->kNaK_K)) * 1.6/(1.5+exp(-(Vm+60)/40));

	}
	
	segments[seg_i].membrane[mem_i].INaK_junc = settings->INaK_junc_scl * (1 - settings->INaKB) * data->Fjunc * inak_junc;
	segments[seg_i].membrane[mem_i].INaK_sl = settings->INaK_sl_scl * (1 - settings->INaKB) * (1.0-data->Fjunc) *  inak_sl;

}

void RAMM::Ion_HH_Nab(int seg_i, int mem_i, double dt)
{

	if(fabs(settings->INabB - 1.0) < 1E-10)
        {
                segments[seg_i].membrane[mem_i].INab_junc = 0;
		segments[seg_i].membrane[mem_i].INab_sl = 0;
                return;
        }
	double ENa_junc = log(settings->Na_o/segments[seg_i].membrane[mem_i].Na_junc)/data->frt;
	double ENa_sl = log(settings->Na_o/segments[seg_i].membrane[mem_i].Na_sl)/data->frt;
	segments[seg_i].membrane[mem_i].INab_junc = settings->INab_junc_scl * (1 - settings->INabB) * data->Fjunc * data->GNab * (Vm-ENa_junc);
	segments[seg_i].membrane[mem_i].INab_sl = settings->INab_sl_scl * (1 - settings->INabB) * (1-data->Fjunc) * data->GNab * (Vm-ENa_sl);
}

void RAMM::Ion_HH_Cab(int seg_i, int mem_i, double dt)
{

	if(fabs(settings->ICabB - 1.0) < 1E-10)
        {
                segments[seg_i].membrane[mem_i].ICab_junc = 0;
		segments[seg_i].membrane[mem_i].ICab_sl = 0;
                return;
        }
        double ECa_junc = log(settings->Ca_o/segments[seg_i].membrane[mem_i].Ca_junc)/data->frt;
        double ECa_sl = log(settings->Ca_o/segments[seg_i].membrane[mem_i].Ca_sl)/data->frt;
        segments[seg_i].membrane[mem_i].ICab_junc = settings->ICab_junc_scl * (1 - settings->ICabB) * data->Fjunc * data->GCab * (Vm-ECa_junc);
        segments[seg_i].membrane[mem_i].ICab_sl = settings->ICab_sl_scl * (1 - settings->ICabB) * (1-data->Fjunc) * data->GCab * (Vm-ECa_sl);
}

void RAMM::Ion_HH_Clb(int seg_i, int mem_i, double dt)
{

        if(fabs(settings->IClbB - 1.0) < 1E-10)
        {
                segments[seg_i].membrane[mem_i].IClb = 0;
                return;
        }

//        if(seg_i == 0 && mem_i == 0)
//        {
		IClb = settings->IClb_scl * (1 - settings->IClbB) * data->GClb * (Vm - EB_Cl) * (1 + exp((Vm - E_Cl - 36.95)/74.514));
//	}
        segments[seg_i].membrane[mem_i].IClb = IClb;

}

void RAMM::Ion_HH_ClCa(int seg_i, int mem_i, double dt)
{

	if(fabs(settings->IClCaB - 1.0) < 1E-10)
        {
                segments[seg_i].membrane[mem_i].IClCa_junc = 0;
		segments[seg_i].membrane[mem_i].IClCa_sl = 0;
                return;
        }
        segments[seg_i].membrane[mem_i].IClCa_junc = settings->IClCa_scl * (1 - settings->IClCaB) * data->Fjunc * data->GClCa / (1+data->KdClCa/segments[seg_i].membrane[mem_i].Ca_junc)*(Vm-E_Cl);
        segments[seg_i].membrane[mem_i].IClCa_sl = settings->IClCa_scl* (1 - settings->IClCaB) * (1 - data->Fjunc) * data->GClCa / (1+data->KdClCa/segments[seg_i].membrane[mem_i].Ca_sl)*(Vm-E_Cl);
}


void RAMM::Ion_HH_K1(int seg_i, int mem_i, double dt)
{
	if(fabs(settings->IK1B - 1.0) < 1E-10)
	{
		segments[seg_i].membrane[mem_i].IK1 = 0;
		return;
	}
	
	double EK = log(settings->K_o/segments[seg_i].K_i)/data->frt;
	double k3 = pow(settings->K_o/(settings->K_o+data->KmK1), 3);
	double ke = exp(data->steepK1*(Vm-EK-data->shiftK1)*data->frt);

	segments[seg_i].membrane[mem_i].IK1 = settings->IK1_scl * (1 - settings->IK1B) * data->GK1max * (Vm-EK) * k3/(1 + ke);
}

void RAMM::Ion_HH_sus(int seg_i, int mem_i, double dt)
{
	
	if(fabs(settings->IsusB - 1.0) < 1E-10)
	{
		segments[seg_i].membrane[mem_i].Isus = 0;
		return;
	}

	//if(seg_i == 0 && mem_i == 0)
        //{
		Isus = settings->Isus_scl * (1 - settings->IsusB) * data->Gsus * (Vm + 70);
	//}
	
	segments[seg_i].membrane[mem_i].Isus = Isus;
	
}


void RAMM::Ion_HH_To(int seg_i, int mem_i, double dt)
{

	if(fabs(settings->ITo1B - 1.0) < 1E-10)
	{
		segments[seg_i].membrane[mem_i].ITo = 0;
		return;
	}
	
	double EK = log(settings->K_o/segments[seg_i].K_i)/data->frt;
	double S1 = 0.59*segments[seg_i].membrane[mem_i].to_s1*segments[seg_i].membrane[mem_i].to_s1*segments[seg_i].membrane[mem_i].to_s1;
	double S2 = 0.41*segments[seg_i].membrane[mem_i].to_s2*segments[seg_i].membrane[mem_i].to_s2*segments[seg_i].membrane[mem_i].to_s2;
	double S3 = 0.6*segments[seg_i].membrane[mem_i].to_s3*segments[seg_i].membrane[mem_i].to_s3*segments[seg_i].membrane[mem_i].to_s3*segments[seg_i].membrane[mem_i].to_s3*segments[seg_i].membrane[mem_i].to_s3*segments[seg_i].membrane[mem_i].to_s3;

	segments[seg_i].membrane[mem_i].ITo = settings->ITo_scl * (1 - settings->ITo1B) * data->GToFast * segments[seg_i].membrane[mem_i].to_r * (S1+S2) * (S3+0.4) * (Vm-EK);

	if(seg_i == 0 && mem_i == 0)
        {
	 	// Ito r gate
		alpha_r = 386.6 * exp(Vm/12);
		beta_r = 8.011 * exp(Vm/-7.2);
		rss = 1/(1 + exp((Vm+15)/-5.633));
		taur = 1000 * (1/(alpha_r+beta_r) + 4E-4);

		// Ito s1 gate
		s1ss = 1/(1 + exp((Vm + 28.29)/7.06));
		taus1 = 1000 * (0.5466/(1 + exp((Vm + 32.8)/0.1)) + 0.0204);

		// Ito s2 gate
		s2ss = 1/(1 + exp((Vm + 28.29)/7.06));
		//taus2 = 1000 * (5.75/(1+exp((Vm + 32.8)/0.1)) + 0.45/(1+exp((Vm-13.54)/-13.97)));
		if(settings->Ito_model == 2){ // Aslanidi formulation
			taus2 = 1000 * (0.189/(1+exp((Vm + 32.8)/0.1)) + 0.45*exp((Vm - 13.54)/-13.97));}
		else{ //Lindblad formulation
			taus2 = 1000 * (0.57/(1+exp((Vm + 32.8)/0.1)) + 0.45/pow(exp((Vm - 13.54)/-13.97),2) + 0.020);
		}
		// Ito s3 gate
		s3ss = ((1/(1 + exp((Vm + 50.67)/27.38))) + 0.666)/1.666;
		taus3 = 1000 * (7.5/(1+exp((Vm + 23)/0.5)) + 0.5);
	}

	segments[seg_i].membrane[mem_i].to_r = check_State_Change(segments[seg_i].membrane[mem_i].to_r, (rss - segments[seg_i].membrane[mem_i].to_r) / taur, dt);
	segments[seg_i].membrane[mem_i].to_s1 = check_State_Change(segments[seg_i].membrane[mem_i].to_s1, (s1ss - segments[seg_i].membrane[mem_i].to_s1) / taus1, dt);
	segments[seg_i].membrane[mem_i].to_s2 = check_State_Change(segments[seg_i].membrane[mem_i].to_s2, (s2ss - segments[seg_i].membrane[mem_i].to_s2) / taus2, dt);
	segments[seg_i].membrane[mem_i].to_s3 = check_State_Change(segments[seg_i].membrane[mem_i].to_s3, (s3ss - segments[seg_i].membrane[mem_i].to_s3) / taus3, dt); 

}

void RAMM::Ion_HH_Ks(int seg_i, int mem_i, double dt)
{
	
	if(fabs(settings->IKsB - 1.0) < 1E-10)
	{
		segments[seg_i].membrane[mem_i].IKs = 0;
		return;
	}
	
	double EK = log(settings->K_o/K_i)/data->frt;
	//double EK = (1/data->frt) * log((settings->K_o + data->pNaK * settings->Na_o)/(K_i + data->pNaK * segments[seg_i].membrane[mem_i].Na_sl));
	segments[seg_i].membrane[mem_i].IKs = settings->IKs_scl * (1 - settings->IKsB) * data->GKs * segments[seg_i].membrane[mem_i].z_Ks * (Vm-EK);
	
	if(seg_i == 0 && mem_i == 0) 
        {
		//z gate
		alpha_z = 1.66 * exp(Vm/69.452);
	        beta_z = 0.3 * exp(Vm/-21.826);
		zss_Ks = 1/(1 + exp((Vm-0.9)/-13.8));
		tauz_Ks = 1000*(1/(alpha_z + beta_z) + 0.06);
	}
	segments[seg_i].membrane[mem_i].z_Ks = check_State_Change(segments[seg_i].membrane[mem_i].z_Ks, (zss_Ks - segments[seg_i].membrane[mem_i].z_Ks) / tauz_Ks, dt);
		
}

void RAMM::Ion_HH_Kr(int seg_i, int mem_i, double dt)
{

	if(fabs(settings->IKrB - 1.0) < 1E-10)
        {
                segments[seg_i].membrane[mem_i].IKr = 0;
                return;
        }

        double EK = log(settings->K_o/segments[seg_i].K_i)/data->frt;
        segments[seg_i].membrane[mem_i].IKr = settings->IKr_scl * (1 - settings->IKrB) * data->GKr * segments[seg_i].membrane[mem_i].pa_Kr * segments[seg_i].membrane[mem_i].pi_Kr * (Vm-EK);

        //if(seg_i == 0 && mem_i == 0)
        //{
		//pa gate
                alpha_pa = 9 * exp(Vm/25.371);
                beta_pa = 1.3 * exp(Vm/-13.026);
                pass_Kr = 1/(1 + exp((Vm+5.1)/-7.4));
                taupa_Kr = 1000/(alpha_pa + beta_pa);

		//pi gate
		alpha_pi = 100 * exp(Vm/-54.645);
		beta_pi = 656 * exp(Vm/106.157);
		//piss_Kr = 1/(1 + exp((Vm+47.3921)/18.6603));
		piss_Kr = alpha_pi/(alpha_pi + beta_pi);
		taupi_Kr = 1000/(alpha_pi + beta_pi);


		//cout << "seg_i = " << seg_i << "   mem_i = " << mem_i << "   Vm = " << Vm << "   pass_Kr = " << pass_Kr << "  piss_kr = " << piss_Kr << "  taupa_kr = " << taupa_Kr << "   taupi_Kr = " << taupi_Kr << "   pa_Kr = " << segments[seg_i].membrane[mem_i].pa_Kr << "   pi_Kr = " << segments[seg_i].membrane[mem_i].pi_Kr << endl;
		//GenPurpose::wait(1);

        //}

        segments[seg_i].membrane[mem_i].pa_Kr = check_State_Change(segments[seg_i].membrane[mem_i].pa_Kr, (pass_Kr - segments[seg_i].membrane[mem_i].pa_Kr) / taupa_Kr, dt);
	segments[seg_i].membrane[mem_i].pi_Kr = check_State_Change(segments[seg_i].membrane[mem_i].pi_Kr, (piss_Kr - segments[seg_i].membrane[mem_i].pi_Kr) / taupi_Kr, dt);

}

void RAMM::Ion_HH_SK(int seg_i, int mem_i, double dt)
{
	if(fabs(settings->ISKB - 1.0) < 1E-10)
        {
                segments[seg_i].membrane[mem_i].ISK_junc = 0;
                segments[seg_i].membrane[mem_i].ISK_sl = 0;
		return;
        }

        double EK = log(settings->K_o/segments[seg_i].K_i)/data->frt;

	ISK_junc = (1 - settings->ISKB) * settings->ISK_scl * data->GKCa * data->Fjunc * segments[seg_i].membrane[mem_i].O_KCa_junc*(Vm - EK)/(1 + exp((Vm - EK + 120)/45));
	ISK_sl = (1 - settings->ISKB) * settings->ISK_scl * data->GKCa * (1.0-data->Fjunc) * segments[seg_i].membrane[mem_i].O_KCa_sl*(Vm - EK)/(1 + exp((Vm - EK + 120)/45));

	segments[seg_i].membrane[mem_i].ISK_junc = ISK_junc;
	segments[seg_i].membrane[mem_i].ISK_sl = ISK_sl;

        segments[seg_i].membrane[mem_i].O_KCa_junc = check_State_Change(segments[seg_i].membrane[mem_i].O_KCa_junc, (1-segments[seg_i].membrane[mem_i].O_KCa_junc)*data->KCa_on * segments[seg_i].membrane[mem_i].Ca_junc*segments[seg_i].membrane[mem_i].Ca_junc - segments[seg_i].membrane[mem_i].O_KCa_junc * data->KCa_off, dt);
	segments[seg_i].membrane[mem_i].O_KCa_sl = check_State_Change(segments[seg_i].membrane[mem_i].O_KCa_sl, (1-segments[seg_i].membrane[mem_i].O_KCa_sl)*data->KCa_on * segments[seg_i].membrane[mem_i].Ca_sl*segments[seg_i].membrane[mem_i].Ca_sl - segments[seg_i].membrane[mem_i].O_KCa_sl * data->KCa_off, dt);

}

void RAMM::Ion_HH_KAch(int seg_i, int mem_i, double dt)
{

	double PM2_KAch;

	if(fabs(settings->IKAchB - 1.0) < 1E-10)
        {
                segments[seg_i].membrane[mem_i].IKAch = 0;
                return;
        }

	 double EK = log(settings->K_o/segments[seg_i].K_i)/data->frt;
	
        //if(seg_i == 0 && mem_i == 0)
        //{
		if((settings->fvagal<100) || (settings->fvagal>25)){  
			PM2_KAch = 1.026/(1 + exp((settings->fvagal+11.05)/-7.5095)) - 0.99;
		}
		else{
			PM2_KAch = 0.0006;
		}
		alpha_a = 17 * exp(0.0133*(Vm + 40))/1000;  //[per millisecond]
		beta_a = 12.32/(1 + 0.0042/settings->Ach)/1000; //[per millisecond]
	//}

	IKAch = (1 - settings->IKAchB) * settings->IKAch_scl * data->GKAch * PM2_KAch * (Vm - EK) * segments[seg_i].membrane[mem_i].a_KAch;

	segments[seg_i].membrane[mem_i].IKAch = IKAch;

        segments[seg_i].membrane[mem_i].a_KAch = check_State_Change(segments[seg_i].membrane[mem_i].a_KAch, (beta_a*(1 - segments[seg_i].membrane[mem_i].a_KAch) - alpha_a*segments[seg_i].membrane[mem_i].a_KAch), dt);
}

void RAMM::Ion_HH_CaL(int seg_i, int mem_i, double dt)
{
 
	if(fabs(settings->ICaLB - 1.0) < 1E-10)
        {
                segments[seg_i].membrane[mem_i].ICaL_junc = 0;
                segments[seg_i].membrane[mem_i].ICaL_sl = 0;
                return;
        }

	//double fact_AF = 1.0;
	if(settings->clamp_ICaL_CDI < 0) settings->clamp_ICaL_CDI = 1; //This means no clamping
	
	//double Ca_junc_observed;
	//if(mem_i == 0) Ca_junc_observed = segments[seg_i].caUnits[0].Ca_srs;
	//if(mem_i == 1) Ca_junc_observed = segments[seg_i].caUnits[settings->N_CaDomains-1].Ca_srs;
		
	// NP Currents	
	if(settings->Ca_o > 0)
	{
		double d_prime = 1/(1 + exp((Vm-23)/-12));
		double ical = (1 - settings->ICaLB) * data->GCaL * (segments[seg_i].membrane[mem_i].d_L*segments[seg_i].membrane[mem_i].f_L+d_prime) * (Vm-50);

		if(seg_i == 0 && mem_i == 0)
        	{
			//d_L gate
			E0_alpha_d_L = Vm + 45;
			E0_beta_d_L = Vm + 5;
			E10 = Vm + 10;
			alpha_d_L = (-16.72 * E0_alpha_d_L/(exp(E0_alpha_d_L/-2.5)-1)-50 * E10/(exp(E10/-4.808)-1))/1000;
			beta_d_L = (4.48 * E0_beta_d_L/(exp(E0_beta_d_L/2.5)-1))/1000;
			d_L_infinity = 1/(1+exp((E10 + 0.95)/-6.6));
			tau_d_L = 1/(alpha_d_L + beta_d_L);

			//f_L gate
			E0_f_L = Vm + 18;
			alpha_f_L = (8.49 * E0_f_L/(exp(E0_f_L/4)-1))/1000;
			beta_f_L = (67.922/(1+exp(E0_f_L/-4)))/1000;
			f_L_infinity = alpha_f_L/(alpha_f_L + beta_f_L);
			tau_f_L = 1/(alpha_f_L + beta_f_L);
		}

		segments[seg_i].membrane[mem_i].d_L = check_State_Change(segments[seg_i].membrane[mem_i].d_L, (d_L_infinity - segments[seg_i].membrane[mem_i].d_L) / tau_d_L, dt);
		segments[seg_i].membrane[mem_i].f_L = check_State_Change(segments[seg_i].membrane[mem_i].f_L, (f_L_infinity - segments[seg_i].membrane[mem_i].f_L) / tau_f_L, dt);


		//Total currents
		segments[seg_i].membrane[mem_i].ICaL_junc = settings->ICaL_junc_scl * data->Fjunc * ical;
		segments[seg_i].membrane[mem_i].ICaL_sl = settings->ICaL_sl_scl * (1.0-data->Fjunc) * ical;
		segments[seg_i].membrane[mem_i].ICaL = segments[seg_i].membrane[mem_i].ICaL_junc + segments[seg_i].membrane[mem_i].ICaL_sl;
		//segments[seg_i].membrane[mem_i].ICaL_Na_junc = ICaL_Na_junc_NP;
		//segments[seg_i].membrane[mem_i].ICaL_Na_sl = ICaL_Na_sl_NP;
		//segments[seg_i].membrane[mem_i].ICaL_Na = segments[seg_i].membrane[mem_i].ICaL_Na_junc + segments[seg_i].membrane[mem_i].ICaL_Na_sl;
		//segments[seg_i].membrane[mem_i].ICaL_K = ICaL_K_NP;
		segments[seg_i].membrane[mem_i].ICaL_tot = segments[seg_i].membrane[mem_i].ICaL;
	}
	else
	{
		segments[seg_i].membrane[mem_i].ICaL_junc = 0;
		segments[seg_i].membrane[mem_i].ICaL_sl = 0;
		segments[seg_i].membrane[mem_i].ICaL = 0;
		//segments[seg_i].membrane[mem_i].ICaL_Na_junc = 0;
		//segments[seg_i].membrane[mem_i].ICaL_Na_sl = 0;
		//segments[seg_i].membrane[mem_i].ICaL_Na = 0;
		//segments[seg_i].membrane[mem_i].ICaL_K = 0;
		segments[seg_i].membrane[mem_i].ICaL_tot = 0;
	}
	
}


void RAMM::Ion_HH_CaT(int seg_i, int mem_i, double dt)
{

        if(fabs(settings->ICaTB - 1.0) < 1E-10)
        {
                segments[seg_i].membrane[mem_i].ICaT_junc = 0;
                segments[seg_i].membrane[mem_i].ICaT_sl = 0;
                return;
        }

        // NP Currents
        if(settings->Ca_o > 0)
        {
		double icat = (1 - settings->ICaTB) * data->GCaT * segments[seg_i].membrane[mem_i].d_T * segments[seg_i].membrane[mem_i].f_T * (Vm-38);

		//Total currents
                segments[seg_i].membrane[mem_i].ICaT_junc = settings->ICaT_junc_scl * data->Fjunc * icat;
                segments[seg_i].membrane[mem_i].ICaT_sl = settings->ICaT_sl_scl * (1.0-data->Fjunc) * icat;
                segments[seg_i].membrane[mem_i].ICaT = segments[seg_i].membrane[mem_i].ICaT_junc + segments[seg_i].membrane[mem_i].ICaT_sl;
                segments[seg_i].membrane[mem_i].ICaT_tot = segments[seg_i].membrane[mem_i].ICaT;

		if(seg_i == 0 && mem_i == 0)
                {
                	//d_T gate
                	E0_d_T = Vm + 23.3;
                	alpha_d_T = (674.173 * exp(E0_d_T/30))/1000;
                	beta_d_T = (674.173 * exp(E0_d_T/-30))/1000;
                	d_T_infinity = 1/(1 + exp((E0_d_T-0.3)/-6.1));
                	tau_d_T = 1/(alpha_d_T + beta_d_T);

                	//f_T gate
                	E0_f_T = Vm + 75;
                	alpha_f_T = (9.637*exp(E0_f_T/-83.3))/1000;
                	beta_f_T = (9.637*exp(E0_f_T/15.38))/1000;
                	f_T_infinity = f_T_infinity = alpha_f_T/(alpha_f_T+beta_f_T);
                	tau_f_T = 1/(alpha_f_T+beta_f_T);
		}

                segments[seg_i].membrane[mem_i].d_T = check_State_Change(segments[seg_i].membrane[mem_i].d_T, (d_T_infinity - segments[seg_i].membrane[mem_i].d_T) / tau_d_T, dt);
                segments[seg_i].membrane[mem_i].f_T = check_State_Change(segments[seg_i].membrane[mem_i].f_T, (f_T_infinity - segments[seg_i].membrane[mem_i].f_T) / tau_f_T, dt);
        }
        else
        {
                segments[seg_i].membrane[mem_i].ICaT_junc = 0;
                segments[seg_i].membrane[mem_i].ICaT_sl = 0;
                segments[seg_i].membrane[mem_i].ICaT = 0;
                segments[seg_i].membrane[mem_i].ICaT_tot = 0;
        }

}

void RAMM::Ion_HH_CaP(int seg_i, int mem_i, double dt)
{
	//double Ca_junc_observed;
        //if(mem_i == 0) Ca_junc_observed = segments[seg_i].caUnits[0].Ca_srs;
        //if(mem_i == 1) Ca_junc_observed = segments[seg_i].caUnits[settings->N_CaDomains-1].Ca_srs;
	segments[seg_i].membrane[mem_i].ICaP_junc = (1 - settings->ICaPB)*settings->ICaP_junc_scl*data->Fjunc*data->ICaP_max*segments[seg_i].membrane[mem_i].Ca_junc/(segments[seg_i].membrane[mem_i].Ca_junc+data->k_CaP);
	segments[seg_i].membrane[mem_i].ICaP_sl = (1 - settings->ICaPB)*settings->ICaP_sl_scl * (1.0-data->Fjunc)*data->ICaP_max*segments[seg_i].membrane[mem_i].Ca_sl/(segments[seg_i].membrane[mem_i].Ca_sl+data->k_CaP);
}

// Original Grandi et al. Circ Res 2011 ICaL formulation - no longer used.
void RAMM::Ion_Ca_Handling(int seg_i, int mem_i, double dt)
{
	
	double Ca_junc_observed;
	if(mem_i == 0) Ca_junc_observed = segments[seg_i].caUnits[0].Ca_srs;
	if(mem_i == 1) Ca_junc_observed = segments[seg_i].caUnits[settings->N_CaDomains-1].Ca_srs;
	
	double temp_IpCa, ECa_junc, ECa_sl;
	
	// Junctional SRS Component
	double ca_junc_pow = pow(Ca_junc_observed,1.6);
	temp_IpCa = data->IbarSLCaP * ca_junc_pow / (data->KmPCaPow1_6+ca_junc_pow);
	ECa_junc = 300;
	if(Ca_junc_observed > 1E-8)
	{
		ECa_junc = (1/data->frt/2) * log(settings->Ca_o/Ca_junc_observed);
	}
	segments[seg_i].membrane[mem_i].IpCa_junc = (1 - settings->IpCaB) * data->Fjunc * temp_IpCa;
	
	// Subsarcolemmal component
	double ca_sl_pow = pow(segments[seg_i].membrane[mem_i].Ca_sl,1.6);
	temp_IpCa = data->IbarSLCaP * ca_sl_pow / (data->KmPCaPow1_6+ca_sl_pow);
	ECa_sl = 300;
	if(segments[seg_i].membrane[mem_i].Ca_sl > 1E-8)
	{
		ECa_sl = (1/data->frt/2) * log(settings->Ca_o/segments[seg_i].membrane[mem_i].Ca_sl);
	}
	segments[seg_i].membrane[mem_i].IpCa_sl = (1 - settings->IpCaB) * (1 - data->Fjunc) * temp_IpCa;
	
	if(settings->Ca_o > 0)
	{
	        //double ECa_junc = log(settings->Ca_o/segments[seg_i].membrane[mem_i].Ca_junc)/data->frt;
        	//double ECa_sl = log(settings->Ca_o/segments[seg_i].membrane[mem_i].Ca_sl)/data->frt;
        	segments[seg_i].membrane[mem_i].ICab_junc = settings->ICab_junc_scl * (1 - settings->ICabB) * data->Fjunc * data->GCab * (Vm-ECa_junc);
        	segments[seg_i].membrane[mem_i].ICab_sl = settings->ICab_sl_scl * (1 - settings->ICabB) * (1-data->Fjunc) * data->GCab * (Vm-ECa_sl);
	}
	else
	{
		segments[seg_i].membrane[mem_i].ICab_junc = 0;
		segments[seg_i].membrane[mem_i].ICab_sl = 0;
	}
	
	if(showlog)
	{
		cout << "  seg_i = " << seg_i << ", mem_i = " << mem_i << ", Ca_sl = " << segments[seg_i].membrane[mem_i].Ca_sl << ", temp_IpCa = " << temp_IpCa << ", ECa_sl = " << ECa_sl << ", ICab_sl = " << segments[seg_i].membrane[mem_i].ICab_sl << endl;
	}

	
}

void RAMM::updateCaDomain(double dt, double t, int segment, int ca_domain, double Ca_i_L_Domain, double Ca_i_R_Domain, double Ca_i_U_Segment, double Ca_i_D_Segment, double Ca_sr_L_Domain, double Ca_sr_R_Domain, double Ca_sr_U_Segment, double Ca_sr_D_Segment, double Ca_srs_L_Domain, double Ca_srs_R_Domain, double Ca_srs_U_Segment, double Ca_srs_D_Segment)
{
	
	// Determine local ionic currents and fluxes
	// ----------------------------------------------------
	if(settings->RyRModelType == 2)	Ca_CICR(segment, ca_domain, dt, t);
	if(settings->RyRModelType == 3)	Ca_CICR_Stoch(segment, ca_domain, dt, t);
	
	// Allow simulation at finer resolution
	int numsteps = max(1, (int)(1 + dt / 0.075)); //settings->PARAM_dTmin);		
	dt = dt / numsteps;
	for(int abc=0; abc<numsteps; abc++)
	{
		
	// Diffusive fluxes between neighboring CaRU's
	// ----------------------------------------------------
		double J_diff_i_subcell = (Ca_i_L_Domain - segments[segment].caUnits[ca_domain].Ca_i + Ca_i_R_Domain - segments[segment].caUnits[ca_domain].Ca_i) / data->tau_diff_Cai_Domain + (Ca_i_U_Segment - segments[segment].caUnits[ca_domain].Ca_i + Ca_i_D_Segment - segments[segment].caUnits[ca_domain].Ca_i) / data->tau_diff_Cai_Segment;
		double J_diff_sr_subcell = (Ca_sr_L_Domain - segments[segment].caUnits[ca_domain].Ca_sr + Ca_sr_R_Domain - segments[segment].caUnits[ca_domain].Ca_sr) / data->tau_diff_Casr_Domain + (Ca_sr_U_Segment - segments[segment].caUnits[ca_domain].Ca_sr + Ca_sr_D_Segment - segments[segment].caUnits[ca_domain].Ca_sr) / data->tau_diff_Casr_Segment;
		
		double J_diff_srs_subcell = (Ca_srs_L_Domain - segments[segment].caUnits[ca_domain].Ca_srs + Ca_srs_R_Domain - segments[segment].caUnits[ca_domain].Ca_srs) / data->tau_diff_Casrs_Domain; 
		J_diff_srs_subcell += (Ca_srs_U_Segment - segments[segment].caUnits[ca_domain].Ca_srs + Ca_srs_D_Segment - segments[segment].caUnits[ca_domain].Ca_srs) / data->tau_diff_Casrs_Segment;
		
		double J_diff_SRS_i = (segments[segment].caUnits[ca_domain].Ca_srs - segments[segment].caUnits[ca_domain].Ca_i) / data->tau_diff_SRS_i;
			
	//Update local ion concentrations
	// ----------------------------------------------------
		// Cytosolic Ca Buffers
		double dBuff_Ca_TnCL = data->kon_tncl * segments[segment].caUnits[ca_domain].Ca_i * (settings->Ca_Buff_Factor * data->Bmax_TnClow - segments[segment].caUnits[ca_domain].Buff_Ca_TnCL) - data->koff_tncl * segments[segment].caUnits[ca_domain].Buff_Ca_TnCL;            		// TnCL      [mM/ms]
		double dBuff_Ca_TnCHc = data->kon_tnchca * segments[segment].caUnits[ca_domain].Ca_i * (settings->Ca_Buff_Factor * data->Bmax_TnChigh - segments[segment].caUnits[ca_domain].Buff_Ca_TnCHc - segments[segment].caUnits[ca_domain].Buff_Ca_TnCHm) - data->koff_tnchca * segments[segment].caUnits[ca_domain].Buff_Ca_TnCHc; 	// TnCHc     [mM/ms]
		double dBuff_Ca_TnCHm = data->kon_tnchmg * settings->Mg_i * (settings->Ca_Buff_Factor * data->Bmax_TnChigh - segments[segment].caUnits[ca_domain].Buff_Ca_TnCHc - segments[segment].caUnits[ca_domain].Buff_Ca_TnCHm) - data->koff_tnchmg * segments[segment].caUnits[ca_domain].Buff_Ca_TnCHm;   // TnCHm     [mM/ms]
		double dBuff_Ca_CaM = data->kon_cam * segments[segment].caUnits[ca_domain].Ca_i * (settings->Ca_Buff_Factor * data->Bmax_CaM - segments[segment].caUnits[ca_domain].Buff_Ca_CaM) - data->koff_cam * segments[segment].caUnits[ca_domain].Buff_Ca_CaM;                 // CaM       [mM/ms]
		double dBuff_Ca_Myosin_ca = data->kon_myoca * segments[segment].caUnits[ca_domain].Ca_i * (settings->Ca_Buff_Factor * data->Bmax_myosin - segments[segment].caUnits[ca_domain].Buff_Ca_Myosin_ca - segments[segment].caUnits[ca_domain].Buff_Ca_Myosin_mg) - data->koff_myoca * segments[segment].caUnits[ca_domain].Buff_Ca_Myosin_ca;    // Myosin_ca [mM/ms]
		double dBuff_Ca_Myosin_mg = data->kon_myomg * settings->Mg_i * (settings->Ca_Buff_Factor * data->Bmax_myosin - segments[segment].caUnits[ca_domain].Buff_Ca_Myosin_ca - segments[segment].caUnits[ca_domain].Buff_Ca_Myosin_mg) - data->koff_myomg * segments[segment].caUnits[ca_domain].Buff_Ca_Myosin_mg;      // Myosin_mg [mM/ms]
		double dBuff_Ca_SRB = data->kon_sr * segments[segment].caUnits[ca_domain].Ca_i * (settings->Ca_Buff_Factor * data->Bmax_SR - segments[segment].caUnits[ca_domain].Buff_Ca_SRB) - data->koff_sr * segments[segment].caUnits[ca_domain].Buff_Ca_SRB;                    // SRB       [mM/ms]
		double dBuff_Ca_EGTA_cyt = (0.005 / dt) * (data->kon_EGTA * segments[segment].caUnits[ca_domain].Ca_i * (data->EGTA - segments[segment].caUnits[ca_domain].Buff_Ca_EGTA_cyt) - data->koff_EGTA * segments[segment].caUnits[ca_domain].Buff_Ca_EGTA_cyt);
		double J_CaB_cytosol = dBuff_Ca_TnCL + dBuff_Ca_TnCHc + dBuff_Ca_TnCHm + dBuff_Ca_CaM + dBuff_Ca_Myosin_ca + dBuff_Ca_Myosin_mg + dBuff_Ca_SRB + dBuff_Ca_EGTA_cyt;

		//cout << "Jup = " << segments[segment].caUnits[ca_domain].Jup << ", J_CaB_cytosol = " << J_CaB_cytosol << ", J_ca_slmyo = " << data->J_ca_slmyo << ", J_diff_i_subcell = " << J_diff_i_subcell << ", J_diff_SRS_i = " << J_diff_SRS_i << endl;
		//cout << "Ca_i = " << segments[segment].caUnits[ca_domain].Ca_i << ", Ca_sr = " << segments[segment].caUnits[ca_domain].Ca_sr << ", Ca_srs = " << segments[segment].caUnits[ca_domain].Ca_srs << endl;
		//cout << "dBuff_Ca_TnCL = " << dBuff_Ca_TnCL << ", dBuff_Ca_TnCHc = " << dBuff_Ca_TnCHc << ", dBuff_Ca_TnCHm = " << dBuff_Ca_TnCHm << ", dBuff_Ca_CaM = " << dBuff_Ca_CaM << ", dBuff_Ca_Myosin_ca = " << dBuff_Ca_Myosin_ca << ", dBuff_Ca_Myosin_mg = " << dBuff_Ca_Myosin_mg << endl; 	
		//cout << "dBuff_Ca_SRB = " << dBuff_Ca_SRB << ", dBuff_Ca_EGTA_cyt = " << dBuff_Ca_EGTA_cyt << endl;
 		double dCa_i, dCa_srs;
		dCa_srs = 0; dCa_i = 0;
		
		if(ca_domain == 0) 
		{
			//Leftmost unit, no release, but diffusion from membrane 0
			dCa_i = -segments[segment].caUnits[ca_domain].Jup * data->vsr/data->vmyo - J_CaB_cytosol + data->J_ca_slmyo*0.5/(data->vmyo/settings->N_CaDomains)*(segments[segment].membrane[0].Ca_sl - segments[segment].caUnits[ca_domain].Ca_i) + J_diff_i_subcell + J_diff_SRS_i;
			
			double I_Ca_tot_junc = segments[segment].membrane[0].ICaL_junc + segments[segment].membrane[0].ICaT_junc + segments[segment].membrane[0].ICab_junc - 2*segments[segment].membrane[0].INaCa_junc;		// [uA/uF]
			dCa_srs = -I_Ca_tot_junc * data->CmemF/(data->vjunc*2) + data->J_ca_juncsl*0.5/(data->vjunc/settings->N_CaDomains)*(segments[segment].membrane[0].Ca_sl - segments[segment].caUnits[ca_domain].Ca_srs) + segments[segment].caUnits[ca_domain].Jrel*data->vsr/data->vjunc + segments[segment].caUnits[ca_domain].Jleak*data->vmyo/data->vjunc + J_diff_srs_subcell - data->vmyo/data->vjunc * J_diff_SRS_i;
		}
		else
		{
			if(ca_domain == settings->N_CaDomains - 1)
			{
				//Rightmost unit, no release, but diffusion from membrane 1
				dCa_i = -segments[segment].caUnits[ca_domain].Jup * data->vsr/data->vmyo - J_CaB_cytosol + data->J_ca_slmyo*0.5/(data->vmyo/settings->N_CaDomains)*(segments[segment].membrane[1].Ca_sl - segments[segment].caUnits[ca_domain].Ca_i) + J_diff_i_subcell + J_diff_SRS_i;

				double I_Ca_tot_junc = segments[segment].membrane[1].ICaL_junc + segments[segment].membrane[1].ICaT_junc + segments[segment].membrane[1].ICab_junc  - 2*segments[segment].membrane[1].INaCa_junc;		// [uA/uF]
				dCa_srs = -I_Ca_tot_junc * data->CmemF/(data->vjunc*2) + data->J_ca_juncsl*0.5/(data->vjunc/settings->N_CaDomains)*(segments[segment].membrane[1].Ca_sl - segments[segment].caUnits[ca_domain].Ca_srs) + segments[segment].caUnits[ca_domain].Jrel*data->vsr/data->vjunc + segments[segment].caUnits[ca_domain].Jleak*data->vmyo/data->vjunc + J_diff_srs_subcell - data->vmyo/data->vjunc * J_diff_SRS_i;
			}
			else
			{
				//Other unit, release from corbular SR, no diffusion from membrane
				dCa_i = -segments[segment].caUnits[ca_domain].Jup * data->vsr/data->vmyo - J_CaB_cytosol + J_diff_i_subcell + J_diff_SRS_i;
				dCa_srs = segments[segment].caUnits[ca_domain].Jrel*data->vsr/data->vjunc + segments[segment].caUnits[ca_domain].Jleak*data->vmyo/data->vjunc + J_diff_srs_subcell - data->vmyo/data->vjunc * J_diff_SRS_i;
			}
		}
		//cout << "Ca_sl = " << Ca_sl << ", dCa_i = " << dCa_i << ", J_diff_i_subcell = " << J_diff_i_subcell << ", J_diff_SRS_i = " << J_diff_SRS_i  << endl;
		//Update states
		double cor_fact = 1.0;
		segments[segment].caUnits[ca_domain].Ca_i = check_State_Change(segments[segment].caUnits[ca_domain].Ca_i, dCa_i, dt);	
				
		segments[segment].caUnits[ca_domain].J_diff_i_subcell = J_diff_i_subcell;
		
		double Km_SLL = data->koff_sll / data->kon_sll;
		double Km_SLH = data->koff_slh / data->kon_slh;
			
		double Ca_srs_tot = segments[segment].caUnits[ca_domain].Ca_srs + segments[segment].caUnits[ca_domain].Buff_Ca_SLLsrs + segments[segment].caUnits[ca_domain].Buff_Ca_SLHsrs;
		Ca_srs_tot = check_State_Change(Ca_srs_tot, dCa_srs, dt);
		segments[segment].caUnits[ca_domain].Ca_srs = Buff_Ca_DoubleBuff(Ca_srs_tot, data->Bmax_SLlowj, data->Bmax_SLhighj, Km_SLL, Km_SLH);
		segments[segment].caUnits[ca_domain].Buff_Ca_SLLsrs = segments[segment].caUnits[ca_domain].Ca_srs * data->Bmax_SLlowj / (segments[segment].caUnits[ca_domain].Ca_srs + Km_SLL);
		segments[segment].caUnits[ca_domain].Buff_Ca_SLHsrs = Ca_srs_tot - segments[segment].caUnits[ca_domain].Ca_srs - segments[segment].caUnits[ca_domain].Buff_Ca_SLLsrs;
		
		double Ca_sr_tot = segments[segment].caUnits[ca_domain].Ca_sr + segments[segment].caUnits[ca_domain].Buff_Ca_CSQN;
		double dCa_sr_tot = segments[segment].caUnits[ca_domain].Jup - (segments[segment].caUnits[ca_domain].Jleak * data->vmyo/data->vsr + segments[segment].caUnits[ca_domain].Jrel) + J_diff_sr_subcell;         // Ca_sr     [mM/ms] %Ratio 3 leak current
		Ca_sr_tot = check_State_Change(Ca_sr_tot, dCa_sr_tot, dt);
		double b_Casr_buff = Ca_sr_tot + data->Bmax_Csqn + data->koff_csqn / data->kon_csqn;
		double c_Casr_buff = Ca_sr_tot * data->Bmax_Csqn;
		segments[segment].caUnits[ca_domain].Buff_Ca_CSQN = 0.5 * (b_Casr_buff - sqrt(b_Casr_buff*b_Casr_buff - 4*c_Casr_buff));
		segments[segment].caUnits[ca_domain].Ca_sr = Ca_sr_tot - segments[segment].caUnits[ca_domain].Buff_Ca_CSQN;
				
		segments[segment].caUnits[ca_domain].Buff_Ca_TnCL = check_State_Change(segments[segment].caUnits[ca_domain].Buff_Ca_TnCL, dBuff_Ca_TnCL*cor_fact, dt);
		segments[segment].caUnits[ca_domain].Buff_Ca_TnCHc = check_State_Change(segments[segment].caUnits[ca_domain].Buff_Ca_TnCHc, dBuff_Ca_TnCHc*cor_fact, dt);
		segments[segment].caUnits[ca_domain].Buff_Ca_TnCHm = check_State_Change(segments[segment].caUnits[ca_domain].Buff_Ca_TnCHm, dBuff_Ca_TnCHm*cor_fact, dt);
		segments[segment].caUnits[ca_domain].Buff_Ca_CaM = check_State_Change(segments[segment].caUnits[ca_domain].Buff_Ca_CaM, dBuff_Ca_CaM*cor_fact, dt);
		segments[segment].caUnits[ca_domain].Buff_Ca_Myosin_ca = check_State_Change(segments[segment].caUnits[ca_domain].Buff_Ca_Myosin_ca, dBuff_Ca_Myosin_ca*cor_fact, dt);
		segments[segment].caUnits[ca_domain].Buff_Ca_Myosin_mg = check_State_Change(segments[segment].caUnits[ca_domain].Buff_Ca_Myosin_mg, dBuff_Ca_Myosin_mg*cor_fact, dt);
		segments[segment].caUnits[ca_domain].Buff_Ca_SRB = check_State_Change(segments[segment].caUnits[ca_domain].Buff_Ca_SRB, dBuff_Ca_SRB*cor_fact, dt);
		segments[segment].caUnits[ca_domain].Buff_Ca_EGTA_cyt = check_State_Change(segments[segment].caUnits[ca_domain].Buff_Ca_EGTA_cyt, dBuff_Ca_EGTA_cyt*cor_fact, dt);
	}

}

void RAMM::UpdateConcentrations(int seg_i, int mem_i, double dt, double t, double Na_sl_U_Segment, double Na_sl_D_Segment, double Ca_sl_U_Segment, double Ca_sl_D_Segment, double Ca_i_neighbor, double Ca_srs_neighbor)
{
	
	// Sodium Buffering and Concentrations
	double dBuff_Na_junc = data->kon_na * segments[seg_i].membrane[mem_i].Na_junc * (data->Bmax_Naj-segments[seg_i].membrane[mem_i].Buff_Na_junc) - data->koff_na*segments[seg_i].membrane[mem_i].Buff_Na_junc; // NaBj      [mM/ms]
	double dBuff_Na_sl = data->kon_na * segments[seg_i].membrane[mem_i].Na_sl * (data->Bmax_Nasl-segments[seg_i].membrane[mem_i].Buff_Na_sl) - data->koff_na*segments[seg_i].membrane[mem_i].Buff_Na_sl;        // NaBsl     [mM/ms]
	
	double INa_tot_junc = segments[seg_i].membrane[mem_i].INa_junc + segments[seg_i].membrane[mem_i].INab_junc + 3*segments[seg_i].membrane[mem_i].INaCa_junc + 3*segments[seg_i].membrane[mem_i].INaK_junc;
	// + segments[seg_i].membrane[mem_i].ICaL_Na_junc;
	double INa_tot_sl = segments[seg_i].membrane[mem_i].INa_sl + segments[seg_i].membrane[mem_i].INab_sl + 3*segments[seg_i].membrane[mem_i].INaCa_sl + 3*segments[seg_i].membrane[mem_i].INaK_sl;
	// + segments[seg_i].membrane[mem_i].ICaL_Na_sl;	
			
	double J_diff_sl_subcell = (Na_sl_U_Segment - segments[seg_i].membrane[mem_i].Na_sl + Na_sl_D_Segment - segments[seg_i].membrane[mem_i].Na_sl) / data->tau_diff_Casl_Segment;
	double dNa_junc = -INa_tot_junc * data->CmemF/data->vjunc + data->J_na_juncsl/data->vjunc * (segments[seg_i].membrane[mem_i].Na_sl - segments[seg_i].membrane[mem_i].Na_junc) - dBuff_Na_junc;
	double dNa_sl = -INa_tot_sl * data->CmemF/data->vsl + data->J_na_juncsl/data->vsl * (segments[seg_i].membrane[mem_i].Na_junc - segments[seg_i].membrane[mem_i].Na_sl) + data->J_na_slmyo/(0.5 * data->vsl) * (segments[seg_i].Na_i - segments[seg_i].membrane[mem_i].Na_sl) - dBuff_Na_sl + J_diff_sl_subcell;
	
	segments[seg_i].membrane[mem_i].Buff_Na_junc = check_State_Change(segments[seg_i].membrane[mem_i].Buff_Na_junc, dBuff_Na_junc, dt);
	segments[seg_i].membrane[mem_i].Buff_Na_sl = check_State_Change(segments[seg_i].membrane[mem_i].Buff_Na_sl, dBuff_Na_sl, dt);
	segments[seg_i].membrane[mem_i].Na_junc = check_State_Change(segments[seg_i].membrane[mem_i].Na_junc, dNa_junc, dt);
	segments[seg_i].membrane[mem_i].Na_sl = check_State_Change(segments[seg_i].membrane[mem_i].Na_sl, dNa_sl, dt);
	
	double I_Ca_tot_junc = segments[seg_i].membrane[mem_i].ICaL_junc + segments[seg_i].membrane[mem_i].ICab_junc + segments[seg_i].membrane[mem_i].ICaT_junc - 2*segments[seg_i].membrane[mem_i].INaCa_junc;// [uA/uF]
	double I_Ca_tot_sl = segments[seg_i].membrane[mem_i].ICaL_sl + segments[seg_i].membrane[mem_i].ICaT_sl + segments[seg_i].membrane[mem_i].ICab_sl + segments[seg_i].membrane[mem_i].IpCa_sl - 2*segments[seg_i].membrane[mem_i].INaCa_sl;            		// [uA/uF]
	double J_diff_Casl_subcell = (Ca_sl_U_Segment - segments[seg_i].membrane[mem_i].Ca_sl + Ca_sl_D_Segment - segments[seg_i].membrane[mem_i].Ca_sl) / data->tau_diff_Casl_Segment;

	double dBuff_Ca_EGTA_sl = (0.005 / dt) * (data->kon_EGTA * segments[seg_i].membrane[mem_i].Ca_sl * (data->EGTA * (data->vmyo / data->vsl) - segments[seg_i].membrane[mem_i].Buff_Ca_EGTA_sl) - data->koff_EGTA * segments[seg_i].membrane[mem_i].Buff_Ca_EGTA_sl);
	double Km_SLL = data->koff_sll / data->kon_sll;
	double Km_SLH = data->koff_slh / data->kon_slh;
		
	double Ca_junc_tot = segments[seg_i].membrane[mem_i].Ca_junc + segments[seg_i].membrane[mem_i].Buff_Ca_SLLj + segments[seg_i].membrane[mem_i].Buff_Ca_SLHj;
	double Ca_sl_tot = segments[seg_i].membrane[mem_i].Ca_sl + segments[seg_i].membrane[mem_i].Buff_Ca_SLLsl + segments[seg_i].membrane[mem_i].Buff_Ca_SLHsl;
	
	double dCa_junc_tot, dCa_sl_tot;
	if(mem_i == 0) dCa_junc_tot = data->J_ca_juncsl/data->vjunc*(segments[seg_i].membrane[mem_i].Ca_sl - segments[seg_i].membrane[mem_i].Ca_junc) + (segments[seg_i].caUnits[0].Ca_srs - segments[seg_i].membrane[mem_i].Ca_junc) / data->tau_diff_Casrs_Domain;  // Ca_j
	if(mem_i == 1) dCa_junc_tot = data->J_ca_juncsl/data->vjunc*(segments[seg_i].membrane[mem_i].Ca_sl - segments[seg_i].membrane[mem_i].Ca_junc) + (segments[seg_i].caUnits[settings->N_CaDomains-1].Ca_srs - segments[seg_i].membrane[mem_i].Ca_junc) / data->tau_diff_Casrs_Domain;  // Ca_j		
	dCa_sl_tot = -I_Ca_tot_sl * data->CmemF/(data->vsl*2) + data->J_ca_juncsl/data->vsl*(segments[seg_i].membrane[mem_i].Ca_junc - segments[seg_i].membrane[mem_i].Ca_sl) + data->J_ca_slmyo/(data->vsl)*(Ca_i_neighbor - segments[seg_i].membrane[mem_i].Ca_sl) + J_diff_Casl_subcell - dBuff_Ca_EGTA_sl;   // Ca_sl
		
	dCa_junc_tot = 0;			
	dCa_sl_tot = -I_Ca_tot_sl * data->CmemF/(data->vsl*2) + data->J_ca_juncsl/data->vsl*(Ca_srs_neighbor - segments[seg_i].membrane[mem_i].Ca_sl) + data->J_ca_slmyo/(data->vsl)*(Ca_i_neighbor - segments[seg_i].membrane[mem_i].Ca_sl) + J_diff_Casl_subcell - dBuff_Ca_EGTA_sl;   // Ca_sl						

	Ca_junc_tot = check_State_Change(Ca_junc_tot, dCa_junc_tot, dt);
	segments[seg_i].membrane[mem_i].Ca_junc = Buff_Ca_DoubleBuff(Ca_junc_tot, data->Bmax_SLlowj, data->Bmax_SLhighj, Km_SLL, Km_SLH);
	segments[seg_i].membrane[mem_i].Buff_Ca_SLLj = segments[seg_i].membrane[mem_i].Ca_junc * data->Bmax_SLlowj / (segments[seg_i].membrane[mem_i].Ca_junc + Km_SLL);
	segments[seg_i].membrane[mem_i].Buff_Ca_SLHj = Ca_junc_tot - segments[seg_i].membrane[mem_i].Ca_junc - segments[seg_i].membrane[mem_i].Buff_Ca_SLLj;
				
	//cout << "I_Ca_tot_sl = " << I_Ca_tot_sl << ", dCa_sl_tot = " << dCa_sl_tot  << ", J_diff_Casl_subcell = " << J_diff_Casl_subcell << ", dBuff_Ca_EGTA_sl = " << dBuff_Ca_EGTA_sl << endl;
	Ca_sl_tot = check_State_Change(Ca_sl_tot, dCa_sl_tot, dt);
	segments[seg_i].membrane[mem_i].Ca_sl = Buff_Ca_DoubleBuff(Ca_sl_tot, data->Bmax_SLlowsl, data->Bmax_SLhighsl, Km_SLL, Km_SLH);
	segments[seg_i].membrane[mem_i].Buff_Ca_SLLsl = segments[seg_i].membrane[mem_i].Ca_sl * data->Bmax_SLlowsl / (segments[seg_i].membrane[mem_i].Ca_sl + Km_SLL);
	segments[seg_i].membrane[mem_i].Buff_Ca_SLHsl = Ca_sl_tot - segments[seg_i].membrane[mem_i].Ca_sl - segments[seg_i].membrane[mem_i].Buff_Ca_SLLsl;
	segments[seg_i].membrane[mem_i].Buff_Ca_EGTA_sl = check_State_Change(segments[seg_i].membrane[mem_i].Buff_Ca_EGTA_sl, dBuff_Ca_EGTA_sl, dt);

}

//Stochastic Ryanodine Receptor formulation
void RAMM::Ca_CICR_Stoch(int seg_i, int cadom_i, double dt, double t)
{
	
	double kCaSR;
	
	double Ca_cis_observed = segments[seg_i].caUnits[cadom_i].Ca_i;
	double Ca_sr_observed = segments[seg_i].caUnits[cadom_i].Ca_sr;
	
	Ca_cis_observed = segments[seg_i].caUnits[cadom_i].Ca_srs;
	segments[seg_i].caUnits[cadom_i].Jrel = (1 - settings->IrelB) * settings->RyRParams_NP[11] * settings->N_CaDomains * segments[seg_i].caUnits[cadom_i].RyRState_NP[1] / ((1.0 * settings->numRyRs) / (settings->N_Segments * settings->N_CaDomains)) * (Ca_sr_observed-Ca_cis_observed);          // [mM/ms]
	segments[seg_i].caUnits[cadom_i].Jleak = 0; //(1 - settings->IrelB) * Jleakmax * (segments[seg_i].caUnits[cadom_i].Ca_sr-Ca_cis_observed);           //   [mM/ms]	
	
	if (settings->RyRParams_NP[21] > 0 && settings->RyRParams_NP[21] < 100 && Ca_sr_observed > settings->RyRParams_NP[21]) Ca_sr_observed = settings->RyRParams_NP[21];
	
	// RyR Gating
	double exp_Ca_cis, Act_ss, tau_Act, Act_ss_low;
	
	double deltaf = 1; // + 2 * (settings->RyRParams_NP[19] - 1) / (1 + pow(Ca_sr_observed / settings->RyRParams_NP[1], settings->RyRParams_NP[2]));	
	exp_Ca_cis = 1.0 / (1 + exp(-(Ca_cis_observed - settings->RyRParams_NP[4]*deltaf)/settings->RyRParams_NP[5]));
	Act_ss = (1 / (1 + pow(segments[seg_i].caUnits[cadom_i].RyR_Integr / 0.5, 8))) * (settings->RyRParams_NP[6] + settings->RyRParams_NP[3] * exp_Ca_cis);	
	Act_ss_low = (1 / (1 + pow(segments[seg_i].caUnits[cadom_i].RyR_Integr / 0.5, 8))) * (settings->RyRParams_NP[6] + settings->RyRParams_NP[3] * (1.0 / (1 + exp(-(Ca_cis_observed - settings->RyRParams_NP[20] * settings->RyRParams_NP[4]*deltaf)/settings->RyRParams_NP[5]))));	
	tau_Act = settings->RyRParams_NP[7] + settings->RyRParams_NP[8] * (1 - 1.0 / (1 + exp(-(Ca_cis_observed - settings->RyRParams_NP[18])/settings->RyRParams_NP[5])));
		
	double midlev = settings->RyRParams_NP[0]; // + (1 - settings->RyRParams_NP[0]) / (1 + pow(0.15 / segments[seg_i].caUnits[cadom_i].RyR_Integr, 8));
	double Inact_ss = midlev / (1 + pow(Ca_sr_observed / 4,4)) + (settings->RyRParams_NP[9] - midlev) / (1 + pow(Ca_sr_observed / settings->RyRParams_NP[1], settings->RyRParams_NP[2]));
	double tau_Inact = settings->RyRParams_NP[10] + settings->RyRParams_NP[12] / (1 + pow(settings->RyRParams_NP[1] / Ca_sr_observed, settings->RyRParams_NP[2]));
	
	double tau_Flec = settings->RyRParams_NP[24];
	double Flec_Block_ss = 0;
	if(settings->Flec > 0)
	{
		Flec_Block_ss = 1.0 / (1 + pow(settings->RyRParams_NP[22] / settings->Flec, settings->RyRParams_NP[23]));
	}
	
	//Methodology to clamp RyR2
	if(settings->clamp_RyR_Po >= -0.1 && settings->clamp_RyR_Po <= 1) 
	{
		if(settings->numRyR_Po_Indices < 0 || GenPurpose::occursInList(settings->RyR_Po_Indices, settings->numRyR_Po_Indices, (seg_i - 1) * settings->N_CaDomains + cadom_i))
		{
			if(settings->clamp_RyR_Po >= 0)
			{
				Act_ss = settings->clamp_RyR_Po;
			
				if(settings->clamp_RyR_Po >= 0.3) 
				{
					Inact_ss = 0;
					Flec_Block_ss = 0;
					segments[seg_i].caUnits[cadom_i].Jrel = segments[seg_i].caUnits[cadom_i].Jrel * 2;
				}	
				exp_Ca_cis = 1.0 / (1 + exp(-(Ca_cis_observed - settings->RyRParams_NP[4])/settings->RyRParams_NP[5]));
				tau_Act = settings->RyRParams_NP[7] + settings->RyRParams_NP[8] * kCaSR * (1 - exp_Ca_cis);
			}
			else
			{
				Act_ss = (-10 * settings->clamp_RyR_Po) * Act_ss;
			}
		}
	}
	
	double k_12 = Act_ss / tau_Act;
	double k_21 = (1 - Act_ss) / tau_Act;
	double k_43 = Act_ss_low / tau_Act;
	double k_34 = (1 - Act_ss_low) / tau_Act;
	double k_23 = Inact_ss / tau_Inact;
	double k_32 = (1 - Inact_ss) / tau_Inact;
	double k_41 = 0.1 * k_32;
	double k_14 = (k_12 * k_23 * k_34 * k_41) / (k_43 * k_32 * k_21); 
	
	double k_26 = Flec_Block_ss / tau_Flec;
	double k_62 = (1 - Flec_Block_ss) / tau_Flec;
	
	// Used to ensure stepsize adaptation
	double RyRStateSum = segments[seg_i].caUnits[cadom_i].RyRState_NP[0] + segments[seg_i].caUnits[cadom_i].RyRState_NP[1] + segments[seg_i].caUnits[cadom_i].RyRState_NP[2] + segments[seg_i].caUnits[cadom_i].RyRState_NP[3];
	double dR = ((k_41*segments[seg_i].caUnits[cadom_i].RyRState_NP[3]-k_14*segments[seg_i].caUnits[cadom_i].RyRState_NP[0])-(k_12*segments[seg_i].caUnits[cadom_i].RyRState_NP[0]-k_21*segments[seg_i].caUnits[cadom_i].RyRState_NP[1])) / RyRStateSum;   // R
	double dO = ((k_12*segments[seg_i].caUnits[cadom_i].RyRState_NP[0]-k_21*segments[seg_i].caUnits[cadom_i].RyRState_NP[1])-(k_23*segments[seg_i].caUnits[cadom_i].RyRState_NP[1]-k_32*segments[seg_i].caUnits[cadom_i].RyRState_NP[2])) / RyRStateSum; // O
	double dI = ((k_23*segments[seg_i].caUnits[cadom_i].RyRState_NP[1]-k_32*segments[seg_i].caUnits[cadom_i].RyRState_NP[2])-(k_34*segments[seg_i].caUnits[cadom_i].RyRState_NP[2]-k_43*segments[seg_i].caUnits[cadom_i].RyRState_NP[3])) / RyRStateSum;   // I
		check_State_Change(segments[seg_i].caUnits[cadom_i].RyRState_NP[0] / RyRStateSum, dR, dt);	
		check_State_Change(segments[seg_i].caUnits[cadom_i].RyRState_NP[1] / RyRStateSum, dO, dt);
		check_State_Change(segments[seg_i].caUnits[cadom_i].RyRState_NP[2] / RyRStateSum, dI, dt);
	
		// Build probability matrix
		// For methodology see Heijman et al. PLoS Comput Biol 2013.
		double RyR_Matrix[8][8];
		double dt_used = dt;
			RyR_Matrix[0][0] = 1 - dt_used * (k_12 + k_14); RyR_Matrix[0][1] = dt_used * k_12; RyR_Matrix[0][2] = 0; RyR_Matrix[0][3] = dt_used * k_14; RyR_Matrix[0][4] = 0; RyR_Matrix[0][5] = 0; RyR_Matrix[0][6] = 0; RyR_Matrix[0][7] = 0;
			RyR_Matrix[1][0] = dt_used * k_21; RyR_Matrix[1][1] = 1 - dt_used * (k_21+k_23+k_26); RyR_Matrix[1][2] = dt_used * k_23; RyR_Matrix[1][3] = 0; RyR_Matrix[1][4] = 0; RyR_Matrix[1][5] = dt_used * k_26; RyR_Matrix[1][6] = 0; RyR_Matrix[1][7] = 0;
			RyR_Matrix[2][0] = 0; RyR_Matrix[2][1] = dt_used * k_32; RyR_Matrix[2][2] = 1 - dt_used * (k_32 + k_34); RyR_Matrix[2][3] = dt_used * k_34; RyR_Matrix[2][4] = 0; RyR_Matrix[2][5] = 0; RyR_Matrix[2][6] = 0; RyR_Matrix[2][7] = 0;
			RyR_Matrix[3][0] = dt_used * k_41; RyR_Matrix[3][1] = 0; RyR_Matrix[3][2] = dt_used * k_43; RyR_Matrix[3][3] = 1 - dt_used * (k_41 + k_43); RyR_Matrix[3][4] = 0; RyR_Matrix[3][5] = 0; RyR_Matrix[3][6] = 0; RyR_Matrix[3][7] = 0;
		
			RyR_Matrix[4][4] = 1 - dt_used * (k_12 + k_14); RyR_Matrix[4][5] = dt_used * k_12; RyR_Matrix[4][6] = 0; RyR_Matrix[4][7] = dt_used * k_14; RyR_Matrix[4][0] = 0; RyR_Matrix[4][1] = 0; RyR_Matrix[4][2] = 0; RyR_Matrix[4][3] = 0;
			RyR_Matrix[5][4] = dt_used * k_21; RyR_Matrix[5][5] = 1 - dt_used * (k_21+k_23+k_62); RyR_Matrix[5][6] = dt_used * k_23; RyR_Matrix[5][7] = 0; RyR_Matrix[5][0] = 0; RyR_Matrix[5][1] = dt_used * k_62; RyR_Matrix[5][2] = 0; RyR_Matrix[5][3] = 0;
			RyR_Matrix[6][4] = 0; RyR_Matrix[6][5] = dt_used * k_32; RyR_Matrix[6][6] = 1 - dt_used * (k_32 + k_34); RyR_Matrix[6][7] = dt_used * k_34; RyR_Matrix[6][0] = 0; RyR_Matrix[6][1] = 0; RyR_Matrix[6][2] = 0; RyR_Matrix[6][3] = 0;
			RyR_Matrix[7][4] = dt_used * k_41; RyR_Matrix[7][5] = 0; RyR_Matrix[7][6] = dt_used * k_43; RyR_Matrix[7][7] = 1 - dt_used * (k_41 + k_43); RyR_Matrix[7][0] = 0; RyR_Matrix[7][1] = 0; RyR_Matrix[7][2] = 0; RyR_Matrix[7][3] = 0;			
		
		
		// Compute state transitions		
			const int nrOfStatesInModel = 8;
			int newRyRState[nrOfStatesInModel];		
			for(int i=0; i<nrOfStatesInModel; i++) { newRyRState[i] = 0;	}
					
			for(int i=0; i<nrOfStatesInModel; i++)
			{
				if(segments[seg_i].caUnits[cadom_i].RyRState_NP[i] >= 1 && segments[seg_i].caUnits[cadom_i].RyRState_NP[i] <= 10) //Use uniform dist approach
				{
					for(int k=0; k<segments[seg_i].caUnits[cadom_i].RyRState_NP[i]; k++)
					{
						double cumulprob = 0.0;
						double rv = rvgen.rand();
						
						for(int j=0; j<nrOfStatesInModel; j++)
						{
							cumulprob += RyR_Matrix[i][j];
							if(cumulprob > rv)
							{
								newRyRState[j]++;
								break;
							}
						}
					}
				}
				else
				{
					if(segments[seg_i].caUnits[cadom_i].RyRState_NP[i] > 10) // Use binomial approach				
					{
						int numchannels = (int)segments[seg_i].caUnits[cadom_i].RyRState_NP[i];
						for(int dd=0; dd<nrOfStatesInModel; dd++)
						{
							if(dd != i && RyR_Matrix[i][dd] > 0) //Only take states other than self with nonzero chance
							{
								int numchannelsout = GenPurpose::applyBinomialInverseTransform(&rvgen, (int)segments[seg_i].caUnits[cadom_i].RyRState_NP[i], RyR_Matrix[i][dd]);
								if(numchannelsout > numchannels) numchannelsout = numchannels;
								newRyRState[dd] += numchannelsout;
								numchannels -= numchannelsout;
							}
						}
						newRyRState[i] += numchannels; // The remaining channels stay where they are
					}
				}		
			}

			//Ability to set channel state to a specific configuration
			if(settings->set_RyR_State != NULL) 
			{
				int dom_index = (seg_i * settings->N_CaDomains + cadom_i);
				double* curDomRyRState = settings->set_RyR_State[dom_index];
				
				if(curDomRyRState[0] >= 0)
				{
					newRyRState[0] = curDomRyRState[0];
					newRyRState[1] = curDomRyRState[1];
					newRyRState[2] = curDomRyRState[2];
					newRyRState[3] = curDomRyRState[3];
				}
			}
			
			for(int i=0; i<nrOfStatesInModel; i++) { segments[seg_i].caUnits[cadom_i].RyRState_NP[i] = newRyRState[i];	}
		
	
	double RyR_SEP_SS = 1.0 / (1.0 + exp((Ca_cis_observed - 1.8E-3) / 1E-6)); //1 / (1 + pow(Ca_cis_observed / 2.75E-4, 20));
	double tau_RyR_SEP = 60;
	segments[seg_i].caUnits[cadom_i].RyR_SEP = check_State_Change(segments[seg_i].caUnits[cadom_i].RyR_SEP, (RyR_SEP_SS - segments[seg_i].caUnits[cadom_i].RyR_SEP) / tau_RyR_SEP, dt);

	//double dRyR_Int = 1 * segments[seg_i].caUnits[cadom_i].Jrel - 0.01 * segments[seg_i].caUnits[cadom_i].RyR_Integr;
	double dRyR_Int = 1 * segments[seg_i].caUnits[cadom_i].Jrel - 0.1 / (1 + exp((segments[seg_i].caUnits[cadom_i].Jrel - 1.5E-3) / 1E-4)) * segments[seg_i].caUnits[cadom_i].RyR_Integr;
	segments[seg_i].caUnits[cadom_i].RyR_Integr = check_State_Change(segments[seg_i].caUnits[cadom_i].RyR_Integr, dRyR_Int, dt);
	
	// SERCA
	double tval_Kmr = pow(segments[seg_i].caUnits[cadom_i].Ca_sr / data->Kmr, data->hillSRCaP);
	double tval_Kmf_NP = pow(segments[seg_i].caUnits[cadom_i].Ca_i / data->Kmf_NP, data->hillSRCaP);
	double Jup_NP = (1 - settings->IupB) * data->Vmax_SRCaP * (tval_Kmf_NP - tval_Kmr) / (1 + tval_Kmf_NP + tval_Kmr);
	segments[seg_i].caUnits[cadom_i].Jup = Jup_NP;
	
}

// Deterministic RyR2 simulation
void RAMM::Ca_CICR(int seg_i, int cadom_i, double dt, double t)
{
	
	double kCaSR;
		
	double Ca_cis_observed = segments[seg_i].caUnits[cadom_i].Ca_srs;
	if(segments[seg_i].caUnits[cadom_i].RyRState_NP[1] > 2 || segments[seg_i].caUnits[cadom_i].RyRState_NP[0] > 2)
	{
		cout << "Stochastic RyR state being used in deterministic simulations!! - t = " << t << ", dt = " << dt << ", seg_i = " << seg_i << ", cadom_i = " << cadom_i << endl;
		GenPurpose::wait(2.0);
	}
	
	segments[seg_i].caUnits[cadom_i].Jrel = (1 - settings->IrelB) * settings->RyRParams_NP[11] * settings->N_CaDomains * segments[seg_i].caUnits[cadom_i].RyRState_NP[1] * (segments[seg_i].caUnits[cadom_i].Ca_sr-Ca_cis_observed);          // [mM/ms]
	segments[seg_i].caUnits[cadom_i].Jleak = 0; //(1 - settings->IrelB) * Jleakmax * (segments[seg_i].caUnits[cadom_i].Ca_sr-Ca_cis_observed);           //   [mM/ms]

	double k_12, k_21, k_14, k_41;
	double exp_Ca_cis, Act_ss, tau_Act, Act_ss_low;
	double Inact_ss, tau_Inact;

	double Ca_sr_observed = segments[seg_i].caUnits[cadom_i].Ca_sr;
	if (settings->RyRParams_NP[21] > 0 && settings->RyRParams_NP[21] < 100 && Ca_sr_observed > settings->RyRParams_NP[21]) Ca_sr_observed = settings->RyRParams_NP[21];
		
	double deltaf = 1; //1 + 2 * (settings->RyRParams_NP[19] - 1) / (1 + pow(Ca_sr_observed / settings->RyRParams_NP[1], settings->RyRParams_NP[2]));
	exp_Ca_cis = 1.0 / (1 + exp(-(Ca_cis_observed - settings->RyRParams_NP[4]*deltaf)/settings->RyRParams_NP[5]));
	Act_ss =  (1 / (1 + pow(segments[seg_i].caUnits[cadom_i].RyR_Integr / 0.5, 8))) * (settings->RyRParams_NP[6] + settings->RyRParams_NP[3] * exp_Ca_cis);		
	Act_ss_low =  (1 / (1 + pow(segments[seg_i].caUnits[cadom_i].RyR_Integr / 0.5, 8))) * (settings->RyRParams_NP[6] + settings->RyRParams_NP[3] * (1.0 / (1 + exp(-(Ca_cis_observed - settings->RyRParams_NP[20] * settings->RyRParams_NP[4]*deltaf)/settings->RyRParams_NP[5]))));
	tau_Act = settings->RyRParams_NP[7] + settings->RyRParams_NP[8] * (1 - 1.0 / (1 + exp(-(Ca_cis_observed - settings->RyRParams_NP[18])/settings->RyRParams_NP[5])));
	double midlev = settings->RyRParams_NP[0]; // + (1 - settings->RyRParams_NP[0]) / (1 + pow(0.008 / Ca_cis_observed, 6));
	Inact_ss = (midlev / (1 + pow(Ca_sr_observed / 4,4)) + (settings->RyRParams_NP[9] - midlev) / (1 + pow(Ca_sr_observed / settings->RyRParams_NP[1], settings->RyRParams_NP[2])));
	tau_Inact = settings->RyRParams_NP[10] + settings->RyRParams_NP[12] / (1 + pow(settings->RyRParams_NP[1] / Ca_sr_observed, settings->RyRParams_NP[2]));
	
	double tau_Flec = settings->RyRParams_NP[24];
	double Flec_Block_ss = 0;
	if(settings->Flec > 0)
	{
		Flec_Block_ss = 1.0 / (1 + pow(settings->RyRParams_NP[22] / settings->Flec, settings->RyRParams_NP[23]));
	}
	
	//Methodology to clamp RyR2
	if(settings->clamp_RyR_Po >= -0.1 && settings->clamp_RyR_Po <= 1) 
	{
		if(settings->numRyR_Po_Indices < 0 || GenPurpose::occursInList(settings->RyR_Po_Indices, settings->numRyR_Po_Indices, (seg_i - 1) * settings->N_CaDomains + cadom_i))
		{
			if(settings->clamp_RyR_Po >= 0)
			{
				Act_ss = settings->clamp_RyR_Po;
			
				if(settings->clamp_RyR_Po >= 0.3) 
				{
					Inact_ss = 0;
					Flec_Block_ss = 0;
					segments[seg_i].caUnits[cadom_i].Jrel = segments[seg_i].caUnits[cadom_i].Jrel * 2;
				}	
				exp_Ca_cis = 1.0 / (1 + exp(-(Ca_cis_observed - settings->RyRParams_NP[4])/settings->RyRParams_NP[5]));
				tau_Act = settings->RyRParams_NP[7] + settings->RyRParams_NP[8] * kCaSR * (1 - exp_Ca_cis);
			}
			else
			{
				Act_ss = (-10 * settings->clamp_RyR_Po) * Act_ss;
			}
		}
	}

	k_12 = Act_ss / tau_Act;
	k_21 = (1 - Act_ss) / tau_Act;
	double k_43 = Act_ss_low / tau_Act;
	double k_34 = (1 - Act_ss_low) / tau_Act;
	double k_23 = Inact_ss / tau_Inact;
	double k_32 = (1 - Inact_ss) / tau_Inact;
	k_41 = 0.1 * k_32;
	k_14 = (k_12 * k_23 * k_34 * k_41) / (k_43 * k_32 * k_21);
	
	double k_26 = Flec_Block_ss / tau_Flec;
	double k_62 = (1 - Flec_Block_ss) / tau_Flec;
	
	double dR = (k_41*segments[seg_i].caUnits[cadom_i].RyRState_NP[3]-k_14*segments[seg_i].caUnits[cadom_i].RyRState_NP[0])-(k_12*segments[seg_i].caUnits[cadom_i].RyRState_NP[0]-k_21*segments[seg_i].caUnits[cadom_i].RyRState_NP[1]);   // R
	double dO = (k_12*segments[seg_i].caUnits[cadom_i].RyRState_NP[0]-k_21*segments[seg_i].caUnits[cadom_i].RyRState_NP[1])-(k_23*segments[seg_i].caUnits[cadom_i].RyRState_NP[1]-k_32*segments[seg_i].caUnits[cadom_i].RyRState_NP[2]) + k_62 * segments[seg_i].caUnits[cadom_i].RyRState_NP[5] - k_26 * segments[seg_i].caUnits[cadom_i].RyRState_NP[1]; // O
	double dI = (k_23*segments[seg_i].caUnits[cadom_i].RyRState_NP[1]-k_32*segments[seg_i].caUnits[cadom_i].RyRState_NP[2])-(k_34*segments[seg_i].caUnits[cadom_i].RyRState_NP[2]-k_43*segments[seg_i].caUnits[cadom_i].RyRState_NP[3]);   // I
	double dCI = k_14*segments[seg_i].caUnits[cadom_i].RyRState_NP[0] - k_41*segments[seg_i].caUnits[cadom_i].RyRState_NP[3] + k_34*segments[seg_i].caUnits[cadom_i].RyRState_NP[2] - k_43*segments[seg_i].caUnits[cadom_i].RyRState_NP[3];
	
	double dR_B = (k_41*segments[seg_i].caUnits[cadom_i].RyRState_NP[7]-k_14*segments[seg_i].caUnits[cadom_i].RyRState_NP[4])-(k_12*segments[seg_i].caUnits[cadom_i].RyRState_NP[4]-k_21*segments[seg_i].caUnits[cadom_i].RyRState_NP[5]);   // R
	double dO_B = (k_12*segments[seg_i].caUnits[cadom_i].RyRState_NP[4]-k_21*segments[seg_i].caUnits[cadom_i].RyRState_NP[5])-(k_23*segments[seg_i].caUnits[cadom_i].RyRState_NP[5]-k_32*segments[seg_i].caUnits[cadom_i].RyRState_NP[6]) - k_62 * segments[seg_i].caUnits[cadom_i].RyRState_NP[5] + k_26 * segments[seg_i].caUnits[cadom_i].RyRState_NP[1]; // O
	double dI_B = (k_23*segments[seg_i].caUnits[cadom_i].RyRState_NP[5]-k_32*segments[seg_i].caUnits[cadom_i].RyRState_NP[6])-(k_34*segments[seg_i].caUnits[cadom_i].RyRState_NP[6]-k_43*segments[seg_i].caUnits[cadom_i].RyRState_NP[7]);   // I
	//double dCI = k_14*segments[seg_i].caUnits[cadom_i].RyRState_NP[0] - k_41*segments[seg_i].caUnits[cadom_i].RyRState_NP[3] + k_34*segments[seg_i].caUnits[cadom_i].RyRState_NP[2] - k_43*segments[seg_i].caUnits[cadom_i].RyRState_NP[3];
		
	segments[seg_i].caUnits[cadom_i].RyRState_NP[0] = check_State_Change(segments[seg_i].caUnits[cadom_i].RyRState_NP[0], dR, dt);	
	segments[seg_i].caUnits[cadom_i].RyRState_NP[1] = check_State_Change(segments[seg_i].caUnits[cadom_i].RyRState_NP[1], dO, dt);
	segments[seg_i].caUnits[cadom_i].RyRState_NP[2] = check_State_Change(segments[seg_i].caUnits[cadom_i].RyRState_NP[2], dI, dt);
	segments[seg_i].caUnits[cadom_i].RyRState_NP[3] = check_State_Change(segments[seg_i].caUnits[cadom_i].RyRState_NP[3], dCI, dt);
	segments[seg_i].caUnits[cadom_i].RyRState_NP[4] = check_State_Change(segments[seg_i].caUnits[cadom_i].RyRState_NP[4], dR_B, dt);
	segments[seg_i].caUnits[cadom_i].RyRState_NP[5] = check_State_Change(segments[seg_i].caUnits[cadom_i].RyRState_NP[5], dO_B, dt);
	segments[seg_i].caUnits[cadom_i].RyRState_NP[6] = check_State_Change(segments[seg_i].caUnits[cadom_i].RyRState_NP[6], dI_B, dt);
	
	if(segments[seg_i].caUnits[cadom_i].RyRState_NP[0] < 0) { segments[seg_i].caUnits[cadom_i].RyRState_NP[3] += segments[seg_i].caUnits[cadom_i].RyRState_NP[0]; segments[seg_i].caUnits[cadom_i].RyRState_NP[0] = 0; } // Compensated by State 3
	if(segments[seg_i].caUnits[cadom_i].RyRState_NP[1] < 0) { segments[seg_i].caUnits[cadom_i].RyRState_NP[2] += segments[seg_i].caUnits[cadom_i].RyRState_NP[1]; segments[seg_i].caUnits[cadom_i].RyRState_NP[1] = 0; }
	if(segments[seg_i].caUnits[cadom_i].RyRState_NP[5] < 0) { segments[seg_i].caUnits[cadom_i].RyRState_NP[6] += segments[seg_i].caUnits[cadom_i].RyRState_NP[5]; segments[seg_i].caUnits[cadom_i].RyRState_NP[5] = 0; }
	
	segments[seg_i].caUnits[cadom_i].RyRState_NP[7] = 1.0-segments[seg_i].caUnits[cadom_i].RyRState_NP[0]-segments[seg_i].caUnits[cadom_i].RyRState_NP[1]-segments[seg_i].caUnits[cadom_i].RyRState_NP[2]-segments[seg_i].caUnits[cadom_i].RyRState_NP[3]-segments[seg_i].caUnits[cadom_i].RyRState_NP[4]-segments[seg_i].caUnits[cadom_i].RyRState_NP[5]-segments[seg_i].caUnits[cadom_i].RyRState_NP[6];
	if(segments[seg_i].caUnits[cadom_i].RyRState_NP[7] < 0)
	{		
		segments[seg_i].caUnits[cadom_i].RyRState_NP[3] += segments[seg_i].caUnits[cadom_i].RyRState_NP[7];
		segments[seg_i].caUnits[cadom_i].RyRState_NP[7] = 0;
	}
	
	double RyR_SEP_SS = 1.0 / (1.0 + exp((Ca_cis_observed - 1.8E-3) / 1E-6)); // 1 / (1 + pow(Ca_cis_observed / 2.75E-4, 20));
	double tau_RyR_SEP = 60;
	segments[seg_i].caUnits[cadom_i].RyR_SEP = check_State_Change(segments[seg_i].caUnits[cadom_i].RyR_SEP, (RyR_SEP_SS - segments[seg_i].caUnits[cadom_i].RyR_SEP) / tau_RyR_SEP, dt);
	
	//double dRyR_Int = 1 * segments[seg_i].caUnits[cadom_i].Jrel - 0.01 * segments[seg_i].caUnits[cadom_i].RyR_Integr;
	double dRyR_Int = 1 * segments[seg_i].caUnits[cadom_i].Jrel - 0.1 / (1 + exp((segments[seg_i].caUnits[cadom_i].Jrel - 1.5E-3) / 1E-4)) * segments[seg_i].caUnits[cadom_i].RyR_Integr;
	segments[seg_i].caUnits[cadom_i].RyR_Integr = check_State_Change(segments[seg_i].caUnits[cadom_i].RyR_Integr, dRyR_Int, dt);
	
	// SERCA
	double tval_Kmr = pow(segments[seg_i].caUnits[cadom_i].Ca_sr / data->Kmr, data->hillSRCaP);
	double tval_Kmf_NP = pow(segments[seg_i].caUnits[cadom_i].Ca_i / data->Kmf_NP, data->hillSRCaP);
	double Jup_NP = (1 - settings->IupB) * data->Vmax_SRCaP * (tval_Kmf_NP - tval_Kmr) / (1 + tval_Kmf_NP + tval_Kmr);	
	//double Jup_P = (1 - settings->IupB) * data->Vmax_SRCaP * (tval_Kmf_P - tval_Kmr) / (1 + tval_Kmf_P + tval_Kmr);
	segments[seg_i].caUnits[cadom_i].Jup = Jup_NP;
	
}

double RAMM::Buff_Ca_DoubleBuff(double ca_t, double buff1_tot, double buff2_tot, double Kmbuff1, double Kmbuff2)
{
	
	double val;
	double b_buff = buff1_tot+buff2_tot-ca_t+Kmbuff1+Kmbuff2;
	double c_buff = Kmbuff1*Kmbuff2 -ca_t*(Kmbuff1+Kmbuff2)+buff2_tot*Kmbuff1+buff1_tot*Kmbuff2;
	double d_buff = -Kmbuff1*Kmbuff2*ca_t;      
    
	double sq_bbc_val = 2 * sqrt(b_buff*b_buff-3*c_buff);
	val = (sq_bbc_val/3.0)*cos(acos((9*b_buff*c_buff-2*b_buff*b_buff*b_buff-27*d_buff)/((b_buff*b_buff-3*c_buff)*sq_bbc_val))/3)-(b_buff/3);
	

	return val;
}

RAMM* RAMM::clone()
{
	RAMM* cl = new RAMM();
		cl->init(data, settings);
		cl->Vm = Vm;
		
		for(int ii=0; ii<settings->N_Segments; ii++)
		{	
			// Loop through the two membrane sides
			for(int jj=0; jj<2; jj++)
			{
				cl->segments[ii].membrane[jj].m_Na = segments[ii].membrane[jj].m_Na; cl->segments[ii].membrane[jj].h1_Na = segments[ii].membrane[jj].h1_Na; cl->segments[ii].membrane[jj].h2_Na = segments[ii].membrane[jj].h2_Na; 
				cl->segments[ii].membrane[jj].to_r = segments[ii].membrane[jj].to_r; cl->segments[ii].membrane[jj].to_s1 = segments[ii].membrane[jj].to_s1; cl->segments[ii].membrane[jj].to_s2 = segments[ii].membrane[jj].to_s2; cl->segments[ii].membrane[jj].to_s3 = segments[ii].membrane[jj].to_s3;
				cl->segments[ii].membrane[jj].pa_Kr = segments[ii].membrane[jj].pa_Kr; cl->segments[ii].membrane[jj].pi_Kr = segments[ii].membrane[jj].pi_Kr; 
				cl->segments[ii].membrane[jj].z_Ks = segments[ii].membrane[jj].z_Ks; 
				cl->segments[ii].membrane[jj].d_L = segments[ii].membrane[jj].d_L; cl->segments[ii].membrane[jj].f_L = segments[ii].membrane[jj].f_L; 
				cl->segments[ii].membrane[jj].d_T = segments[ii].membrane[jj].d_T; cl->segments[ii].membrane[jj].f_T = segments[ii].membrane[jj].f_T; 
		
				//for(int abcd = 0; abcd<8; abcd++)
				//{
				//	cl->segments[ii].membrane[jj].ICaLState_NP[abcd] = segments[ii].membrane[jj].ICaLState_NP[abcd];
				//}
		
				cl->segments[ii].membrane[jj].Buff_Na_junc = segments[ii].membrane[jj].Buff_Na_junc; cl->segments[ii].membrane[jj].Buff_Na_sl = segments[ii].membrane[jj].Buff_Na_sl;
				cl->segments[ii].membrane[jj].Buff_Ca_SLLj = segments[ii].membrane[jj].Buff_Ca_SLLj; cl->segments[ii].membrane[jj].Buff_Ca_SLLsl = segments[ii].membrane[jj].Buff_Ca_SLLsl;
				cl->segments[ii].membrane[jj].Buff_Ca_SLHj = segments[ii].membrane[jj].Buff_Ca_SLHj; cl->segments[ii].membrane[jj].Buff_Ca_SLHsl = segments[ii].membrane[jj].Buff_Ca_SLHsl;
				cl->segments[ii].membrane[jj].Buff_Ca_EGTA_sl = segments[ii].membrane[jj].Buff_Ca_EGTA_sl;
								
				cl->segments[ii].membrane[jj].Ca_junc = segments[ii].membrane[jj].Ca_junc;
				cl->segments[ii].membrane[jj].Ca_sl = segments[ii].membrane[jj].Ca_sl;
				cl->segments[ii].membrane[jj].Na_junc = segments[ii].membrane[jj].Na_junc;
				cl->segments[ii].membrane[jj].Na_sl = segments[ii].membrane[jj].Na_sl;
			}
			
			// Loop through the cytosolic Ca2+ domains
			for(int jj=0; jj<settings->N_CaDomains; jj++)
			{
				cl->segments[ii].caUnits[jj].RyRState_NP[0] = segments[ii].caUnits[jj].RyRState_NP[0]; cl->segments[ii].caUnits[jj].RyRState_NP[1] = segments[ii].caUnits[jj].RyRState_NP[1]; cl->segments[ii].caUnits[jj].RyRState_NP[2] = segments[ii].caUnits[jj].RyRState_NP[2]; 
				cl->segments[ii].caUnits[jj].RyRState_NP[3] = segments[ii].caUnits[jj].RyRState_NP[3]; cl->segments[ii].caUnits[jj].RyRState_NP[4] = segments[ii].caUnits[jj].RyRState_NP[4]; cl->segments[ii].caUnits[jj].RyRState_NP[5] = segments[ii].caUnits[jj].RyRState_NP[5]; 
				cl->segments[ii].caUnits[jj].RyRState_NP[6] = segments[ii].caUnits[jj].RyRState_NP[6]; cl->segments[ii].caUnits[jj].RyRState_NP[7] = segments[ii].caUnits[jj].RyRState_NP[7];
				
				cl->segments[ii].caUnits[jj].RyR_SEP = segments[ii].caUnits[jj].RyR_SEP;
				cl->segments[ii].caUnits[jj].RyR_Integr = segments[ii].caUnits[jj].RyR_Integr;
								
				cl->segments[ii].caUnits[jj].Buff_Ca_TnCL = segments[ii].caUnits[jj].Buff_Ca_TnCL;
				cl->segments[ii].caUnits[jj].Buff_Ca_TnCHc = segments[ii].caUnits[jj].Buff_Ca_TnCHc;
				cl->segments[ii].caUnits[jj].Buff_Ca_TnCHm = segments[ii].caUnits[jj].Buff_Ca_TnCHm;
				cl->segments[ii].caUnits[jj].Buff_Ca_CaM = segments[ii].caUnits[jj].Buff_Ca_CaM;
				cl->segments[ii].caUnits[jj].Buff_Ca_Myosin_ca = segments[ii].caUnits[jj].Buff_Ca_Myosin_ca;
				cl->segments[ii].caUnits[jj].Buff_Ca_Myosin_mg = segments[ii].caUnits[jj].Buff_Ca_Myosin_mg;
				cl->segments[ii].caUnits[jj].Buff_Ca_SRB = segments[ii].caUnits[jj].Buff_Ca_SRB;
				cl->segments[ii].caUnits[jj].Buff_Ca_CSQN = segments[ii].caUnits[jj].Buff_Ca_CSQN;
				cl->segments[ii].caUnits[jj].Buff_Ca_EGTA_cyt = segments[ii].caUnits[jj].Buff_Ca_EGTA_cyt;
		
				cl->segments[ii].caUnits[jj].Ca_i = segments[ii].caUnits[jj].Ca_i;
				cl->segments[ii].caUnits[jj].Ca_sr = segments[ii].caUnits[jj].Ca_sr;
				
				cout << "Clone RyR_State = [" << cl->segments[ii].caUnits[jj].RyRState_NP[0] << "," << cl->segments[ii].caUnits[jj].RyRState_NP[1] << "," << cl->segments[ii].caUnits[jj].RyRState_NP[2] << "," << cl->segments[ii].caUnits[jj].RyRState_NP[3] << "," << cl->segments[ii].caUnits[jj].RyRState_NP[4] << "," << cl->segments[ii].caUnits[jj].RyRState_NP[5] << "," << cl->segments[ii].caUnits[jj].RyRState_NP[6] << "," << cl->segments[ii].caUnits[jj].RyRState_NP[7] << "]" << endl;
			}
			
			// Other stuff
			cl->segments[ii].Na_i = segments[ii].Na_i;
			cl->segments[ii].K_i = segments[ii].K_i;
		}
		
	return cl;
}

void RAMM::outputData(FILE* outputfile, int levelOfElectroDetail, int levelOfSignalingDetail, bool binaryOutput)
{     
	// Minimum Level of Detail:
	//   Vm 
	if(levelOfElectroDetail == 1)
	{
		fprintf(outputfile, "%-8e\t", Vm);		
	}

	// Average Level of Detail (1):
	if(levelOfElectroDetail == 2)
	{		
		double Irel_adjust = Jrel + Jleak;
		fprintf(outputfile,"%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t",Vm,Ca_i,Ca_sr,(ICaL_junc+ICaL_sl),Irel_adjust,INaCa_junc+INaCa_sl,Jup,Ca_sl,ITo,Na_i,Ca_junc,INaK_junc+INaK_sl, IKr, IK1, K_i);
	}

	if(levelOfElectroDetail == 21)
	{
		double kiont_alt = IKr+IKs+IK1+Isus-2*(INaK_junc+INaK_sl)+ITo;
		if(!settings->applyVoltageClamp) kiont_alt +=Istim;	
		double dVdt = caiont + naiont + kiont_alt + clont;
		double Avg_Ca_sr_tot = 0;
		double totnumdoms = (1.0 * settings->N_Segments * settings->N_CaDomains);
		for(int i=0; i<settings->N_Segments; i++) 
		{
			for(int j=0; j<settings->N_CaDomains; j++)
			{
				Avg_Ca_sr_tot += (1.0 / totnumdoms) * (segments[i].caUnits[j].Ca_sr + segments[i].caUnits[j].Buff_Ca_CSQN);
			}
		}
		double Irel_adjust = Jrel + Jleak;
		
		fprintf(outputfile,"%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t%-8e\t",Vm,1000*Ca_i,Ca_sr,dVdt,Irel_adjust,ICab_junc+ICab_sl,INaCa_junc+INaCa_sl,Jup,Ca_sl,Na_i,Ca_junc,INaK_junc+INaK_sl, Avg_Ca_sr_tot, ICaL_junc+ICaL_sl, INaCa_junc, K_i);
	}
	
	if(levelOfElectroDetail == 3)
	{
		fprintf(outputfile, "%-4e\t%-4e\t%-4e\t", Vm, Ca_i, Na_i);
	}
	
	if(levelOfElectroDetail == 41)
	{
		fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t", Vm, 1000*Ca_i,(ICaL_junc+ICaL_sl), (ICaT_junc+ICaT_sl));
	}
	
	if(levelOfElectroDetail == 42)
	{
		fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", Vm, Ca_i, Ca_sr, (ICaL_junc+ICaL_sl), (ICaT_junc+ICaT_sl), (INaK_junc+INaK_sl), (INa_junc+INa_sl), IK1, ITo, (INaCa_junc+INaCa_sl), (INab_junc+INab_sl),IKr, IKs, Isus, IClb, (IClCa_junc+IClCa_sl), Na_i, K_i, Jup, Jleak, Jrel, Istim);
	}
	
	if(levelOfElectroDetail == 43)
	{
                double Irel_adjust = Jleak + Jrel;

                fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", Vm, Ca_i, Ca_sr, (ICaL_junc+ICaL_sl), (ICaT_junc+ICaT_sl), (INaK_junc+INaK_sl), (INa_junc+INa_sl), IK1, ITo, (INaCa_junc+INaCa_sl), IKr, IKs, Isus, (INab_junc+INab_sl), Na_i, K_i, Jup, Istim);

                for(int i=0; i<settings->N_Segments; i+=settings->segOutputRes) 
                {
                        for(int j=0; j<settings->N_CaDomains; j+=settings->domOutputRes)
                        {
                                fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t", segments[i].caUnits[j].Ca_i, segments[i].caUnits[j].Ca_sr, segments[i].caUnits[j].Jrel, segments[i].caUnits[j].Jup);
                        }
                }
	}
	
	if(levelOfElectroDetail == 44)
	{

		double Irel_adjust = Jleak + Jrel;
                double m3 = segments[0].membrane[0].m_Na*segments[0].membrane[0].m_Na*segments[0].membrane[0].m_Na;
                double h_tot = 0.635*segments[0].membrane[0].h1_Na + 0.365*segments[0].membrane[0].h2_Na;

                double S1 = segments[0].membrane[0].to_s1;
                double S2 = segments[0].membrane[0].to_s2;
                double S3 = segments[0].membrane[0].to_s3;

                fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", Vm, Ca_i, Ca_sr, (ICaL_junc+ICaL_sl), (ICaT_junc+ICaT_sl), (INaK_junc+INaK_sl), (INa_junc+INa_sl), IK1, ITo, (INaCa_junc+INaCa_sl), IKr, IKs, Isus, (INab_junc+INab_sl), IClb, (IClCa_junc+IClCa_sl), (ISK_junc+ISK_sl),IKAch);
                fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t", m3, h_tot, tauh1, tauh2);
                fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", segments[0].membrane[0].d_L, segments[0].membrane[0].f_L, tau_d_L, tau_f_L, d_L_infinity, f_L_infinity);
                fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", segments[0].membrane[0].to_r,S1, S2,S3, taur, taus1, taus2, taus3);
                fprintf(outputfile, "%-4e\t%-4e\t%-4e\t", segments[0].membrane[0].z_Ks, segments[0].membrane[0].pa_Kr, segments[0].membrane[0].pi_Kr);
                fprintf(outputfile, "%-4e\t%-4e\t", Na_i, K_i);
		fprintf(outputfile, "%-4e\t", Istim);

		for(int i=0; i<settings->N_Segments; i+=settings->segOutputRes) 
                {
                	for(int j=0; j<settings->N_CaDomains; j+=settings->domOutputRes)
                        {
                        	//fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", segments[i].caUnits[j].Ca_i, segments[i].caUnits[j].Ca_sr, segments[i].caUnits[j].Jrel, segments[i].caUnits[j].Ca_srs, segments[i].caUnits[j].Jup, segments[i].caUnits[j].J_diff_i_subcell);
                                fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t", segments[i].caUnits[j].Ca_i, segments[i].caUnits[j].Ca_sr, segments[i].caUnits[j].Jrel, segments[i].caUnits[j].Ca_srs);
                        }
                }
	}

	if(levelOfElectroDetail == 45)
	{
		double Irel_adjust = Jleak + Jrel;
		fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", Vm, Ca_i, Ca_sr, Ca_sl, Ca_junc, Irel_adjust, (ICaL_junc+ICaL_sl));
		for(int i=0; i<settings->N_Segments; i++) 
		{
			for(int j=0; j<settings->N_CaDomains; j++)
			{
				fprintf(outputfile, "%-4e\t%-4e\t%-4e\t", segments[i].caUnits[j].Ca_i, segments[i].caUnits[j].Ca_sr, segments[i].caUnits[j].Jrel);
			}
		}
	}
	
	if(levelOfElectroDetail == 46)
	{
		fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", Vm, 1000*Ca_i, INa_junc + INa_sl, INaK_junc + INaK_sl, INaCa_junc+INaCa_sl, ICaL_junc+ICaL_sl, ITo, IK1, IKr, IKs, Isus);
	}
	
	if(levelOfElectroDetail == 47)
	{
		fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", Vm, 1000*Ca_i, INa_junc + INa_sl, segments[0].membrane[0].m_Na, segments[0].membrane[0].h1_Na, segments[0].membrane[0].h2_Na);
	}
	
	if(levelOfElectroDetail == 48)
	{
		double Irel_adjust = Jleak + Jrel;
		fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", Vm, Ca_i, Ca_sr, Ca_sl, Ca_junc, Irel_adjust, Jup, INaCa_junc+INaCa_sl, ICaL_junc+ICaL_sl);
		for(int i=0; i<settings->N_Segments; i++) 
		{
			for(int j=0; j<settings->N_CaDomains; j++)
			{
				fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", segments[i].caUnits[j].Ca_i, segments[i].caUnits[j].Ca_sr, segments[i].caUnits[j].Jrel, segments[i].caUnits[j].RyRState_NP[0], segments[i].caUnits[j].RyRState_NP[1], segments[i].caUnits[j].RyRState_NP[2], segments[i].caUnits[j].Jup, segments[i].caUnits[j].Ca_srs);
			}
		}
	}
	
	if(levelOfElectroDetail == 49)
	{
		double Irel_adjust = Jleak + Jrel;
		fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", Vm, Ca_i, Ca_sr, Ca_sl, Ca_junc, Irel_adjust, Jup);		
	}
	
	if(levelOfElectroDetail == 5)
	{
		double kiont_alt = IKr+IKs+IK1+Isus-2*(INaK_junc+INaK_sl)+ITo;
		if(!settings->applyVoltageClamp) kiont_alt +=Istim;
		double Irel_adjust = Jleak + Jrel;
		double ICaL_tot = ICaL_junc + ICaL_sl;
		
		if(!binaryOutput)
		{
			fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", Vm, 1000*Ca_i, Ca_sr, Ca_sl, Buff_Ca_EGTA_cyt, Irel_adjust, ICaL_tot, caiont + naiont + kiont_alt + clont);
			for(int i=0; i<settings->N_Segments; i+=settings->segOutputRes) 
			{
				for(int j=0; j<settings->N_CaDomains; j+=settings->domOutputRes)
				{
					//fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t", segments[i].caUnits[j].Ca_i, segments[i].caUnits[j].Ca_sr, segments[i].caUnits[j].Jrel, segments[i].caUnits[j].Ca_srs);
					fprintf(outputfile, "%-4e\t%-4e\t%-4e\t", segments[i].caUnits[j].Ca_i, segments[i].caUnits[j].Ca_sr, segments[i].caUnits[j].Jrel);
				}
			}
		}
		else
		{
			size_t doubles = sizeof(double);
			double Cai_x_1000 = 1000 * Ca_i;
			double dVdt = caiont + naiont + kiont_alt + clont;
			fwrite(&(Vm), doubles, 1, outputfile); fwrite(&(Cai_x_1000), doubles, 1, outputfile); fwrite(&(Ca_sr), doubles, 1, outputfile); fwrite(&(Ca_sl), doubles, 1, outputfile); 
			fwrite(&(Buff_Ca_EGTA_cyt), doubles, 1, outputfile); fwrite(&(Irel_adjust), doubles, 1, outputfile); fwrite(&(ICaL_tot), doubles, 1, outputfile); fwrite(&(dVdt), doubles, 1, outputfile); 
			for(int i=0; i<settings->N_Segments; i+=settings->segOutputRes) 
			{
				for(int j=0; j<settings->N_CaDomains; j+=settings->domOutputRes)
				{									
					fwrite(&(segments[i].caUnits[j].Ca_i), doubles, 1, outputfile); //fwrite(&(segments[i].caUnits[j].Ca_sr), doubles, 1, outputfile); fwrite(&(segments[i].caUnits[j].Jrel), doubles, 1, outputfile); //fwrite(&(segments[i].caUnits[j].Ca_srs), doubles, 1, outputfile); 
				}
			}
		}
	}
	
	if(levelOfElectroDetail == 51)
	{
		double kiont_alt = IKr+IKs+IK1+Isus-2*(INaK_junc+INaK_sl)+ITo;
		if(!settings->applyVoltageClamp) kiont_alt +=Istim;
		
		double TotCellCa = 0;
		double totnumdoms = (1.0 * settings->N_Segments * settings->N_CaDomains);
		for(int i=0; i<settings->N_Segments; i++) 
		{
			//TotCellCa += (data->vsl / settings->N_Segments) / data->vcell * (segments[i].membrane[0].Ca_sl + segments[i].membrane[0].Buff_Ca_SLLsl + segments[i].membrane[0].Buff_Ca_SLHsl + segments[i].membrane[1].Ca_sl + segments[i].membrane[1].Buff_Ca_SLLsl + segments[i].membrane[1].Buff_Ca_SLHsl);
			//TotCellCa += (data->vjunc / settings->N_Segments) / data->vcell * (segments[i].membrane[0].Ca_junc + segments[i].membrane[0].Buff_Ca_SLLj + segments[i].membrane[0].Buff_Ca_SLHj + segments[i].membrane[1].Ca_junc + segments[i].membrane[1].Buff_Ca_SLLj + segments[i].membrane[1].Buff_Ca_SLHj);
			
			for(int j=0; j<settings->N_CaDomains; j++)
			{
				TotCellCa += (data->vmyo / totnumdoms) / data->vcell * (segments[i].caUnits[j].Ca_i + segments[i].caUnits[j].Buff_Ca_TnCL + segments[i].caUnits[j].Buff_Ca_TnCHc + segments[i].caUnits[j].Buff_Ca_TnCHm + segments[i].caUnits[j].Buff_Ca_CaM + segments[i].caUnits[j].Buff_Ca_Myosin_ca + segments[i].caUnits[j].Buff_Ca_Myosin_mg + segments[i].caUnits[j].Buff_Ca_SRB);
				//TotCellCa += (data->vsr / totnumdoms) / data->vcell * (segments[i].caUnits[j].Ca_sr + segments[i].caUnits[j].Buff_Ca_CSQN);
				//TotCellCa += (data->vjunc / totnumdoms) / data->vcell * (segments[i].caUnits[j].Ca_srs + segments[i].caUnits[j].Buff_Ca_SLLsrs + segments[i].caUnits[j].Buff_Ca_SLHsrs);
			}
		}
		double Irel_adjust = Jleak + Jrel;
		
		if(!binaryOutput)
		{	
			fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", Vm, 1000*Ca_i, Ca_sr, Ca_sl, Buff_Ca_EGTA_cyt, Irel_adjust, (ICaL_junc+ICaL_sl), caiont + naiont + kiont_alt + clont);
			//fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", segments[0].caUnits[0].Buff_Ca_TnCL, segments[0].caUnits[0].Buff_Ca_TnCHc, segments[0].caUnits[0].Buff_Ca_TnCHm, segments[0].caUnits[0].Buff_Ca_CaM, segments[0].caUnits[0].Buff_Ca_Myosin_ca, segments[0].caUnits[0].Buff_Ca_Myosin_mg, segments[0].caUnits[0].Buff_Ca_SRB);
			//fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", segments[0].membrane[0].ICaLState_NP[0], segments[0].membrane[0].ICaLState_NP[1], segments[0].membrane[0].ICaLState_NP[2], segments[0].membrane[0].ICaLState_NP[3], segments[0].membrane[0].ICaLState_NP[4], segments[0].membrane[0].ICaLState_NP[5], segments[0].membrane[0].ICaLState_NP[6], segments[0].membrane[0].Ca_sl,segments[0].caUnits[0].Ca_srs);
			for(int i=0; i<settings->N_Segments; i+=settings->segOutputRes) 
			{
				for(int j=0; j<settings->N_CaDomains; j+=settings->domOutputRes)
				{
					//fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", segments[i].caUnits[j].Ca_i, segments[i].caUnits[j].Ca_sr, segments[i].caUnits[j].Jrel, segments[i].caUnits[j].Ca_srs, segments[i].caUnits[j].Jup, segments[i].caUnits[j].J_diff_i_subcell);
					fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t", segments[i].caUnits[j].Ca_i, segments[i].caUnits[j].Ca_sr, segments[i].caUnits[j].Jrel, segments[i].caUnits[j].Ca_srs);
				}
			}
		}
		else
		{
			size_t doubles = sizeof(double);
			double Cai_x_1000 = 1000 * Ca_i;
			double ICaL_tot = (ICaL_junc+ICaL_sl);
			double dVdt = caiont + naiont + kiont_alt + clont;
			fwrite(&(Vm), doubles, 1, outputfile); fwrite(&(Cai_x_1000), doubles, 1, outputfile); fwrite(&(Ca_sr), doubles, 1, outputfile); fwrite(&(Ca_sl), doubles, 1, outputfile); 
			fwrite(&(Buff_Ca_EGTA_cyt), doubles, 1, outputfile); fwrite(&(Irel_adjust), doubles, 1, outputfile); fwrite(&(ICaL_tot), doubles, 1, outputfile); fwrite(&(dVdt), doubles, 1, outputfile); 
			//fwrite(segments[0].membrane[0].ICaLState_NP, doubles, 8, outputfile); 
			fwrite(segments[0].caUnits[0].RyRState_NP, doubles, 4, outputfile);
			fwrite(segments[0].caUnits[1].RyRState_NP, doubles, 4, outputfile); fwrite(&(segments[0].caUnits[0].Ca_srs), doubles, 1, outputfile);
			for(int i=0; i<settings->N_Segments; i+=settings->segOutputRes) 
			{
				for(int j=0; j<settings->N_CaDomains; j+=settings->domOutputRes)
				{									
					fwrite(&(segments[i].caUnits[j].Ca_i), doubles, 1, outputfile); fwrite(&(segments[i].caUnits[j].Ca_sr), doubles, 1, outputfile); fwrite(&(segments[i].caUnits[j].Jrel), doubles, 1, outputfile); fwrite(&(segments[i].caUnits[j].Ca_srs), doubles, 1, outputfile); 
				}
			}
		}
	}
	
	if(levelOfElectroDetail == 52)
	{
		double kiont_alt = IKr+IKs+IK1+Isus-2*(INaK_junc+INaK_sl)+ITo;
		if(!settings->applyVoltageClamp) kiont_alt +=Istim;
		
		double TotCellCa = 0;
		double totnumdoms = (1.0 * settings->N_Segments * settings->N_CaDomains);
		for(int i=0; i<settings->N_Segments; i++) 
		{
			//TotCellCa += (data->vsl / settings->N_Segments) / data->vcell * (segments[i].membrane[0].Ca_sl + segments[i].membrane[0].Buff_Ca_SLLsl + segments[i].membrane[0].Buff_Ca_SLHsl + segments[i].membrane[1].Ca_sl + segments[i].membrane[1].Buff_Ca_SLLsl + segments[i].membrane[1].Buff_Ca_SLHsl);
			//TotCellCa += (data->vjunc / settings->N_Segments) / data->vcell * (segments[i].membrane[0].Ca_junc + segments[i].membrane[0].Buff_Ca_SLLj + segments[i].membrane[0].Buff_Ca_SLHj + segments[i].membrane[1].Ca_junc + segments[i].membrane[1].Buff_Ca_SLLj + segments[i].membrane[1].Buff_Ca_SLHj);
			
			for(int j=0; j<settings->N_CaDomains; j++)
			{
				TotCellCa += (data->vmyo / totnumdoms) / data->vcell * (segments[i].caUnits[j].Ca_i + segments[i].caUnits[j].Buff_Ca_TnCL + segments[i].caUnits[j].Buff_Ca_TnCHc + segments[i].caUnits[j].Buff_Ca_TnCHm + segments[i].caUnits[j].Buff_Ca_CaM + segments[i].caUnits[j].Buff_Ca_Myosin_ca + segments[i].caUnits[j].Buff_Ca_Myosin_mg + segments[i].caUnits[j].Buff_Ca_SRB);
				//TotCellCa += (data->vsr / totnumdoms) / data->vcell * (segments[i].caUnits[j].Ca_sr + segments[i].caUnits[j].Buff_Ca_CSQN);
				//TotCellCa += (data->vjunc / totnumdoms) / data->vcell * (segments[i].caUnits[j].Ca_srs + segments[i].caUnits[j].Buff_Ca_SLLsrs + segments[i].caUnits[j].Buff_Ca_SLHsrs);
			}
		}
		double Irel_adjust = Jleak + Jrel;
		
		if(!binaryOutput)
		{	
			fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", Vm, 1000*Ca_i, Ca_sr, Ca_sl, Buff_Ca_EGTA_cyt, Irel_adjust, (ICaL_junc+ICaL_sl), caiont + naiont + kiont_alt + clont);
			//fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", segments[0].membrane[0].ICaLState_NP[0], segments[0].membrane[0].ICaLState_NP[1], segments[0].membrane[0].ICaLState_NP[2], segments[0].membrane[0].ICaLState_NP[3], segments[0].membrane[0].ICaLState_NP[4], segments[0].membrane[0].ICaLState_NP[5], segments[0].membrane[0].ICaLState_NP[6], segments[0].membrane[0].Ca_sl,segments[0].caUnits[0].Ca_srs);
			for(int i=0; i<settings->N_Segments; i+=settings->segOutputRes) 
			{
				for(int j=0; j<settings->N_CaDomains; j+=settings->domOutputRes)
				{
					fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", segments[i].caUnits[j].Ca_i, segments[i].caUnits[j].Ca_sr, segments[i].caUnits[j].Jrel, segments[i].caUnits[j].Ca_srs, segments[i].caUnits[j].RyRState_NP[0], segments[i].caUnits[j].RyRState_NP[1], segments[i].caUnits[j].RyRState_NP[2], segments[i].caUnits[j].RyRState_NP[3]);
				}
			}
		}
		else
		{
			size_t doubles = sizeof(double);
			double Cai_x_1000 = 1000 * Ca_i;
			double ICaL_tot = (ICaL_junc+ICaL_sl);
			double dVdt = caiont + naiont + kiont_alt + clont;
			fwrite(&(Vm), doubles, 1, outputfile); fwrite(&(Cai_x_1000), doubles, 1, outputfile); fwrite(&(Ca_sr), doubles, 1, outputfile); fwrite(&(Ca_sl), doubles, 1, outputfile); 
			fwrite(&(Buff_Ca_EGTA_cyt), doubles, 1, outputfile); fwrite(&(Irel_adjust), doubles, 1, outputfile); fwrite(&(ICaL_tot), doubles, 1, outputfile); fwrite(&(dVdt), doubles, 1, outputfile); 
			//fwrite(segments[0].membrane[0].ICaLState_NP, doubles, 8, outputfile); 
			fwrite(segments[0].caUnits[0].RyRState_NP, doubles, 4, outputfile);
			fwrite(segments[0].caUnits[1].RyRState_NP, doubles, 4, outputfile); fwrite(&(segments[0].caUnits[0].Ca_srs), doubles, 1, outputfile);
			for(int i=0; i<settings->N_Segments; i+=settings->segOutputRes) 
			{
				for(int j=0; j<settings->N_CaDomains; j+=settings->domOutputRes)
				{									
					fwrite(&(segments[i].caUnits[j].Ca_i), doubles, 1, outputfile); fwrite(&(segments[i].caUnits[j].Ca_sr), doubles, 1, outputfile); fwrite(&(segments[i].caUnits[j].Jrel), doubles, 1, outputfile); fwrite(&(segments[i].caUnits[j].Ca_srs), doubles, 1, outputfile); 
					fwrite(&(segments[i].caUnits[j].RyRState_NP), doubles, 4, outputfile);
				}
			}
		}
	}
	
	if(levelOfElectroDetail == 53)
	{
		double kiont_alt = IKr+IKs+IK1+Isus-2*(INaK_junc+INaK_sl)+ITo;
		if(!settings->applyVoltageClamp) kiont_alt +=Istim;
		
		double TotCellCa = 0;
		double totnumdoms = (1.0 * settings->N_Segments * settings->N_CaDomains);
		for(int i=0; i<settings->N_Segments; i++) 
		{
			for(int j=0; j<settings->N_CaDomains; j++)
			{
				TotCellCa += (data->vmyo / totnumdoms) / data->vcell * (segments[i].caUnits[j].Ca_i + segments[i].caUnits[j].Buff_Ca_TnCL + segments[i].caUnits[j].Buff_Ca_TnCHc + segments[i].caUnits[j].Buff_Ca_TnCHm + segments[i].caUnits[j].Buff_Ca_CaM + segments[i].caUnits[j].Buff_Ca_Myosin_ca + segments[i].caUnits[j].Buff_Ca_Myosin_mg + segments[i].caUnits[j].Buff_Ca_SRB);
			}
		}
		double Irel_adjust = Jleak + Jrel;
		
		if(!binaryOutput)
		{	
			fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", Vm, 1000*Ca_i, Ca_sr, Ca_sl, Buff_Ca_EGTA_cyt, Irel_adjust, (ICaL_junc+ICaL_sl), caiont + naiont + kiont_alt + clont);
			//fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", segments[0].membrane[0].ICaLState_NP[0], segments[0].membrane[0].ICaLState_NP[1], segments[0].membrane[0].ICaLState_NP[2], segments[0].membrane[0].ICaLState_NP[3], segments[0].membrane[0].ICaLState_NP[4], segments[0].membrane[0].ICaLState_NP[5], segments[0].membrane[0].ICaLState_NP[6], segments[0].membrane[0].Ca_sl,segments[0].caUnits[0].Ca_srs);
			for(int i=0; i<settings->N_Segments; i+=settings->segOutputRes) 
			{
				for(int j=0; j<settings->N_CaDomains; j+=settings->domOutputRes)
				{
					fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", segments[i].caUnits[j].Ca_i, segments[i].caUnits[j].Ca_sr, segments[i].caUnits[j].Jrel, segments[i].caUnits[j].Ca_srs, segments[i].caUnits[j].RyR_Integr);
				}
			}
		}
		else
		{
			size_t doubles = sizeof(double);
			double Cai_x_1000 = 1000 * Ca_i;
			double ICaL_tot = (ICaL_junc+ICaL_sl);
			double dVdt = caiont + naiont + kiont_alt + clont;
			fwrite(&(Vm), doubles, 1, outputfile); fwrite(&(Cai_x_1000), doubles, 1, outputfile); fwrite(&(Ca_sr), doubles, 1, outputfile); fwrite(&(Ca_sl), doubles, 1, outputfile); 
			fwrite(&(Buff_Ca_EGTA_cyt), doubles, 1, outputfile); fwrite(&(Irel_adjust), doubles, 1, outputfile); fwrite(&(ICaL_tot), doubles, 1, outputfile); fwrite(&(dVdt), doubles, 1, outputfile); 
			fwrite(segments[0].caUnits[0].RyRState_NP, doubles, 4, outputfile);
			fwrite(segments[0].caUnits[1].RyRState_NP, doubles, 4, outputfile); fwrite(&(segments[0].caUnits[0].Ca_srs), doubles, 1, outputfile);
			for(int i=0; i<settings->N_Segments; i+=settings->segOutputRes) 
			{
				for(int j=0; j<settings->N_CaDomains; j+=settings->domOutputRes)
				{									
					fwrite(&(segments[i].caUnits[j].Ca_i), doubles, 1, outputfile); 
					fwrite(&(segments[i].caUnits[j].Ca_sr), doubles, 1, outputfile); 
					fwrite(&(segments[i].caUnits[j].Jrel), doubles, 1, outputfile); 
					fwrite(&(segments[i].caUnits[j].Ca_srs), doubles, 1, outputfile); 
					fwrite(&(segments[i].caUnits[j].RyR_Integr), doubles, 1, outputfile); 
				}
			}
		}
	}
	
	if(levelOfElectroDetail == 9)
	{		
		double Irel_adjust = Jleak + Jrel;
		fprintf(outputfile, "%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t%-4e\t", Vm, 1000*Ca_i, Ca_sr, Ca_sl, Ca_junc, Irel_adjust, INa_junc+INa_sl, (ICaL_junc+ICaL_sl), IKr, IKs, ITo, Isus, IK1, INaK_junc+INaK_sl, INaCa_junc+INaCa_sl, Na_i);
	}

	// Full detail
	//    Also used to store X0, so don't alter this!
	if(levelOfElectroDetail == 10)
	{
		fprintf(outputfile, "%-8e\t", Vm);
		
		for(int ii=0; ii<settings->N_Segments; ii++)
		{	
			// Loop through the two membrane sides
			for(int jj=0; jj<2; jj++)
			{
				fprintf(outputfile, "%-8e\t%-8e\t%-8e\t", segments[ii].membrane[jj].m_Na, segments[ii].membrane[jj].h1_Na, segments[ii].membrane[jj].h2_Na); 
				fprintf(outputfile, "%-8e\t%-8e\t%-8e\t%-8e\t", segments[ii].membrane[jj].to_r, segments[ii].membrane[jj].to_s1, segments[ii].membrane[jj].to_s2, segments[ii].membrane[jj].to_s3); 
				fprintf(outputfile, "%-8e\t%-8e\t", segments[ii].membrane[jj].pa_Kr, segments[ii].membrane[jj].pi_Kr);
				fprintf(outputfile, "%-8e\t", segments[ii].membrane[jj].z_Ks);
				fprintf(outputfile, "%-8e\t,%-8e\t", segments[ii].membrane[jj].d_L, segments[ii].membrane[jj].f_L);
				fprintf(outputfile, "%-8e\t%-8e\t", segments[ii].membrane[jj].d_T, segments[ii].membrane[jj].f_T);
		
				fprintf(outputfile, "%-8e\t%-8e\t", segments[ii].membrane[jj].Buff_Na_junc, segments[ii].membrane[jj].Buff_Na_sl);
				fprintf(outputfile, "%-8e\t%-8e\t", segments[ii].membrane[jj].Buff_Ca_SLLj, segments[ii].membrane[jj].Buff_Ca_SLLsl);
				fprintf(outputfile, "%-8e\t%-8e\t", segments[ii].membrane[jj].Buff_Ca_SLHj, segments[ii].membrane[jj].Buff_Ca_SLHsl);
				fprintf(outputfile, "%-8e\t", segments[ii].membrane[jj].Buff_Ca_EGTA_sl);
				
				fprintf(outputfile, "%-8e\t%-8e\t", segments[ii].membrane[jj].Ca_junc, segments[ii].membrane[jj].Ca_sl);
				fprintf(outputfile, "%-8e\t%-8e\t", segments[ii].membrane[jj].Na_junc, segments[ii].membrane[jj].Na_sl);
				
				//for(int abcd=0; abcd<8; abcd++)
				//{
				//	fprintf(outputfile, "%-8e\t", segments[ii].membrane[jj].ICaLState_NP[abcd]);
				//}
				
				//for(int abcd=0; abcd<4; abcd++)
				//{
				//	fprintf(outputfile, "%-8e\t", segments[ii].membrane[jj].SK_State[abcd]);
				//}
			}
			
			// Loop through the cytosolic Ca2+ domains
			for(int jj=0; jj<settings->N_CaDomains; jj++)
			{
				fprintf(outputfile, "%-8e\t%-8e\t%-8e\t%-8e\t",segments[ii].caUnits[jj].RyRState_NP[0], segments[ii].caUnits[jj].RyRState_NP[1], segments[ii].caUnits[jj].RyRState_NP[2], segments[ii].caUnits[jj].RyRState_NP[3]);
				fprintf(outputfile, "%-8e\t%-8e\t%-8e\t%-8e\t",segments[ii].caUnits[jj].RyRState_NP[4], segments[ii].caUnits[jj].RyRState_NP[5], segments[ii].caUnits[jj].RyRState_NP[6], segments[ii].caUnits[jj].RyRState_NP[7]);
				
				fprintf(outputfile, "%-8e\t%-8e\t%-8e\t", segments[ii].caUnits[jj].Buff_Ca_TnCL, segments[ii].caUnits[jj].Buff_Ca_TnCHc, segments[ii].caUnits[jj].Buff_Ca_TnCHm);
				fprintf(outputfile, "%-8e\t%-8e\t%-8e\t", segments[ii].caUnits[jj].Buff_Ca_CaM, segments[ii].caUnits[jj].Buff_Ca_Myosin_ca, segments[ii].caUnits[jj].Buff_Ca_Myosin_mg);
				fprintf(outputfile, "%-8e\t%-8e\t", segments[ii].caUnits[jj].Buff_Ca_SRB, segments[ii].caUnits[jj].Buff_Ca_CSQN);
		
				fprintf(outputfile, "%-8e\t%-8e\t", segments[ii].caUnits[jj].Ca_i, segments[ii].caUnits[jj].Ca_sr);
				fprintf(outputfile, "%-8e\t", segments[ii].caUnits[jj].Buff_Ca_EGTA_cyt);
				
				fprintf(outputfile, "%-8e\t%-8e\t%-8e\t", segments[ii].caUnits[jj].Ca_srs, segments[ii].caUnits[jj].Buff_Ca_SLLsrs, segments[ii].caUnits[jj].Buff_Ca_SLHsrs);
				fprintf(outputfile, "%-8e\t", segments[ii].caUnits[jj].RyR_Integr);
			}
			
			// Other stuff
			fprintf(outputfile, "%-8e\t%-8e\t", segments[ii].Na_i, segments[ii].K_i);
		}
	}
}

void RAMM::processX0DataForMarkovModel(double* pX0Data, int sindex, double* MM, int num_states, int model_type, int num_channels)
{
	double statesum = 0;	
	for(int sk=0; sk<num_states; sk++) { MM[sk] = pX0Data[sindex+sk]; statesum += MM[sk]; }
	if(model_type == 3 && statesum < 1.1) 
	{
		cout << "Input is deterministic, target stochastic. Adapting to " << num_channels << " channels." << endl;
		statesum = 0; 
		for(int sk=0; sk<num_states-1; sk++) 
		{ 
			MM[sk] = max((int)floor(num_channels * MM[sk]),0); 
			statesum += MM[sk]; 
		}
		MM[num_states-1] = num_channels - statesum;
	}
	else
	{
		if(model_type < 3 && statesum > 2.1) // Deterministic model but stochastic input.
		{
			cout << "Input is stochastic, target deterministic. Scaling by " << statesum << "." << endl;
			for(int sk=0; sk<num_states; sk++)
			{
				MM[sk] = 1.0 * MM[sk] / statesum;
			}
		}
	}
}

void RAMM::processX0DataForMarkovModel(double* pX0Data, int sindex, int* MM, int num_states, int model_type, int num_channels)
{
	double statesum = 0;	
	for(int sk=0; sk<num_states; sk++) { MM[sk] = (int)(pX0Data[sindex+sk]); statesum += MM[sk]; }
	
	if(statesum != num_channels)
	{			
		double fact = (1.0 * num_channels) / statesum;
		int newsum = 0;
		
		cout << "WARNING: Num channels (" << num_channels << ") does not match state sum (" << statesum << "), scaling by " << fact << endl;
		GenPurpose::wait(5.0);
		
		for(int sk=0; sk<num_states-1; sk++)
		{
			MM[sk] = (int)(MM[sk] * fact);
			newsum += MM[sk];
		}
		MM[num_states-1] = num_channels - newsum;
	}
}

int RAMM::lengthMinimumStateVector()
{
	// Total number of variables:
	//		Global electro variables:	  1 + 2 x N_Segments
	//		Membrane electro variables:  37 x 2 x N_Segments
	//		Local CaRU variables:        19 x N_Segments x N_CaDomains
	//		Global signaling variables:  00
	//									-------------
	//									 1 + 76 x N_Segments + 19 x N_Segments x N_CaDomains
	
	return 1 + 76 + 2 * 19;
}

void RAMM::loadX0Data(double* x0data, int x0length)
{
	int startat = -1;
	int numVarsPerMembr = 25;
	int numVarsPerSegment = 2 + 2 * numVarsPerMembr + 23 * settings->N_CaDomains;
	
	Vm = x0data[startat+1];
	
	if(x0length >= (1 + settings->N_Segments*numVarsPerSegment))
	{
		cout << "Using data from every individual Segment to populate " << settings->N_Segments << " Segments" << endl;
	}
	else
	{
		if(x0length >= 1 + numVarsPerSegment)
		{
			cout << "Using data from first Segment to populate all " << settings->N_Segments << " Segments" << endl;
			numVarsPerSegment = 0;
		}
		else
		{
			cout << "ERROR: X0 is of unsufficient length (" << x0length << "), minimum = " << numVarsPerSegment << endl;
			return;
		}
	}
	
	for(int ii=0; ii<settings->N_Segments; ii++)
	{	
		// Loop through the two membrane sides
		for(int jj=0; jj<2; jj++)
		{
			segments[ii].membrane[jj].m_Na = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+2];
			segments[ii].membrane[jj].h1_Na = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+3];
			segments[ii].membrane[jj].h2_Na = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+4];
			segments[ii].membrane[jj].to_r = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+5]; 
			segments[ii].membrane[jj].to_s1 = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+6]; 
			segments[ii].membrane[jj].to_s2 = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+7];
			segments[ii].membrane[jj].to_s3 = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+8];
			segments[ii].membrane[jj].pa_Kr = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+9]; 
			segments[ii].membrane[jj].pi_Kr = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+10];
			segments[ii].membrane[jj].z_Ks = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+11];
			segments[ii].membrane[jj].d_L = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+12];
			segments[ii].membrane[jj].f_L = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+13];
			segments[ii].membrane[jj].d_T = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+14]; 
			segments[ii].membrane[jj].f_T = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+15];

			segments[ii].membrane[jj].Buff_Na_junc = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+16];
			segments[ii].membrane[jj].Buff_Na_sl = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+17];
			segments[ii].membrane[jj].Buff_Ca_SLLj = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+18];
			segments[ii].membrane[jj].Buff_Ca_SLLsl = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+19];
			segments[ii].membrane[jj].Buff_Ca_SLHj = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+20];
			segments[ii].membrane[jj].Buff_Ca_SLHsl = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+21];
			segments[ii].membrane[jj].Buff_Ca_EGTA_sl = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+22];
			
			segments[ii].membrane[jj].Ca_junc = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+23];
			segments[ii].membrane[jj].Ca_sl = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+24];
			segments[ii].membrane[jj].Na_junc = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+25];
			segments[ii].membrane[jj].Na_sl = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+26];
			
			//for(int abcd=0; abcd<8; abcd++)
			//{
			//	segments[ii].membrane[jj].ICaLState_NP[abcd] = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+31+abcd];
			//}
			
			//for(int abcd=0; abcd<4; abcd++)
			//{
			//	segments[ii].membrane[jj].SK_State[abcd] = x0data[ii*numVarsPerSegment+jj*numVarsPerMembr+startat+39+abcd];
			//}
		}
		
		// Loop through the cytosolic Ca2+ domains
		for(int jj=0; jj<settings->N_CaDomains; jj++)
		{
			segments[ii].caUnits[jj].RyRState_NP[0] = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+2]; 
			segments[ii].caUnits[jj].RyRState_NP[1] = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+3];
			segments[ii].caUnits[jj].RyRState_NP[2] = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+4];
			segments[ii].caUnits[jj].RyRState_NP[3] = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+5];
			segments[ii].caUnits[jj].RyRState_NP[4] = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+6]; 
			segments[ii].caUnits[jj].RyRState_NP[5] = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+7];
			segments[ii].caUnits[jj].RyRState_NP[6] = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+8];
			segments[ii].caUnits[jj].RyRState_NP[7] = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+9];
			
			double sumRyR = segments[ii].caUnits[jj].RyRState_NP[0] + segments[ii].caUnits[jj].RyRState_NP[1] + segments[ii].caUnits[jj].RyRState_NP[2] + segments[ii].caUnits[jj].RyRState_NP[3] + segments[ii].caUnits[jj].RyRState_NP[4] + segments[ii].caUnits[jj].RyRState_NP[5] + segments[ii].caUnits[jj].RyRState_NP[6] + segments[ii].caUnits[jj].RyRState_NP[7];
			if(settings->RyRModelType == 3)
			{
				if(fabs(sumRyR * (settings->N_Segments * settings->N_CaDomains) - settings->numRyRs) > 0.001)
				{
					double fact = settings->numRyRs * 1.0 / (settings->N_Segments * settings->N_CaDomains);					
					cout << "Sum RyR in unit [" << ii << "," << jj << "] = " << sumRyR << ", rescaling to " << fact << ", old dist = [" << segments[ii].caUnits[jj].RyRState_NP[0] << "," << segments[ii].caUnits[jj].RyRState_NP[1] << "," << segments[ii].caUnits[jj].RyRState_NP[2] << "," << segments[ii].caUnits[jj].RyRState_NP[3] << "]";
					segments[ii].caUnits[jj].RyRState_NP[0] = (int)((fact * segments[ii].caUnits[jj].RyRState_NP[0]) / sumRyR);
					segments[ii].caUnits[jj].RyRState_NP[1] = (int)((fact * segments[ii].caUnits[jj].RyRState_NP[1]) / sumRyR);
					segments[ii].caUnits[jj].RyRState_NP[2] = (int)((fact * segments[ii].caUnits[jj].RyRState_NP[2]) / sumRyR);
					segments[ii].caUnits[jj].RyRState_NP[3] = (int)((fact * segments[ii].caUnits[jj].RyRState_NP[3]) / sumRyR);					
					segments[ii].caUnits[jj].RyRState_NP[4] = (int)((fact * segments[ii].caUnits[jj].RyRState_NP[4]) / sumRyR);					
					segments[ii].caUnits[jj].RyRState_NP[5] = (int)((fact * segments[ii].caUnits[jj].RyRState_NP[5]) / sumRyR);					
					segments[ii].caUnits[jj].RyRState_NP[6] = (int)((fact * segments[ii].caUnits[jj].RyRState_NP[6]) / sumRyR);					
					segments[ii].caUnits[jj].RyRState_NP[7] = (int)((fact * segments[ii].caUnits[jj].RyRState_NP[7]) / sumRyR);					
					cout << ", new dist = [" << segments[ii].caUnits[jj].RyRState_NP[0] << "," << segments[ii].caUnits[jj].RyRState_NP[1] << "," << segments[ii].caUnits[jj].RyRState_NP[2] << "," << segments[ii].caUnits[jj].RyRState_NP[3] << "]" << endl;
				}
			}
			else
			{
				if(settings->RyRModelType == 2 && sumRyR > 2.1)
				{
					cout << "Input is stochastic, target deterministic. Scaling by " << sumRyR;
					segments[ii].caUnits[jj].RyRState_NP[0] = segments[ii].caUnits[jj].RyRState_NP[0] / sumRyR;
					segments[ii].caUnits[jj].RyRState_NP[1] = segments[ii].caUnits[jj].RyRState_NP[1] / sumRyR;
					segments[ii].caUnits[jj].RyRState_NP[2] = segments[ii].caUnits[jj].RyRState_NP[2] / sumRyR;
					segments[ii].caUnits[jj].RyRState_NP[3] = segments[ii].caUnits[jj].RyRState_NP[3] / sumRyR;
					segments[ii].caUnits[jj].RyRState_NP[4] = segments[ii].caUnits[jj].RyRState_NP[4] / sumRyR;
					segments[ii].caUnits[jj].RyRState_NP[5] = segments[ii].caUnits[jj].RyRState_NP[5] / sumRyR;
					segments[ii].caUnits[jj].RyRState_NP[6] = segments[ii].caUnits[jj].RyRState_NP[6] / sumRyR;
					segments[ii].caUnits[jj].RyRState_NP[7] = segments[ii].caUnits[jj].RyRState_NP[7] / sumRyR;
					cout << ", new dist = [" << segments[ii].caUnits[jj].RyRState_NP[0] << "," << segments[ii].caUnits[jj].RyRState_NP[1] << "," << segments[ii].caUnits[jj].RyRState_NP[2] << "," << segments[ii].caUnits[jj].RyRState_NP[3] << "]" << endl;
				}
			}
			
			segments[ii].caUnits[jj].Buff_Ca_TnCL = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+10];
			segments[ii].caUnits[jj].Buff_Ca_TnCHc = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+11];
			segments[ii].caUnits[jj].Buff_Ca_TnCHm = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+12];
			segments[ii].caUnits[jj].Buff_Ca_CaM = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+13];
			segments[ii].caUnits[jj].Buff_Ca_Myosin_ca = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+14];
			segments[ii].caUnits[jj].Buff_Ca_Myosin_mg = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+15];;
			segments[ii].caUnits[jj].Buff_Ca_SRB = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+16];
			segments[ii].caUnits[jj].Buff_Ca_CSQN = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+17];
	
			segments[ii].caUnits[jj].Ca_i = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+18];
			segments[ii].caUnits[jj].Ca_sr = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+19];
			segments[ii].caUnits[jj].Buff_Ca_EGTA_cyt = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+20];
			
			segments[ii].caUnits[jj].Ca_srs = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+21];
			segments[ii].caUnits[jj].Buff_Ca_SLLsrs = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+22];
			segments[ii].caUnits[jj].Buff_Ca_SLHsrs = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+23];
			
			segments[ii].caUnits[jj].RyR_Integr = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+jj*23+startat+24];
		}
		
		// Other stuff
		segments[ii].Na_i = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+settings->N_CaDomains*23+startat+2];
		segments[ii].K_i = x0data[ii*numVarsPerSegment+2*numVarsPerMembr+settings->N_CaDomains*23+startat+3];
	}
}
