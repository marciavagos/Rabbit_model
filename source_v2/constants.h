#include <math.h>

class Constants
{
	public:		        
		bool cellType, isAF, isRA;
		
    // General parameters
    // ================================================
		double F;   // Faraday constant
		double R;   // Gas constant
		double T;	 // Temperature
		
		double l, a, vcell, ageo, Acap, vmyo, vsr, vsl, vjunc, AF, frt, RTF; 
		double Fjunc, Fjunc_CaL;
		double Cmem, CmemF;
		double N_A, nN_A;
		double Rmyo;
    
    // Sodium
    // ================================================
		double GNa, GNab, IbarNaK, KmNaip_NP, KmNaip_P, KmK1, steepK1, shiftK1, Q10NaK, Q10KmNai,iNaKmax, kNaK_K, kNaK_Na_n, kNaK_Na;
		double J_na_slmyo, J_na_juncsl;

    // Calcium
    // ================================================
		//double p_ICaL_Ca, p_ICaL_Na, p_ICaL_K, Q10CaL, 
		double IbarNCX, KmCai, KmCao, KmNai, KmNao, KmNao_p3, ksat, Kdact, nu;
		double GCaL, GCaT;
		double ICaP_max, k_CaP;
		double kNaCa, dNaCa, n_NaCa, gamma;
		double Q10NCX, IbarSLCaP, KmPCa, KmPCaPow1_6, GCab, Q10SLCaP, Q10SRCaP;
		double Kmf_NP, Kmf_P, Kmr, hillSRCaP, Vmax_SRCaP;
		
		double Bmax_Naj, Bmax_Nasl, koff_na, kon_na, Bmax_TnClow, koff_tncl, kon_tncl, Bmax_TnChigh, koff_tnchca;
		double kon_tnchca, koff_tnchmg, kon_tnchmg, Bmax_CaM, koff_cam, kon_cam, Bmax_myosin, koff_myoca, kon_myoca;
		double koff_myomg, kon_myomg, Bmax_SR, koff_sr, kon_sr, Bmax_SLlowsl, Bmax_SLlowj, koff_sll, kon_sll, Bmax_SLhighsl;
		double Bmax_SLhighj, koff_slh, kon_slh, Bmax_Csqn, koff_csqn, kon_csqn, kon_EGTA, koff_EGTA;
		double EGTA;
		
		double J_ca_slmyo, J_ca_juncsl;
		double tau_diff_Cai_Segment, tau_diff_Cai_Domain, tau_diff_Casr_Segment, tau_diff_Casr_Domain, tau_diff_Casl_Segment;
		double tau_diff_SRS_i, tau_diff_Casrs_Domain, tau_diff_Casrs_Segment;
		
    // Potassium
    // ================================================
		double pNaK, GKp, GK1max, GToFast, Gsus, GKr, GKs, GKCa, GKAch;
		double KCa_on, KCa_off;

    // Chloride
    // ================================================
		double GClb, GClCa, KdClCa; 
		
		Constants()
		{
			Constants(false, false, false);
		}
		
		Constants(bool isEndo, bool l_isAF, bool l_isRA)
		{
			cellType = isEndo;
			isAF = l_isAF;
			isRA = l_isRA;
					
			Rmyo = 150;
			
			R = 8314;
			F = 96485;
			T = 308; //35
			N_A = 6.0221415E23;
			nN_A = 6.0221415E14;

			l = 100; // [um]
			a = 10.25; // [um]
	
			vcell = 1.26E-14; // [L]          
			vmyo = vcell*0.65;     
			vsr = vcell*0.035;      
			vsl = vcell*0.02;   
			vjunc = vcell*0.001;
			
			Fjunc = 0.33;
			Fjunc_CaL = 0.9;
			
			ageo = 2*M_PI*a*a+2*M_PI*a*l; 
			Acap = ageo*2;     
			AF=Acap/F;
			frt=F/T/R;
			RTF=1/frt;
			
			Cmem = 0.5E-10; //[picoF]
			CmemF = Cmem / F;
			
			Gsus = 26E-3; //[nS/pF]
			GK1max = 203.2E-3; //[nS/pF]
			GNa = 2.8E-5; //[microL/(sec*pF)]
			//P_Na = 2.8E-3;  //[nL/(sec*pF)]
			GToFast = 200E-3; //[nS/pF]
                        GKs = 50E-3; //2.5/50  //[nS/pF]
                        GKr = 70E-3; //3.5/50 //[nS/pF]
			GKCa = 0.7E-3; //[nS/pF]
			GKAch = 0.157E-3; //[nS/pF]
			GClb = 2.4E-3; //[nS/pF]
			GClCa = 54.8E-3; 
			KCa_on = 47E6;
			KCa_off = 13;
			KdClCa = 100E-3;

			if(isAF)
			{
				GK1max = 1.8 * GK1max;
				if(isRA)
				{
					Gsus = 0.5*Gsus*1.2; // (1+2*ISO)*
				}
				else
				{
					Gsus = 0.5*Gsus; // (1+2*ISO)*
				}
			}
			else
			{
				if(isRA)
				{
					Gsus = Gsus*1.2; // (1+2*ISO)*
				}
				else
				{
					Gsus = Gsus*1.0; // (1+2*ISO)*
				}
			}
			
			GNab = 0.4E-3;
			IbarNaK = 1.0;
			KmNaip_NP = 11.0;
			KmNaip_P = 0.75 * 11.0;
			KmK1 = 0.59;
			steepK1 = 1.393;
			shiftK1 = -3.6;
			Q10NaK = 1.63;
			Q10KmNai = 1.39;
			iNaKmax = 1.288; //[pA/pF]
			kNaK_K = 1;
			kNaK_Na_n = 1.5;
			kNaK_Na = 11;

			pNaK = 0.01833;
			
			J_ca_juncsl = 0.01 * 1/1.2134e12; // [L/msec] = 8.2413e-13
			J_ca_slmyo = 2 * 1/2.68510e11; // [L/msec] = 3.2743e-12
			J_na_juncsl = 1/(1.6382e12/3*100); // [L/msec] = 6.1043e-13
			J_na_slmyo = 1/(1.8308e10/3*100);  // [L/msec] = 5.4621e-11
			
			GCaL = 144E-3;
			GCaT = 120E-3;
			gamma = 0.45;
			dNaCa = 3e-4;
			kNaCa = 1.41E-4;
			n_NaCa = 3;
			ICaP_max = 190.0E-3; 
			k_CaP = 2e-4;
			Q10NCX = 1.57;
			IbarSLCaP = 0.0471;
			KmPCa = 0.5E-3;
			KmPCaPow1_6 = pow(KmPCa,1.6);
			GCab = 0.4E-3;
			Q10SLCaP = 2.35;
			
			KmCai = 3.59E-3;
                        KmCao = 1.3;
                        KmNai = 12.29;
                        KmNao = 87.5;
                        KmNao_p3 = pow(KmNao,3);
                        ksat = 0.27;
                        IbarNCX = 1.4 * 3.15;
                        Kdact = 0.384E-3;
                        nu = 0.35;
		
			Q10SRCaP = 2.6;
			Vmax_SRCaP = 5.31114E-3;
			Kmf_NP = 0.000625;
			Kmf_P = 1.25 * 0.246E-3;
			Kmr = 1.0;
			hillSRCaP = 1.787;

			// Buffering parameters
			// koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
			Bmax_Naj = 7.561;       // [mM] % Na buffering
			Bmax_Nasl = 1.65;       // [mM]
			koff_na = 1e-3;
			kon_na = 0.1e-3;
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
			Bmax_SLlowsl = 37.4e-3*vmyo/vsl;        // [mM]    % SL buffering
			Bmax_SLlowj = 4.6e-3*vmyo/vjunc*0.1;    // [mM]    %Fei *0.1!!! junction reduction factor
			koff_sll = 1300e-3;
			kon_sll = 100;
			Bmax_SLhighsl = 13.4e-3*vmyo/vsl;
			Bmax_SLhighj = 1.65e-3*vmyo/vjunc*0.1;  // [mM] %Fei *0.1!!! junction reduction factor
			koff_slh = 30e-3;
			kon_slh = 100;
			Bmax_Csqn = 140e-3*vmyo/vsr;            // [mM] % Bmax_Csqn = 2.6;      % Csqn buffering
			koff_csqn = 65;
			kon_csqn = 100;
			kon_EGTA = 5.0; // Smith et al. Analytical Biochemistry 1984
			koff_EGTA = 7.5E-4; // Smith et al. Analytical Biochemistry 1984
			EGTA = 0.0;

			tau_diff_SRS_i = 12;
			tau_diff_Casl_Segment = 3.4; tau_diff_Casrs_Segment = 0.3; tau_diff_Casrs_Domain = 0.125;
			tau_diff_Cai_Segment = 0.6; tau_diff_Cai_Domain = 0.6;
			tau_diff_Casr_Segment = 15.0; tau_diff_Casr_Domain = 15;
		}
		
		~Constants()
		{
		
		}
};
