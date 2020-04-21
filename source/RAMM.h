#include <math.h>
#include <stdio.h>
#include <iostream>
#include "genpurpose.h"
#include "constants.h"
#include "settings.h"
#include "MersenneTwister.h"

class RAMM
{
	public:
		bool showlog;
		double caiont, kiont, naiont, clont;
		
		double fTnIP, fINaP, fINaKP, fIKsP, fICaLP, fRyRP, fPLBP, fIKurP;
		double fPLBP_CaMKII, fRyRP_CaMKII, fIToP_CaMKII, fICaLP_CaMKII, fINaP_CaMKII, fIK1P_CaMKII;

		double Vm;
		
		double Na_i, Na_junc, Na_sl, K_i;
		double Ca_junc, Ca_sl, Ca_sr, Ca_srs, Ca_i;
		double Cl_i, Cl_o;

		double INa, INa_junc, INa_sl, INaK, INaK_junc, INaK_sl, INab, INab_junc, INab_sl, IKr, IKs, IK1, Isus, ITo, ISK_junc, ISK_sl, ICaP_junc, ICaP_sl, IClb, IKAch, IClCa_junc, IClCa_sl;
		double ICaL_junc, ICaL_sl, ICaT_junc, ICaT_sl, INaCa_junc, INaCa_sl, ICab_junc, ICab_sl; // Whole-cell equivalents of CaRU currents
		double IpCa_junc, IpCa_sl; 
		double Jrel, Jleak, Jup, Buff_Ca_EGTA_cyt;
		double Istim;
		double inak_junc, inak_sl, isus, iClb;
		double E_Cl, EB_Cl;

		Constants* data;
		Settings* settings;

		struct CaDomain
		{
			double RyRState_NP[8];
			double RyR_SEP;
			double RyR_Integr;
			double Ca_i, Ca_sr, Ca_srs;
			double Buff_Ca_TnCL, Buff_Ca_TnCHc, Buff_Ca_TnCHm, Buff_Ca_CaM, Buff_Ca_Myosin_ca, Buff_Ca_Myosin_mg, Buff_Ca_SRB, Buff_Ca_CSQN;
			double Buff_Ca_SLLsrs, Buff_Ca_SLHsrs;
			double Buff_Ca_EGTA_cyt;
			
			double Jrel, Jleak, Jup;
			double J_diff_i_subcell;
		};
		
		struct MembrDomain
		{
			double Ca_junc, Ca_sl;
			double Na_junc, Na_sl;
			
			//double d_CaL_NP, d_CaL_P, f_CaL_NP, f_CaL_P, f_CaL_Ca_junc, f_CaL_Ca_sl;
			//double ICaLState_NP[8];
			//double SK_State[4];
			double m_Na, h1_Na, h2_Na, to_r, to_s1, to_s2, to_s3, z_Ks, pa_Kr, pi_Kr, d_L, f_L, d_T, f_T, O_KCa_junc, O_KCa_sl, a_KAch; //gates
			
			double Buff_Ca_SLLj, Buff_Ca_SLLsl, Buff_Ca_SLHj, Buff_Ca_SLHsl;
			double Buff_Na_junc, Buff_Na_sl, Buff_Ca_EGTA_sl;
			
			double INa_junc, INa_sl, INaK_junc, INaK_sl, INab_junc, INaCa_junc, INaCa_sl, INab_sl, IKAch;
			double IK1, Isus, IKs, IKr, ITo, ISK_junc, ISK_sl, IClb, IClCa_junc, IClCa_sl;
			
			double ICab_junc, ICab_sl, ICaP_junc, ICaP_sl;
			double ICaL_junc, ICaL_sl, ICaL, ICaL_tot, ICaT_sl, ICaT_junc, ICaT_tot, ICaT;
			double IpCa_junc, IpCa_sl;
		};
		
		struct Segment
		{
			MembrDomain membrane[2];
			CaDomain* caUnits;
			double Na_i, K_i, Cl_i;
		};
		
		Segment* segments;
		
		double** Copy_Ca_i;
		double** Copy_Ca_sr;
		double** Copy_Ca_srs;
		double** Copy_Ca_sl;
		double** Copy_Na_sl;
		double* Copy_Na_i;
		
		double vrest, vmax;
		double tmaxdvdt, valmaxdvdt, Max_State_Change;
		double APD, CaT_Min, CaT_Max, Vmin, APA;

		void init(Constants *idata, Settings *isettings);
		void update(double t, double dt);
		void update(int cellnumber, double t, double dt, double vleft, double vright);
		void outputData(FILE* outputfile, int levelOfElectroDetail, int levelOfSignalingDetail, bool binaryOutput);
		void update(int cellnumber, double t, double dt, double vleft, double vright, double vup, double vdown);
		RAMM* clone();

		void loadX0Data(double* x0data, int x0length);
		int lengthMinimumStateVector();
		
	private:
		MTRand rvgen;	
		
		double sigma, fnak; // Helper variables so that these only have to be calculated once per timestep
		double h_infinity, tauh1, tauh2, am, bm, ah, bh;
		double alpha_r, beta_r, rss, taur, s1ss, taus1, s2ss, taus2, s3ss, taus3;
		//double inak_junc, inak_sl, icat, ical, dprime, icat;
		double alpha_z, beta_z, zss_Ks, tauz_Ks, alpha_pa, beta_pa, pass_Kr, taupa_Kr, alpha_pi, beta_pi, piss_Kr, taupi_Kr, alpha_a, beta_a;
		double E0_alpha_d_L, E0_beta_d_L, E10, alpha_d_L, beta_d_L, d_L_infinity, tau_d_L, d_L, E0_f_L, alpha_f_L, beta_f_L, f_L_infinity, tau_f_L;
		double E0_d_T, alpha_d_T, beta_d_T, d_T_infinity, tau_d_T, E0_f_T, alpha_f_T,beta_f_T, f_T_infinity, tau_f_T;

		void processX0DataForMarkovModel(double* pX0Data, int sindex, double* MM, int num_states, int model_type, int num_channels);
		void processX0DataForMarkovModel(double* pX0Data, int sindex, int* MM, int num_states, int model_type, int num_channels);
		
		void updateCaDomain(double dt, double t, int segment, int ca_domain, double Ca_i_L_Domain, double Ca_i_R_Domain, double Ca_i_U_Segment, double Ca_i_D_Segment, double Ca_sr_L_Domain, double Ca_sr_R_Domain, double Ca_sr_U_Segment, double Ca_sr_D_Segment, double Ca_srs_L_Domain, double Ca_srs_R_Domain, double Ca_srs_U_Segment, double Ca_srs_D_Segment);
				
		void Ion_HH_Na(int seg_i, int mem_i, double dt);		
		void Ion_HH_NaL(int seg_i, int mem_i, double dt);
		void Ion_HH_NaCa(int seg_i, int mem_i, double dt);
		void Ion_HH_NaK(int seg_i, int mem_i, double dt);
		void Ion_HH_Nab(int seg_i, int mem_i, double dt);
		
		void Ion_HH_K1(int seg_i, int mem_i, double dt);
		void Ion_HH_Kp(int seg_i, int mem_i, double dt);
		//void Ion_HH_KAch(int seg_i, int mem_i, double dt);
		void Ion_HH_sus(int seg_i, int mem_i, double dt);
		void Ion_HH_Kr(int seg_i, int mem_i, double dt);
		void Ion_HH_To(int seg_i, int mem_i, double dt);
		void Ion_HH_SK(int seg_i, int mem_i, double dt);
		void Ion_HH_KAch(int seg_i, int mem_i, double dt);
		void Ion_HH_Clb(int seg_i, int mem_i, double dt);
		void Ion_HH_ClCa(int seg_i, int mem_i, double dt);		
		//void Ion_HH_Cl(int seg_i, int mem_i, double dt);
		void Ion_HH_Ks(int seg_i, int mem_i, double dt);
		void Ion_HH_CaL(int seg_i, int mem_i, double dt);
		void Ion_HH_CaT(int seg_i, int mem_i, double dt);
		void Ion_HH_CaP(int seg_i, int mem_i, double dt);
		void Ion_HH_Cab(int seg_i, int mem_i, double dt);
		//void Ion_HH_CaL_Alt(int seg_i, int mem_i, double dt);
		//void Ion_Markov_CaL(int seg_i, int mem_i, double dt);
		void Ion_Ca_Handling(int seg_i, int mem_i, double dt);
		void Ca_CICR(int seg_i, int cadom_i, double dt, double t);
		void Ca_CICR_Stoch(int seg_i, int cadom_i, double dt, double t);
		void Ca_CICR_Stoch_Sobie(int seg_i, int cadom_i, double dt, double t);
		void UpdateConcentrations(int seg_i, int mem_i, double dt, double t, double Na_sl_U_Segment, double Na_sl_L_Segment, double Ca_sl_U_Segment, double Ca_sl_L_Segment, double Ca_i_neighbor, double Ca_srs_neighbor);

		double Buff_Ca_DoubleBuff(double ca_t, double buff1_tot, double buff2_tot, double Kmbuff1, double Kmbuff2);
		
		double check_State_Change(double oldval, double dval, double dt);
};
