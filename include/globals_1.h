#ifndef GLOBALS_1_H
#define GLOBALS_1_H

#include <string>
#include <vector>

extern bool log_print_mode;
extern int  external_random_seed;

// ===========================================================================
// (1) OPINION DYNAMICS (PART A) PARAMETERS **********************************
// ===========================================================================
// --[Parameters Group 1]- - - - - - - - - - - - - - - - - - - - - - - - - - -
extern int t_0, t_1, t_2; 
// --[Parameters Group 2]- - - - - - - - - - - - - - - - - - - - - - - - - - -
extern int S_1, S_2, t_3, t_4; 
// --[Parameters Group 3]- - - - - - - - - - - - - - - - - - - - - - - - - - -
extern double alpha_le64, alpha_ge65;  
extern double C_p, C_n, C_q, AbsoluteLoss, UPOWV;
// --[Parameters Group 4]- - - - - - - - - - - - - - - - - - - - - - - - - - -
extern double u_c2_P00, u_c2_Ps;  
extern double learn_P00, learn_u00, learn_Ps, learn_us;
extern double u_ind_P00, u_ind_Ps;
extern std::vector<double> c_P00, c_u00, c_Ps, c_us;
// --[Parameters Group 5]- - - - - - - - - - - - - - - - - - - - - - - - - - -
extern double P10, P01, P11;
// --[Parameters Group 6]- - - - - - - - - - - - - - - - - - - - - - - - - - -
extern double sci_1; 
extern int    sci_2a;
extern double sci_2b;
extern int    sci_3a;
extern double sci_3b;

// ===========================================================================
// (2) DISEASE DYNAMICS (PART B) PARAMETERS **********************************
// ===========================================================================
extern double sigma_val, gamma_val, trans_val, EFF, I0v, ChR, omega_val;
extern int    ChT0, ChT1;

// ===========================================================================
// (3) OTHER HARD-CODED VARIABLES
// ===========================================================================
extern std::string inputDirectory, outputDirectory;
extern int n_reps, write_log_flag, N, T;
extern std::vector<std::string> argNames;
extern std::vector<std::string> dirs;
extern std::string argStr; 

#endif

