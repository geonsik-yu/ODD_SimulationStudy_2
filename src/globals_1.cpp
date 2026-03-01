// ===========================================================================
// <<  GLOBAL VARIABLES 1. MODEL PARAMETERS >>
// ===========================================================================
#include "include/globals_1.h"

bool log_print_mode = false;
int external_random_seed = 0; // argv[40] External random seed value
// ===========================================================================
// (1) OPINION DYNAMICS (PART A) PARAMETERS **********************************
// ===========================================================================
// --[Parameters Group 1]- - - - - - - - - - - - - - - - - - - - - - - - - - -
//   >> Timeline line parameters regarding NPI policy
int t_0 = 0;
int t_1; // argv[1] (X1) The date on which penalties for not complying with NPIs began to be imposed.
int t_2; // argv[2] (X2) The date on which penalties for not complying with NPIs stopped to be imposed for **< those who are vaccinated >**.


// --[Parameters Group 2]- - - - - - - - - - - - - - - - - - - - - - - - - - -
//   >> Timeline line parameters regarding vaccine admin capacity & phase timeline
int S_1; // argv[3] (X3) Daily vaccine administration capacity during phase 1
int S_2; // argv[4] (X4) Daily vaccine administration capacity during phase 2
int t_3 = 255; // The date on which the vaccination phase 1 started.    
int t_4 = 343; // The date on which the vaccination phase 2 started. 


// --[Parameters Group 3]- - - - - - - - - - - - - - - - - - - - - - - - - - -
//   >> Recourse level modifier for infectiousness, other cost parameters, and utility power.
//      We assume xi_le64 = alpha_le64  and  xi_ge65 = alpha_ge65
double alpha_le64;  // argv[5]  (X5) Resource level modifier for getting infectious given agent k’s age is smaller than or equal to 64.
double alpha_ge65;  // argv[6]  (X6) Resource level modifier for getting infectious given agent k’s age is greater than or equal to 65.
double C_p;         // argv[7]  (X7) The social penalty cost of non-compliance with NPIs 
double C_n;         // argv[8]  (X8) The cost of NPI compliance, such as inconvenience and expenses for the required items.
double C_q;         // argv[9]  (X9) The cost of experiencing congestion in the vaccine queue.
double AbsoluteLoss;// argv[10] (X10) Absoluted monetary burden 
double UPOWV;       // argv[11] (X11) Shape parameter for the utility function


// --[Parameters Group 4]- - - - - - - - - - - - - - - - - - - - - - - - - - -
//   >> Input file indicators for some opinion dynamics parameters
//   >>> Directly used as string inputs from argv. No need to have variables.
//     0.5 (fixed) Avg of the agents' media usage values: m_k
//     0.1 (fixed) Var of the agents' media usage values: m_k
//     argv[12] (X12) Avg of the agents' channel usage for M_1
//     argv[13] (X13) Var of the agents' channel usage for M_1
//     argv[14] (X14) Avg of P(0,0) distribution of the agents. 
//     argv[15] (X15) Var of P(0,0) distribution of the agents. 
//     argv[16] (X16) Avg of P(s) distribution of the agents. 
//     argv[17] (X17) Var of P(s) distribution of the agents. 
//   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//   >> Other opinion dynamics parameters (calibration)
double u_c2_P00; // argv[18] (X18) Uncertainty value assigned to M_2's messages regarding P(0,0) [after t_0]
double u_c2_Ps;  // argv[19] (X19) Uncertainty value assigned to M_2's messages regarding P(s)   [after t_3]
double learn_P00;// argv[20] (X20) P00 learning rate
double learn_u00;// argv[21] (X21) U00 learning rate
double learn_Ps; // argv[22] (X22) Ps learning rate
double learn_us; // argv[23] (X23) Us learning rate
//   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//   >> Other opinion dynamics parameters (non-calibration)
double u_ind_P00 = 1.2; // 
double u_ind_Ps  = 1.2; // 
std::vector<double> c_P00 = { 0.25, 0.00 };      // {Pro, Anti} Channels' P00 values
std::vector<double> c_u00 = { 0.05, 4.00, 4.00}; // {Pro, Anti, Anti-RD} Channels' P00 uncertainties
std::vector<double> c_Ps  = { 0.00, 0.75 };      // {Pro, Anti} Channels' Ps values
std::vector<double> c_us  = { 0.05, 4.00, 4.00}; // {Pro, Anti, Anti-RD} Channels' Ps uncertainties


// --[Parameters Group 5]- - - - - - - - - - - - - - - - - - - - - - - - - - -
//   >> Fixed perceived risk (probability) values.
double P10 = 0.010; // The probability of the agent getting infected, assuming compliance with NPIs.
double P01;      // argv[24] (X24) The probability of the agent getting infected, assuming vaccination but non-compliance with NPIs.
double P11;      // argv[25] (X25) The probability of the agent getting infected, assuming compliance with NPIs and vaccination.


// --[Parameters Group 6]- - - - - - - - - - - - - - - - - - - - - - - - - - -
//   >> Scenario defining parameters.
// *** (Scenario 1): NPI compliance cost reduction
double sci_1 = 0.0; // argv[26] (X26) Default scenario (0) & ranges betweem 0.0 and 1.0
                    //              - ratio value of NPI compliance cost reduction (default scenario -- no reduction)
// *** (Scenario 2): Informational Intervention
int    sci_2a = 0;  // argv[27] (X27) Default scenario (0) --- Informational Intervention 1 (scenario index)
  // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
  // - sci_2a == 0  (given ratio) of the randomly selected individual are affected
  // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
  // - sci_2a == 1  (given ratio) of the randomly selected individual from Age-Q1 demographic group are affected 
  // - sci_2a == 2  (given ratio) of the randomly selected individual from Age-Q2 demographic group are affected 
  // - sci_2a == 3  (given ratio) of the randomly selected individual from Age-Q3 demographic group are affected 
  // - sci_2a == 4  (given ratio) of the randomly selected individual from Age-Q4 demographic group are affected 
  // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
  // - sci_2a == 5  (given ratio) of the randomly selected individual from RA-Q1 demographic group are affected 
  // - sci_2a == 6  (given ratio) of the randomly selected individual from RA-Q2 demographic group are affected 
  // - sci_2a == 7  (given ratio) of the randomly selected individual from RA-Q3 demographic group are affected 
  // - sci_2a == 8  (given ratio) of the randomly selected individual from RA-Q4 demographic group are affected 
  // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
double sci_2b = 0.0; // argv[28] (X28) Default scenario (0) --- Informational Intervention 2 (ratio)
                     //             - sci_2b defines the ratio of the affected population size: Range (0.0, 0.25)
// *** (Scenario 3): Vaccine Eligibility Management
int    sci_3a = 0;   // argv[29] (X29) Default scenario (0) 
  // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
  // - sci_3a == 0  Default scenario -- Age-based rolling eligibility only
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  //                 ** Age-based rolling eligibility up to (given ratio of admin capacity)
  // - sci_3a == 1    >>  Then￼apply [ Open Eligibility ] for the remaining capacity.
  // - sci_3a == 2    >>  Then￼apply [ Exposure-based Eligibility ] for the remaining capacity. (Contact Net Degree)
  // - sci_3a == 3    >>  Then￼apply [ Targeted access for disadvantaged population ] for the remaining capacity.
  // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
double sci_3b = 1.0; // argv[30] (X30) Default scenario (1.0) 
                     //              defines the capacity ratio for the primany vaccine admin capacity: Range (0.4, 1.0)



// ===========================================================================
// (2) DISEASE DYNAMICS (PART B) PARAMETERS **********************************
// ===========================================================================
double sigma_val = 0.1538; // = 1/6.5
double gamma_val = 0.0553; // = 1/18
double trans_val;  // argv[31] (X31 / D1)  Transmission rate per contact
double EFF;        // argv[32] (X32 / D2)  Effectiveness of compliance to the NPI recommendation
double I0v;        // argv[33] (X33 / D3)  Infectious population at t=0
int    ChT0;       // argv[34] (X34 / D4)  Start time of temporary adjustment
int    ChT1;       // argv[35] (X35 / D5)  End time of temporary adjustment
double ChR;        // argv[36] (X36 / D6)  Adjustment rate
double omega_val;  // argv[37] (X37 / D7)  Waning immunity coefficient




// ===========================================================================
// (3) OTHER HARD-CODED VARIABLES ********************************************
// ===========================================================================
// Input and output directories 
std::string inputDirectory = "./";
std::string outputDirectory = "./res/"; 
std::string ccode = "usa";  // Country code (usa:> USA; swe:> Sweden)
int n_reps;        // argv[38]  >>  Number of simulation runs
int write_log_flag;// argv[39]  >>  if 0, this source code does not write log files for 0th instance.
                   //           >>  if 1, this source code writes log files for 0th intance. 
int N = 20000; // Simulation's population size;
int T = 450;   // Simulated number of days (day 0 = Apr 1st, 2020); >> Changed from 550 to 450

// Argument names
std::vector<std::string> argNames =
      { "Simulator",
        //====================================================================
        // [Calibration Parameters from TABLE 1]
        "argv[1] (X1) The date on which NPI-penalties began to be imposed",
        "argv[2] (X2) The date on which NPI-penalties stopped to be imposed",
        "argv[3] (X3) Daily vaccine administration capacity during phase 1",
        "argv[4] (X4) Daily vaccine administration capacity during phase 2",
        //====================================================================        
        // [Calibration Parameters from TABLE 2]
        "argv[5] (X5) Resource level modifier for getting infectious (<65)",
        "argv[6] (X6) Resource level modifier for getting infectious (>=65)", 
        "argv[7] (X7) Cost of (0,0) or (0,1)",
        "argv[8] (X8) Cost of NPI compliance",
        "argv[9] (X9) Cost of congestion in the queue",
        "argv[10] (X10) Absolute loss by infection (value at default scenario = 0",
        "argv[11] (X11) Utility function power value",
        //-------------------------------------------------------------------        
        "argv[12] (X12) Avg of the channel usage for M_1",
        "argv[13] (X13) Var of the channel usage for M_1",
        "argv[14] (X14) Avg of P(0,0) distribution",
        "argv[15] (X15) Var of P(0,0) distribution",
        "argv[16] (X16) Avg of P(s) distribution",
        "argv[17] (X17) Var of P(s) distribution",
        // [Calibration Parameters from TABLE 3]
        //====================================================================
        "argv[18] (X18) Uncertainty - M_2 messages regarding P(0,0)",
        "argv[19] (X19) Uncertainty - M_2 messages regarding P(s)",
        "argv[20] (X20) P00 learning rate",
        "argv[21] (X11) U00 learning rate",
        "argv[22] (X22) Ps learning rate",
        "argv[23] (X23) Us learning rate",
        "argv[24] (X24) P01",
        "argv[25] (X25) P11",
        //====================================================================
        "argv[26] (X26) Scenario Input 1",
        "argv[27] (X27) Scenario Input 2a",
        "argv[28] (X28) Scenario Input 2b",
        "argv[29] (X29) Scenario Input 3a",
        "argv[30] (X30) Scenario Input 3b",
        //====================================================================
        "argv[31] (X31 / D1)  Transmission rate per contact",
        "argv[32] (X32 / D2)  Effectiveness of compliance to the NPI recommendation",
        "argv[33] (X33 / D3)  Infectious population at t=0",
        "argv[34] (X34 / D4)  Start time of temporary adjustment",
        "argv[35] (X35 / D5)  End time of temporary adjustment",
        "argv[36] (X36 / D6)  Adjustment rate",
        "argv[37] (X37 / D7)  1/{expected immunity duration} = Waning immunity coeff",
        "argv[38] Number of simulation runs",
        "argv[39] write_log_flag",
        "argv[40] random seed"};

std::vector<std::string> dirs;
std::string argStr; 

