// ===========================================================================
// C++ Simulator 
// This simulator utilizes nlohmann.json library for json io, 
// which was developed by Niels Lohmann (https://github.com/nlohmann/json).
// 
// This section explains the version information:
// ===========================================================================
// Update Log -- Version 2/2/2026
// - Add input for external seed value assignment for randomization.
// ===========================================================================
// Update Log -- Version 2/22/2025
// ===========================================================================
// Update Log -- Version 2/12/2025
// - Change from static contact network to dynamic network setups
//   >> Load cnets dynamically every iteration.
// - (Memory vs Input File Loading) Consideration Update -- Further Optim ??
//    >> More pre-loaded input files & less file reading procedure calls
//    >> Less pre-lodaed input files & more file reading procedure calls.
// - Add a new epi stat: "# of the new infection among 65+"
// ---------------------------------------------------------------------------
// Update Log -- Version 2/5/2025
// - "Informational intervention" to "reactive devaluation".
// - Add a new epi stat: "# of the vaccinated so far among 65+".
// ---------------------------------------------------------------------------
// Update Log -- Version 12/15/2024
// - Control randomness function update. (split into pieces.)
// ---------------------------------------------------------------------------
// Update Log -- Version 1.3
// - Changed P01: from 0.050 >>> to 0.055
// - Opinion dynamics logic changed (to follow ov-opt_1.3)
//   >> Decision functions updated (all three).
// ---------------------------------------------------------------------------
// Update Log -- Version 1.2
// - "Delta Variant Peak" Implementation
//   Additional input (>> argv[9])  "transmission rate amplifier"
//   Additional input (>> argv[10]) "delta variant period start time"
// - Population characteristic (year, ...) string 
//   >> Modified from argv[9] to argv[11]
// ---------------------------------------------------------------------------
// Update Log -- Version 1.1
// - Implement random initial infections (population with "I" status)
//   randomly select the given number of agents while trying to satisfy
//   the following conditions:
//   (1) Selected agents' average contact degree is close enough to the population's.
//   (2) Selected agents' average income is close enough to the population's.
//   (3) Selected agents' average age is close enough to the population's.
//   >> "Close enough" is defined by the percentage error smaller than 1.0%.
//   >> This process is slow if the threshold percentage error is too small.
// - To test different durations of "waning immunity", another input is now 
//   available >> argv[8].
// - Additional population characteristics can be specified in the input arguments. 
//   >> argv[9]. 


// ===========================================================================
// <<  REQUIRED LIBRARIES  >>
// ===========================================================================
// #include <iomanip>
// #include <numeric>
#include <chrono>
#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <fmt/core.h>

#include "include/globals_1.h"
#include "include/globals_2.h"
#include "include/utils.h"
#include "include/simul-logic.h"

#include "include/json.hpp"  // For json i/o
using json = nlohmann::json;


// ===========================================================================
// <<  PARAMETER INITIALIZATION  >>
// ===========================================================================
bool initialize_globals(int argc, char* argv[]) {
    if (argc < static_cast<int>(argNames.size()) ) {
        std::cerr << "Insufficient arguments.\n";
        return false;
    }

    // Argument string for output file name.
    std::vector<std::string> argVec(argv + 1, argv + argc);
    join(argVec, '_', argStr);

    // Input directory string vector.
    dirs.push_back( fmt::format("{}Data_Calibration/Data_Calibration_usa.json", inputDirectory) );
    dirs.push_back( fmt::format("{}Data_Synthetic/[1,2]_AgeResource/Feature_Age-Resource_usa_2019.json", inputDirectory) );
    dirs.push_back( fmt::format("{}Data_Synthetic/[3]_MediaUsage/Feature_MedUsage_mean=0.5_variance=0.01.json", inputDirectory) );
    dirs.push_back( fmt::format("{}Data_Synthetic/[4]_ChannelUsage/NAD_mean={}_variance={}.json", inputDirectory, argv[12], argv[13]) );
    dirs.push_back( fmt::format("{}Data_Synthetic/[5]_Initial_P00/P00_mean={}_variance={}.json", inputDirectory, argv[14], argv[15]) );
    dirs.push_back( fmt::format("{}Data_Synthetic/[6]_Initial_Ps/Ps_mean={}_variance={}.json", inputDirectory, argv[16], argv[17]) );
    dirs.push_back( fmt::format("{}Data_Synthetic/[7]_OpinionNetwork/OpinionNet_NeiList_10.json", inputDirectory) );

    //------------------------------------------------------------------------
    // [Parameters from TABLE 1]
    t_1 = atoi(argv[1]);// (X1)
    t_2 = atoi(argv[2]);// (X2)
    S_1 = atoi(argv[3]);// (X3)
    S_2 = atoi(argv[4]);// (X4)
    //------------------------------------------------------------------------
    // [Parameters from TABLE 2]
    alpha_le64 = atof(argv[5]);// (X5) 
    alpha_ge65 = atof(argv[6]);// (X6)
    // (Note that) xi_le64 = alpha_le64;  xi_ge65 = alpha_ge65;
    C_p = atof(argv[7]);// (X7) 
    C_n = atof(argv[8]);// (X8) 
    C_q = atof(argv[9]);// (X9) 
    AbsoluteLoss = atof(argv[10]);// (X10)
    UPOWV = atof(argv[11]);// (X11) 
    //------------------------------------------------------------------------
    // [Parameters from TABLE 3]
    u_c2_P00  = atof(argv[18]);// (X18)
    u_c2_Ps   = atof(argv[19]);// (X19)
    c_u00[1] = u_c2_P00;   // {Pro, Anti, Anti-RD} Channels' P00 uncertainties
    c_us[1]  = u_c2_Ps;    // {Pro, Anti, Anti-RD} Channels' Ps uncertainties
    c_u00[2] = u_c2_P00/2; // {Pro, Anti, Anti-RD} Channels' P00 uncertainties
    c_us[2]  = u_c2_Ps/2;  // {Pro, Anti, Anti-RD} Channels' Ps uncertainties
    learn_P00 = atof(argv[20]);// (X20)
    learn_u00 = atof(argv[21]);// (X21)
    learn_Ps  = atof(argv[22]);// (X22)
    learn_us  = atof(argv[23]);// (X23)
    P01   = atof(argv[24]);
    P11   = atof(argv[25]);
    //------------------------------------------------------------------------
    // [Scenario Inputs]
    // < Scenario Type 1) NPI compliance cost reduction >
    sci_1  = atof(argv[26]);// (X26) Default scenario (0) & ranges betweem 0.0 and 1.0
    //  < sci_1 >  - ratio value of NPI compliance cost reduction (default scenario -- no reduction)
    // Apply the NPI compliance cost reduction 
    C_n = C_n * (1.0 - sci_1);
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
    // < Scenario Type 2) Informational Intervention >
    sci_2a = atoi(argv[27]);// (X27) Default scenario (0) --- Informational Intervention 1 (scenario index)
    sci_2b = atof(argv[28]);// (X28) Default scenario (0) --- Informational Intervention 2 (ratio)
    /*
      < sci_2a >
        - sci_2a == 0  (given ratio) of the randomly selected individual are affected
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        - sci_2a == 1  (given ratio) of the randomly selected individual from Age-Q1 demographic group are affected 
        - sci_2a == 2  (given ratio) of the randomly selected individual from Age-Q2 demographic group are affected 
        - sci_2a == 3  (given ratio) of the randomly selected individual from Age-Q3 demographic group are affected 
        - sci_2a == 4  (given ratio) of the randomly selected individual from Age-Q4 demographic group are affected 
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        - sci_2a == 5  (given ratio) of the randomly selected individual from RA-Q1 demographic group are affected 
        - sci_2a == 6  (given ratio) of the randomly selected individual from RA-Q2 demographic group are affected 
        - sci_2a == 7  (given ratio) of the randomly selected individual from RA-Q3 demographic group are affected 
        - sci_2a == 8  (given ratio) of the randomly selected individual from RA-Q4 demographic group are affected 
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        - sci_2a == 9  (given ratio) of the randomly selected individual from RA-Q1 demographic group are affected 
        - sci_2a == 10 (given ratio) of the randomly selected individual from RA-Q2 demographic group are affected 
        - sci_2a == 11 (given ratio) of the randomly selected individual from RA-Q3 demographic group are affected 
        - sci_2a == 12 (given ratio) of the randomly selected individual from RA-Q4 demographic group are affected 
      < sci_2b >
        - sci_2b defines the ratio value & ranges between 0.0 and 0.25
    */
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
    // < Scenario Type 3) Vaccination Eligibility Management >
    sci_3a = atoi(argv[29]);// (X29) Default scenario (0) --- Vaccine Eligibility Management
    sci_3b = atof(argv[30]);// (X30) Default scenario (1.0) 
    /*
      < sci_3a >
        - sci_3a == 0  Default scenario -- Age-based Rolling Eligibility
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
         << Age-based rolling eligibility up to (given ratio of admin capacity) >>
        - sci_3a == 1   >>  Then￼apply [ Open Eligibility ] for the remaining capacity.
        - sci_3a == 2   >>  Then￼apply [ Exposure-based Eligibility ] for the remaining capacity.
        - sci_3a == 3   >>  Then￼apply [ Targeted access for disadvantaged population ] for the remaining capacity.

      < sci_3b >
        - sci_3b defines the ratio value & ranges between 0.4 and 1.0 (Capacity ratio for the primany vaccine admin capacity)
    */
    //------------------------------------------------------------------------
    trans_val = std::stof(argv[31]); // (X31 / D1) Transmission rate per contact
    EFF       = std::stof(argv[32]); // (X32 / D2) Effectiveness of compliance to the NPI recommendation
    I0v       = std::stof(argv[33]); // (X33 / D3) Infectious population at t=0
    ChT0      = std::stoi(argv[34]); // (X34 / D4) Start time of temporary adjustment
    ChT1      = std::stoi(argv[35]); // (X35 / D1) End time of temporary adjustment
    ChR       = std::stof(argv[36]); // (X36 / D5) Adjustment rate
    omega_val = std::stof(argv[37]); // (X37 / D6) 1/{expected immunity duration} = Waning immunity coeff
    //------------------------------------------------------------------------
    n_reps    = std::stoi(argv[38]); // (X38)   
    write_log_flag = std::atoi(argv[39]);
    external_random_seed = std::stoi(argv[40]); // argv[40] External random seed value
    //************************************************************************
    // Print arguments
    /*
    if (true){
        for (int idx = 1; idx <= 11; idx++)
            std::cout << argNames[idx] << ":= " << argv[idx] << std::endl;;
        std::cout << std::endl;
    }
    */
    return true;
}


// ===========================================================================
// <<  MAIN FUNCTION  >>
// ===========================================================================
int main(int argc, char* argv[]) {
    //========================================================================
    // (SETUP 1) Library nlohmann.json setup initialization
    using json = nlohmann::basic_json<std::map, std::vector, std::string, bool, std::int64_t, std::uint64_t, float >;

    // (SETUP 2) Process input arguments.
    if (!initialize_globals(argc, argv)) { return 1; }// Exit if init fails

    // -----------------------------------------------------------------------
    std::discrete_distribution<int> ddist {0.0, 0.0}; // Discrete Random 

    // (SETUP 3) Setup objects for i/o process
    // 1. Standard I/O objects
    std::ofstream ofs;   std::ifstream ifs; 
    // 2. Json objects
    // Temp obj // Compliance(%) save of instances // Vaccination(%) logs of instances // Disease stats logs of instances
    json jObj, JO_clog, JO_vlog, JO_dlog;
    // -----------------------------------------------------------------------
    json JO_com10_age; // for "compliance_per_10perc_age"
    json JO_com10_rl;  // for "compliance_per_10perc_rl"
    json JO_com10_cdeg;// for "compliance_per_10perc_cdeg"
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    json JO_vac10_age; // for "vaccination_per_10perc_age"
    json JO_vac10_rl;  // for "vaccination_per_10perc_rl"
    json JO_vac10_cdeg;// for "vaccination_per_10perc_cdeg"
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    json JO_ne10_age;  // for "newly_exposed_per_10perc_age"
    json JO_ne10_rl;   // for "newly_exposed_per_10perc_rl"
    json JO_ne10_cdeg; // for "newly_exposed_per_10perc_cdeg"
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    json JO_dev10_age; // for "deviation_per_10perc_age"
    json JO_dev10_rl;  // for "deviation_per_10perc_rl"
    json JO_dev10_cdeg;// for "deviation_per_10perc_cdeg"
    //========================================================================
    // (SIM START) Iteratively run the requested number of simul instances
    //========================================================================
    //std::vector<double> scoreVector_V(n_reps, 0.0); // Vaccination fitting score  
    //std::vector<double> scoreVector_C(n_reps, 0.0); // Compliance fitting score
    std::vector<double> scoreVector_D(n_reps, 0.0); // Disease dynamics fitting score
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    double score = 0.0;  // Score for all the instances;
    for (int r_idx = 0; r_idx < n_reps; r_idx++) {
        int temp_add = r_idx + external_random_seed;
        // (SETUP 4) Set a random number generating object // Remove contamination from multiple iterations
        // std::default_random_engine generator{static_cast<long unsigned int>(time(0))}; // Time dependent seed.
        // Fixed seed for reproducibility  ---------------------------------------
        std::default_random_engine gen_00{static_cast<unsigned int>(2024 + temp_add)};  // Initializer (initial infection)
        std::default_random_engine gen_OD{static_cast<unsigned int>(2124 + temp_add)};  // Opinion dynamics section
        std::default_random_engine gen_II{static_cast<unsigned int>(2224 + temp_add)};  // Informational intervention section
        std::default_random_engine gen_VM{static_cast<unsigned int>(2324 + temp_add)};  // Vaccination section
        std::default_random_engine gen_DD{static_cast<unsigned int>(2424 + temp_add)};  // Disease dynamics section        
        std::default_random_engine gen_CN{static_cast<unsigned int>(2524 + temp_add)};  // Dynamic contact network selection  
        //====================================================================
        // (SIMUL INIT 1) Vector Initialization: "STATIC" characteristics (mostly)
        //====================================================================
        // -- 5 variations in a file      
        int r_idx_5 = r_idx % 5;   
        std::vector<double> age(N, 0.0);        // [1] of [1,2]
        std::vector<double> resource(N, 0.0);   // [2] of [1,2]
        std::vector<double> mediaUsg(N, 0.0);   // [3]
        std::vector<double> chan1Usg(N, 0.0);   // [4]
        std::vector<std::vector<int>> oNeiList; // [7] Neighbor lists (onet)
        std::vector<std::vector<int>> cNeiList; // [8] Neighbor lists (cnet)
        std::vector<double> cDegList(N, 0.0);   //   >> initial cnet's degree list.
        //--------------------------------------------------------------------
        // Dynamic Contact Networks' ID Sequence Generation:
        // >> Current Cnet ID = 0th item of the generated random sequence.
        // >> Pre-load 20 different networks and use them 
        std::uniform_int_distribution<int> cnet_id_dist(0, 99); // Total 100 independent network samples are ready -- Number 0 ~ Number 99 (including both ends)
        std::vector<int> cnet_id_sequence;
        cnet_id_sequence.reserve(T);
        for (int i = 0; i < T; ++i) { cnet_id_sequence.push_back(cnet_id_dist(gen_CN)); }
        int current_cnet_id = cnet_id_sequence[0];
        //--------------------------------------------------------------------
        bool cout_flag = false;
        // Initial Contact Network Loading.
        if (cout_flag){ std::cout << "Attempt to access the contact network file. -- " << current_cnet_id << std::endl;}
        ifs.open( fmt::format("{}Data_Synthetic/[8]_ContactNetwork/ContactNet_NeiList_{}.json", inputDirectory, current_cnet_id) ); 
        jObj = json::parse(ifs);   ifs.close(); 
        for (int i = 0; i < N; i++) {  cNeiList.push_back(jObj[i]); }
        for (int i = 0; i < N; i++) {  cDegList[i] = static_cast<double> (cNeiList[i].size());  }            
        //--------------------------------------------------------------------
        if (cout_flag){ std::cout << "Attempt to access " << dirs[1] << std::endl;}
        ifs.open(dirs[1]);   jObj = json::parse(ifs);   ifs.close(); 
        for (int i = 0; i < N; i++) {  age[i] = jObj[r_idx_5][i][0];  resource[i] = jObj[r_idx_5][i][1];  }
        //--------------------------------------------------------------------
        if (cout_flag){ std::cout << "Attempt to access " << dirs[2] << std::endl;}
        ifs.open(dirs[2]);   jObj = json::parse(ifs);   ifs.close(); 
        for (int i = 0; i < N; i++) {  mediaUsg[i] = jObj[r_idx_5][i];  }
        //--------------------------------------------------------------------
        if (cout_flag){ std::cout << "Attempt to access " << dirs[3] << std::endl;}
        ifs.open(dirs[3]);   jObj = json::parse(ifs);   ifs.close(); 
        for (int i = 0; i < N; i++) {  chan1Usg[i] = jObj[r_idx_5][i];  }
        //--------------------------------------------------------------------
        if (cout_flag){ std::cout << "Attempt to access " << dirs[6] << std::endl;}
        ifs.open(dirs[6]);   jObj = json::parse(ifs);   ifs.close(); 
        for (int i = 0; i < N; i++) {  oNeiList.push_back(jObj[r_idx_5][i]);  }

        //====================================================================
        // Define the subpopulation with "Reactive Devaluation" 
        //====================================================================
        std::vector<int> info_intervention_indices; 
        int num_to_select = static_cast<int>(N * sci_2b);
        if (num_to_select > 0){
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Randomizer ( gen_II ) is limited to this section. 
            if (sci_2a == 0){ // Case { 0 }
                // (given ratio) of the randomly selected individual are affected from the entire population
                std::vector<int> numbers(N);
                for (int i = 0; i < N; ++i) {numbers[i] = i;}
                // Shuffle the numbers
                std::shuffle(numbers.begin(), numbers.end(), gen_II);
                info_intervention_indices = std::vector<int>(numbers.begin(), numbers.begin() + num_to_select);
            } else if (sci_2a <= 4) { // Case { 1, 2, 3, 4 }
                //(given ratio) of the randomly selected individual are affected from Age-based demographic group
                /*
                    Case 1) Q0 (Min) ~ Q1       in age
                    Case 2) Q1       ~ Q2       in age
                    Case 3) Q2       ~ Q3       in age
                    Case 4) Q3       ~ Q4 (Max) in age
                */
                int qt_start = sci_2a - 1;
                int qt_end   = sci_2a;
                std::vector<int> qt_vect = find_quantile_indices(age, qt_start, qt_end);
                std::shuffle(qt_vect.begin(), qt_vect.end(), gen_II);
                info_intervention_indices = std::vector<int>(qt_vect.begin(), qt_vect.begin() + num_to_select);
            } else if (sci_2a <= 8) { // Case { 5, 6, 7, 8 }
                //(given ratio) of the randomly selected individual are affected from Resource-availability-based demographic group
                /*
                    Case 5) Q0 (Min) ~ Q1       in resource level
                    Case 6) Q1       ~ Q2       in resource level
                    Case 7) Q2       ~ Q3       in resource level
                    Case 8) Q3       ~ Q4 (Max) in resource level
                */
                int qt_start = sci_2a - 4 - 1;
                int qt_end   = sci_2a - 4;
                std::vector<int> qt_vect = find_quantile_indices(resource, qt_start, qt_end);
                std::shuffle(qt_vect.begin(), qt_vect.end(), gen_II);
                info_intervention_indices = std::vector<int>(qt_vect.begin(), qt_vect.begin() + num_to_select);
            } else { // Case { 9, 10, 11, 12 }
                //(given ratio) of the randomly selected individual are affected from exposure-degree-based demographic group
                /*
                    Case 5) Q0 (Min) ~ Q1       in contact degree
                    Case 6) Q1       ~ Q2       in contact degree
                    Case 7) Q2       ~ Q3       in contact degree
                    Case 8) Q3       ~ Q4 (Max) in contact degree
                */
                int qt_start = sci_2a - 8 - 1;
                int qt_end   = sci_2a - 8;
                std::vector<int> qt_vect = find_quantile_indices(cDegList, qt_start, qt_end);
                std::shuffle(qt_vect.begin(), qt_vect.end(), gen_II);
                info_intervention_indices = std::vector<int>(qt_vect.begin(), qt_vect.begin() + num_to_select);
            }
        }
        std::unordered_set<int> iii_set(info_intervention_indices.begin(), info_intervention_indices.end());  
        // For faster inclusion checking (O(1) is faster than std::find's O(n))   set.find(target) != set.end())

        //====================================================================
        // (SIMUL INIT 2) Vector Initialization: Initial values of "Dynamic" attributes
        //====================================================================
        std::vector<double> P00_vector(N, 0.0); // P00 vector init
        std::vector<double> Ps_vector (N, 0.0); // Ps vector init   
        if (cout_flag){ std::cout << "Attempt to access " << dirs[4] << std::endl;}
        ifs.open(dirs[4]); jObj = json::parse(ifs); ifs.close(); 
        for (int i = 0; i < N; i++) { P00_vector[i] = jObj[r_idx_5][i]; }
        if (cout_flag){ std::cout << "Attempt to access " << dirs[5] << std::endl;}
        ifs.open(dirs[5]); jObj = json::parse(ifs); ifs.close(); 
        for (int i = 0; i < N; i++) { Ps_vector[i]  = jObj[r_idx_5][i]; }
        //--------------------------------------------------------------------
        // Uncertainty vectors initialization
        std::vector<double> uncertainty_P00(N, u_ind_P00); // Init individuals' uncertainty values for their P00
        std::vector<double> uncertainty_Ps(N, u_ind_Ps);   // Init individuals' uncertainty values for their Ps
        //--------------------------------------------------------------------
        // Status vectors of individual initialization  
        std::vector<bool> compliance(N, false); // Compliance status (initialized as "non-compliant (false)")
        std::vector<bool> vaccinated(N, false); // Vaccination status (initialized as "not vaccinated (false)")
        std::vector<char> dstatus(N, 'S');      // Disease status (initialized as "susceptible")
        std::vector<bool> deviation(N, false);  // Experienced Deviation    

        //====================================================================
        // (SIMUL INIT 3) Vector Initialization: temporary vectors for updates
        //====================================================================
        std::vector<double> new_probability_vector(N, 0.0);
        std::vector<double> new_uncertainty_vector(N, 0.0);
        std::vector<char>   new_dstatus(N, 'S');
        // for the "compliance" & "vaccinated" vectors, we will update directly.

        std::vector<double> vaccinated_percentage_hist(T, 0.0); 
        std::vector<double> compliance_percentage_hist(T, 0.0); 
        std::vector<std::vector<int>> epidemic_stats_hist(T, std::vector<int>(7, 0)); 
        // {|S_t|, |E_t|, |I_t|, |R_t|, # of new infections at t, vaccinated so far among 65+, # of new infections at t among 65+}


        std::vector<int> codes_age_10_percentile(N, 0);
        std::vector<int> codes_rl_10_percentile(N, 0);
        std::vector<int> codes_cdeg_10_percentile(N, 0);
        for(int idx = 0; idx < 10; idx++){
            std::vector<int> selected_indices;
            //----------------------------------------------------------------
            selected_indices = find_10_percentile_indices(age, idx);
            for (int& jdx: selected_indices){  codes_age_10_percentile[jdx] = idx; }
            //----------------------------------------------------------------
            selected_indices = find_10_percentile_indices(resource, idx);
            for (int& jdx: selected_indices){  codes_rl_10_percentile[jdx] = idx; }
            //----------------------------------------------------------------
            selected_indices = find_10_percentile_indices(cDegList, idx);
            for (int& jdx: selected_indices){  codes_cdeg_10_percentile[jdx] = idx; }
        }
        // Decile-based result (1) new infection
        std::vector<std::vector<int>> newly_exposed_per_10perc_age(T, std::vector<int>(10, 0)); 
        std::vector<std::vector<int>> newly_exposed_per_10perc_rl(T, std::vector<int>(10, 0)); 
        std::vector<std::vector<int>> newly_exposed_per_10perc_cdeg(T, std::vector<int>(10, 0)); 

        // Decile-based result (2) deviation in decision making 
        //             -- pair (# deviation while not vaccinated, # deviation given vaccinated)
        //std::vector<std::vector<int>> deviation_per_10perc_age(T, std::vector<int>(10, 0)); 
        //std::vector<std::vector<int>> deviation_per_10perc_rl(T, std::vector<int>(10, 0)); 
        //std::vector<std::vector<int>> deviation_per_10perc_cdeg(T, std::vector<int>(10, 0)); 
        std::vector<std::vector<std::pair<int, int>>> deviation_per_10perc_age(T, std::vector<std::pair<int, int>>(10, {0, 0})); 
        std::vector<std::vector<std::pair<int, int>>> deviation_per_10perc_rl(T, std::vector<std::pair<int, int>>(10, {0, 0})); 
        std::vector<std::vector<std::pair<int, int>>> deviation_per_10perc_cdeg(T, std::vector<std::pair<int, int>>(10, {0, 0})); 

        // Decile-based result (3) behavior stats. 
        std::vector<std::vector<std::pair<int, int>>> compliance_per_10perc_age(T, std::vector<std::pair<int, int>>(10, {0, 0})); // Compliance {among unvaccinated, among vaccinated} 
        std::vector<std::vector<std::pair<int, int>>> compliance_per_10perc_rl(T, std::vector<std::pair<int, int>>(10, {0, 0}));  // Compliance {among unvaccinated, among vaccinated
        std::vector<std::vector<std::pair<int, int>>> compliance_per_10perc_cdeg(T, std::vector<std::pair<int, int>>(10, {0, 0}));// Compliance {among unvaccinated, among vaccinated 
        std::vector<std::vector<int>> vaccinated_per_10perc_age(T, std::vector<int>(10, 0)); 
        std::vector<std::vector<int>> vaccinated_per_10perc_rl(T, std::vector<int>(10, 0)); 
        std::vector<std::vector<int>> vaccinated_per_10perc_cdeg(T, std::vector<int>(10, 0)); 

        //  >>>  Sizes of ( (1) Susceptible, (2) Exposed, (3) Infectious, (4) Recovered) population
        //                and (5) The number of new infection (E -> I) for each time step.

        //=======================================================================
        // (SIMUL INIT 4) Initialize the infectious agents at t = 0
        // >>> Pick the initial set of agents to have the average number of neighbors.
        //=======================================================================

        // Fixed setup -- initial sets with avg degree/income/age similar to the entire population
        std::vector<int> candidateIds;  
        std::vector<int> selected;
        for (int aidx = 0; aidx < N; aidx++) { candidateIds.push_back(aidx);}
        // 
        double t_avgDeg = 0.0, t_avgInc = 0.0, t_avgAge = 0.0;
        for (int aidx = 0; aidx < N; aidx++) { 
            t_avgDeg += cDegList[aidx];
            t_avgInc += resource[aidx];
            t_avgAge += age[aidx];
        }
        t_avgDeg = t_avgDeg/N;  t_avgInc = t_avgInc/N;  t_avgAge = t_avgAge/N;

        double e1 = 0.0, e2 = 0.0, e3 = 0.0;
        do{
            double s_avgDeg = 0.0,  s_avgInc = 0.0,  s_avgAge = 0.0;
            selected = randomSelectElements(candidateIds, I0v, gen_00); // Independent randomizer gen_00 is used.
            for (int sidx : selected){
                s_avgDeg += cDegList[sidx];
                s_avgInc += resource[sidx];
                s_avgAge += age[sidx];
            }
            s_avgDeg = s_avgDeg/I0v;  s_avgInc = s_avgInc/I0v;  s_avgAge = s_avgAge/I0v;
            e1 = std::abs(s_avgDeg-t_avgDeg)/t_avgDeg;
            e2 = std::abs(s_avgInc-t_avgInc)/t_avgInc;
            e3 = std::abs(s_avgAge-t_avgAge)/t_avgAge;
            
            if ((e1 < 0.01) && (e2 < 0.01) && (e3 < 0.01)){
                //std::cout << "[Initial infection ready]:" << e1 << " " << e2 << " " << e3 << std::endl;
                break;
            }
        } while (true);
        for (int& idx: selected){ dstatus[idx] = 'I'; }

        //====================================================================
        // Start iterations of one simulation instance (for the timespan T)
        //====================================================================

        int current_vaccine_idx_age = 0;
        int current_vaccine_idx_rl = 0;
        int current_vaccine_idx_cdeg = 0;
        json Ins0_p00log, Ins0_pslog, Ins0_u00log, Ins0_uslog; // Instance 0's logs
        json Ins0_clog, Ins0_vlog, Ins0_dlog, Ins0_dev;
        //********************************************************************
        //********************************************************************
        //***************** << RUN SIMULATION ITERATIONS >> ******************
        //********************************************************************
        //********************************************************************

        for (int current_time = 0; current_time < T; current_time++) {
            //std::cout << current_time << std::endl;
            //================================================================
            //[SETUP] Time-dependent binary indicators
            //================================================================
            // X_p: The binary variable of penalty activation for (0,0) behavior  
            //      at time t. (1: penalty activated; 0: no penalty.)
            int prev_cnet_id = current_cnet_id;
            current_cnet_id = cnet_id_sequence[current_time];
            if (prev_cnet_id != current_cnet_id){
                if (cout_flag){ std::cout << "Attempt to access contact network " << current_cnet_id << std::endl;}
                ifs.open( fmt::format("{}Data_Synthetic/[8]_ContactNetwork/ContactNet_NeiList_{}.json", inputDirectory, current_cnet_id) ); 
                jObj = json::parse(ifs);   ifs.close(); 
                for (int i = 0; i < N; i++) {  cNeiList.push_back(jObj[i]); }
            }

            if ((r_idx == 0) && (current_time%7 == 3) && (write_log_flag == 1)){
                Ins0_p00log.push_back(P00_vector); 
                Ins0_pslog.push_back(Ps_vector);  
                Ins0_u00log.push_back(uncertainty_P00); 
                Ins0_uslog.push_back(uncertainty_Ps);  
                Ins0_clog.push_back(compliance);   
                Ins0_vlog.push_back(vaccinated);   
                Ins0_dlog.push_back(dstatus);
                Ins0_dev.push_back(deviation);
            }
            //
            //================================================================
            //[SETUP] Time-dependent binary indicators
            //================================================================
            // X_p: The binary variable of penalty activation for (0,0) behavior  
            //      at time t. (1: penalty activated; 0: no penalty.)
            int X_p  = 0; 
            if (current_time >= t_1 && current_time < t_2) { X_p = 1; }  // C_n is applied to all non-compliant agents
            else if (current_time >= t_2) { X_p = 2; }           // C_n is applied to non-compliant agents who are not vaccinated.
            //----------------------------------------------------------------
            if (current_time > vaccine_update_times_age [current_vaccine_idx_age]) { current_vaccine_idx_age++; }
            if (current_time > vaccine_update_times_rl  [current_vaccine_idx_rl])  { current_vaccine_idx_rl++; }
            if (current_time > vaccine_update_times_cdeg[current_vaccine_idx_cdeg]){ current_vaccine_idx_cdeg++; }

            //================================================================
            //[PHASE #1] PROBABILITY/UNCERTAINTY UPDATING PHASE
            //================================================================    
            double rndMedia;    //double rndComp;    double rndChng; 
            int nidx;   double h_ij;
            double* cown = NULL;
            double* oown = NULL, * uown = NULL, * onei = NULL, * unei = NULL;
            bool updateFlag = false;
            //================================================================    
            // - Before vaccine (time < t_3)  
            // >>>> Update P00 
            // - Vaccination phase 1 & 2 (time >= t_3) 
            // >>>> Update Ps
            //================================================================
            //[BEGIN OPINION UPDATE PROCEDURE]================================
            // Randomizer for opinion dynamics is used in this section: gen_OD
            for (int aidx = 0; aidx < N; aidx++) {
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                updateFlag = false;
                rndMedia = rand01(gen_OD);
                cown = &chan1Usg[aidx];
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                double effective_learn_P00 = learn_P00;
                double effective_learn_u00 = learn_u00;
                double effective_learn_Ps  = learn_Ps;
                double effective_learn_us  = learn_us;
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                bool from_media = (rndMedia < mediaUsg[aidx]); // Determine info source between media & neighbor
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                // Setup values for opinion update:
                if (current_time < t_3){  oown = &P00_vector[aidx];  uown = &uncertainty_P00[aidx];  } 
                else                   {  oown = &Ps_vector[aidx];   uown = &uncertainty_Ps[aidx];   }

                if (from_media) { // Info source is the media channels.
                    // [ Agent updates attributes using mass media's message ]
                    ddist.param({ *cown, 1 - *cown });
                    nidx = ddist(gen_OD);
                    if (current_time < t_3){ onei = &c_P00[nidx];  unei = &c_u00[nidx]; }
                    else                   { onei = &c_Ps[nidx];   unei = &c_us[nidx];  }
                    h_ij = std::min(*oown + *uown, *onei + *unei) - std::max(*oown - *uown, *onei - *unei);
                    updateFlag = true;
                }
                else {
                    // [ Agent updates attributes using neighbors' message ]
                    std::vector<int> activeNei = oNeiList[aidx];
                    if ( (activeNei.size() != 0) ) {
                        std::vector<double> weights(activeNei.size(), 1.0 / activeNei.size());
                        ddist.param({ begin(weights), end(weights) });
                        nidx = activeNei[ddist(gen_OD)];
                        //----------------------------------------------------
                        if (current_time < t_3){ onei = &P00_vector[nidx];  unei = &uncertainty_P00[nidx]; }
                        else                   { onei = &Ps_vector[nidx];   unei = &uncertainty_Ps[nidx];  }
                        h_ij = std::min(*oown + *uown, *onei + *unei) - std::max(*oown - *uown, *onei - *unei);
                        updateFlag = true;
                    }
                }
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                // < Reactive Devaluation -- Learning Rate Updates >
                if( !info_intervention_indices.empty() && from_media && iii_set.find(aidx) != iii_set.end()){
                    if (nidx == 0){ // When pro-intervention channel (channel 1 with idx 0) is selected
                        effective_learn_P00 = 0.1 * learn_P00;
                        effective_learn_u00 = 0.1 * learn_u00;
                        effective_learn_Ps  = 0.1 * learn_Ps;
                        effective_learn_us  = 0.1 * learn_us;
                    } else if (nidx == 1){ // When anti-intervention channel (channel 2 with idx 1) is selected
                        effective_learn_P00 = 2.0 * learn_P00;
                        effective_learn_u00 = 2.0 * learn_u00;
                        effective_learn_Ps  = 2.0 * learn_Ps;
                        effective_learn_us  = 2.0 * learn_us;
                        // Update the uncertainty too. 
                        if (current_time < t_3){ unei = &c_u00[2]; } else { unei = &c_us[2]; }
                        h_ij = std::min(*oown + *uown, *onei + *unei) - std::max(*oown - *uown, *onei - *unei);
                    } else {
                        std::cout << "Error: learning rate adjustment." << std::endl;
                    }
                } 

                // << Update attribute values >>
                if (current_time < t_3){  
                    new_probability_vector[aidx] = P00_vector[aidx];
                    new_uncertainty_vector[aidx] = uncertainty_P00[aidx];
                    if ( (updateFlag == true) && (h_ij > *unei) && (*unei < 4.0) ) {
                        new_probability_vector[aidx] += effective_learn_P00 * ((h_ij / *unei) - 1) * (*onei - *oown);
                        new_uncertainty_vector[aidx] += effective_learn_u00 * ((h_ij / *unei) - 1) * (*unei - *uown);
                    }
                }
                else{  
                    new_probability_vector[aidx] = Ps_vector[aidx];
                    new_uncertainty_vector[aidx] = uncertainty_Ps[aidx];
                    if ( (updateFlag == true) && (h_ij > *unei) && (*unei < 4.0) ) {
                        new_probability_vector[aidx] += effective_learn_Ps * ((h_ij / *unei) - 1) * (*onei - *oown);
                        new_uncertainty_vector[aidx] += effective_learn_us * ((h_ij / *unei) - 1) * (*unei - *uown);
                    }
                }
            }
            //[Transfer from temp storage]------------------------------------
            if (current_time < t_3){  
                P00_vector      = new_probability_vector;
                uncertainty_P00 = new_uncertainty_vector;
            }
            else{
                Ps_vector      = new_probability_vector;
                uncertainty_Ps = new_uncertainty_vector;
            }
            //[END OPINION UPDATE PROCEDURE]==================================
            //
            //================================================================
            //[PHASE #2] VACCINATION PHASE
            //================================================================    
            // (1) Vaccine Administration Capacity Setting:
            //    -- Depending on the proportion assigned
            int capacity_1 = 0;
            int capacity_2 = 0;
            if (current_time >= t_4) { 
                // After t_4 (The date on which the vaccination phase 2 started.) 
                // S_2 (X4) Daily vaccine administration capacity during phase 2
                capacity_1 = std::round(S_2 * sci_3b); 
                capacity_2 = S_2 - capacity_1;
            } else if (current_time >= t_3) { 
                // After t_3 (The date on which the vaccination phase 1 started.) but before t_4
                // S_1 (X3) Daily vaccine administration capacity during phase 1
                capacity_1 = std::round(S_1 * sci_3b); 
                capacity_2 = S_1 - capacity_1;
            } 
            // -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - ---
            // (2) Threshold value update:
            int age_bar = vaccine_thresholds_age [current_vaccine_idx_age];
            int rl_bar  = vaccine_thresholds_rl  [current_vaccine_idx_rl];
            int cdeg_bar= vaccine_thresholds_cdeg[current_vaccine_idx_cdeg];
            // -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - ---
            // (3-1) First part of the vaccine administration (age-based eligibility)
            std::vector<int> willing_vect_qTrue;  // A set of IDs representing agents willing to get vaccinated despite the queue congestion cost.
            std::vector<int> willing_vect_qFalse; // A set of IDs representing agents willing to get vaccinated only in the absence of congestion.
            if (capacity_1 > 0){
                for (int idx = 0; idx < N; idx++) {
                    // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    if (age[idx] < age_bar){ continue; }       // Pass -- If not eligible
                    if (vaccinated[idx] == true){ continue; }  // Pass -- If already vaccinated.
                    if (dstatus[idx] == 'I'){ continue; }
                    // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    int fill_1 = 1; int fill_0 = 0;
                    bool w_true  = v_willing(age[idx], resource[idx], P00_vector[idx], Ps_vector[idx], X_p, fill_1); // With congestion -- value "1"
                    bool w_false = v_willing(age[idx], resource[idx], P00_vector[idx], Ps_vector[idx], X_p, fill_0); // W\O  congestion -- value "0"
                    // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // - Collect agent ids who are willing to get vaccinated assuming queue congestion
                    if (w_true  == true){ willing_vect_qTrue.push_back(idx);}
                    // - Collect agent ids who cannot tolerate the queue congestion, but otherwise willing to get vaccinated.
                    else if (w_false == true){ willing_vect_qFalse.push_back(idx);}
                }
                // Shuffle for randomness if needed.
                int qTrue_size = static_cast<int> (willing_vect_qTrue.size());
                int qFalse_size = static_cast<int> (willing_vect_qFalse.size());
                if ( capacity_1 <= qTrue_size ) {  
                    // Capacity cannot cover all qTrue;
                    std::shuffle(std::begin(willing_vect_qTrue),  std::end(willing_vect_qTrue),  gen_VM);
                    for (int idx = 0; idx < capacity_1; idx++){ 
                        vaccinated[willing_vect_qTrue[idx]] = true; 
                        dstatus[willing_vect_qTrue[idx]] = 'R';
                    }                        
                } else if (capacity_1 < (qTrue_size + qFalse_size)) { 
                    // Capacity can cover qTrue, but not all qFalse.
                    std::shuffle(std::begin(willing_vect_qFalse), std::end(willing_vect_qFalse), gen_VM);
                    for (int& idx: willing_vect_qTrue){ vaccinated[idx] = true; dstatus[idx] = 'R';}
                    for (int idx = 0; idx < (capacity_1 - qTrue_size); idx++){ 
                        vaccinated[willing_vect_qFalse[idx]] = true;  
                        dstatus[willing_vect_qFalse[idx]] = 'R';
                    }
                } else {  
                    // Capacity can cover qTrue, and also all qFalse.
                    for (int& idx: willing_vect_qTrue){  vaccinated[idx] = true; dstatus[idx] = 'R';}
                    for (int& idx: willing_vect_qFalse){ vaccinated[idx] = true; dstatus[idx] = 'R';}
                }
            }
            // -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - ---
            // (3-2) Second part of the vaccine administration (secondary eligibility rule)
            willing_vect_qTrue.clear(); 
            willing_vect_qFalse.clear();
            if (sci_3a != 0 && capacity_2 != 0){
                std::vector<bool> eligibility_flags(N, false);
                switch (sci_3a) { 
                    case 1: {// Open Eligibility
                        for (int idx = 0; idx < N; idx++) { eligibility_flags[idx] = true; } 
                        break;
                    } 
                    case 2: {// Exposure-based Eligibility {3 - 2 - 1 - 0, LB};
                        for (int idx = 0; idx < N; idx++) { 
                            if (cDegList[idx] < cdeg_bar){ continue; } // Pass -- If not eligible or If already vaccinated.
                            eligibility_flags[idx] = true; 
                        } 
                        break;
                    } 
                    case 3: {// Resorce-availability-based Eligibility {UB} 
                        for (int idx = 0; idx < N; idx++) { 
                            if (resource[idx] > rl_bar){ continue; } // Pass -- If not eligible or If already vaccinated.
                            eligibility_flags[idx] = true; 
                        } 
                        break;
                    } 
                    default:
                        break; 
                }
                for (int idx = 0; idx < N; idx++) {
                    // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    if (vaccinated[idx] == true){ continue; }  // Pass -- If already vaccinated.
                    if (eligibility_flags[idx] == false){ continue; }
                    if (dstatus[idx] == 'I'){ continue; }
                    // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    int fill_1 = 1; int fill_0 = 0;
                    bool w_true  = v_willing(age[idx], resource[idx], P00_vector[idx], Ps_vector[idx], X_p, fill_1); // With congestion -- value "1"
                    bool w_false = v_willing(age[idx], resource[idx], P00_vector[idx], Ps_vector[idx], X_p, fill_0); // W\O  congestion -- value "0"
                    // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // - Collect agent ids who are willing to get vaccinated assuming queue congestion
                    if (w_true  == true){ willing_vect_qTrue.push_back(idx);}
                    // - Collect agent ids who cannot tolerate the queue congestion, but otherwise willing to get vaccinated.
                    else if (w_false == true){ willing_vect_qFalse.push_back(idx);}
                }
                // Shuffle for randomness if needed.
                int qTrue_size = static_cast<int> (willing_vect_qTrue.size());
                int qFalse_size = static_cast<int> (willing_vect_qFalse.size());
                if ( capacity_2 <= qTrue_size ) {  
                    // Capacity cannot cover all qTrue;
                    std::shuffle(std::begin(willing_vect_qTrue),  std::end(willing_vect_qTrue),  gen_VM);
                    for (int idx = 0; idx < capacity_2; idx++){ 
                        vaccinated[willing_vect_qTrue[idx]] = true; 
                        dstatus[willing_vect_qTrue[idx]] = 'R';                        
                    }                        
                } else if (capacity_2 < (qTrue_size + qFalse_size)) { 
                    // Capacity can cover qTrue, but not all qFalse.
                    std::shuffle(std::begin(willing_vect_qFalse), std::end(willing_vect_qFalse), gen_VM);
                    for (int& idx: willing_vect_qTrue){ vaccinated[idx] = true;  dstatus[idx] = 'R';}
                    for (int idx = 0; idx < (capacity_2 - qTrue_size); idx++){ 
                        vaccinated[willing_vect_qFalse[idx]] = true;  
                        dstatus[willing_vect_qFalse[idx]] = 'R';
                    }
                } else {  
                    // Capacity can cover qTrue, and also all qFalse.
                    for (int& idx: willing_vect_qTrue){  vaccinated[idx] = true; dstatus[idx] = 'R';}
                    for (int& idx: willing_vect_qFalse){ vaccinated[idx] = true; dstatus[idx] = 'R';}
                }
            }

            //
            //================================================================
            //[PHASE #3] NPI COMPLIANCE PHASE
            //================================================================    
            int X_0 = 0;   
            bool null_decision = false;
            for (int idx = 0; idx < N; idx++) {
                // 
                int code_age  = codes_age_10_percentile[idx];
                int code_rl   = codes_rl_10_percentile[idx];
                int code_cdeg = codes_cdeg_10_percentile[idx];                
                // 
                if (vaccinated[idx] == false){
                    // Not vaccinated yet >> B0 decision function
                    compliance[idx] = decision_B0(age[idx], resource[idx], P00_vector[idx], X_p);
                    null_decision   = decision_B0(age[idx], resource[idx], P00_vector[idx], X_0);
                    deviation[idx]  = (compliance[idx] != null_decision);
                    if(deviation[idx]){
                        deviation_per_10perc_age[current_time][code_age].first++; 
                        deviation_per_10perc_rl[current_time][code_rl].first++; 
                        deviation_per_10perc_cdeg[current_time][code_cdeg].first++;
                    }
                }else{
                    // Already vaccinated >> B1 decision function
                    compliance[idx] = decision_B1(age[idx], resource[idx], X_p);
                    null_decision   = decision_B1(age[idx], resource[idx], X_0);
                    deviation[idx]  = (compliance[idx] != null_decision);
                    if(deviation[idx]){
                        deviation_per_10perc_age[current_time][code_age].second++; 
                        deviation_per_10perc_rl[current_time][code_rl].second++; 
                        deviation_per_10perc_cdeg[current_time][code_cdeg].second++;
                    }
                }
                /*
                if(deviation[idx]){
                    deviation_per_10perc_age[current_time][code_age]++; 
                    deviation_per_10perc_rl[current_time][code_rl]++; 
                    deviation_per_10perc_cdeg[current_time][code_cdeg]++;
                }
                */
            }

            //================================================================
            // Gather Additional Behavior Stats
            for (int idx = 0; idx < N; idx++) {
                int code_age  = codes_age_10_percentile[idx];
                int code_rl   = codes_rl_10_percentile[idx];
                int code_cdeg = codes_cdeg_10_percentile[idx];    
                if (compliance[idx] == true){ 
                    if (vaccinated[idx] == false){ // Not vaccinated compliance
                        compliance_per_10perc_age[current_time][code_age].first++; 
                        compliance_per_10perc_rl[current_time][code_rl].first++; 
                        compliance_per_10perc_cdeg[current_time][code_cdeg].first++; 
                    } else {
                        compliance_per_10perc_age[current_time][code_age].second++; 
                        compliance_per_10perc_rl[current_time][code_rl].second++; 
                        compliance_per_10perc_cdeg[current_time][code_cdeg].second++; 
                    }
                }
                if (vaccinated[idx] == true){ 
                    vaccinated_per_10perc_age[current_time][code_age]++; 
                    vaccinated_per_10perc_rl[current_time][code_rl]++; 
                    vaccinated_per_10perc_cdeg[current_time][code_cdeg]++; 
                }
            }

            //================================================================
            //[PHASE #4] DISEASE STATUS UPDATING PHASE
            //================================================================    
            double current_EFF = ((current_time < ChT0) || (current_time > ChT1)) ? EFF : (EFF *(1.0 - ChR));
            if (current_EFF > 1.0){ current_EFF = 1.0; }
            //-------------------------------------------------------------------
            double effectiveBeta;
            double cRnd;
            for (int aidx = 0; aidx < N; aidx++) {
                if (dstatus[aidx] == 'S') {
                    for (int nidx : cNeiList[aidx]) {
                        if (dstatus[nidx] == 'I') {
                            cRnd = rand01(gen_DD);
                            effectiveBeta = (compliance[aidx] || compliance[nidx]) ? trans_val * (1-current_EFF) : trans_val; 
                            if (cRnd < effectiveBeta) {
                                new_dstatus[aidx] = 'E';
                                break;
                            }
                        }
                    }
                }
                else {
                    cRnd = rand01(gen_DD);
                    if (dstatus[aidx] == 'E') {
                        if (cRnd < sigma_val) {
                            new_dstatus[aidx] = 'I';
                            epidemic_stats_hist[current_time][4]++; // Count new infection
                            if ( age[aidx] >= 65 ) { epidemic_stats_hist[current_time][6]++; } // Count new infection among 65+
                            //---------------------------------------------------
                            // Gather 10-percentile chunk based new infection stats
                            int code_age  = codes_age_10_percentile[aidx];
                            int code_rl   = codes_rl_10_percentile[aidx];
                            int code_cdeg = codes_cdeg_10_percentile[aidx];
                            newly_exposed_per_10perc_age[current_time][code_age]++;
                            newly_exposed_per_10perc_rl[current_time][code_rl]++;
                            newly_exposed_per_10perc_cdeg[current_time][code_cdeg]++;
                            //---------------------------------------------------
                        }
                        else {
                            new_dstatus[aidx] = 'E';
                        }
                    }
                    else if (dstatus[aidx] == 'I'){ new_dstatus[aidx] = (cRnd < gamma_val) ? 'R' : 'I';}
                    else if (dstatus[aidx] == 'R'){ new_dstatus[aidx] = (cRnd < omega_val) ? 'S' : 'R';}
                }
            }

            //===================================================================
            // (SIMUL PART 3-2) Save log of disease dynamics subprocess
            //===================================================================
            epidemic_stats_hist[current_time][0] = countOccurrences(dstatus, 'S');
            epidemic_stats_hist[current_time][1] = countOccurrences(dstatus, 'E');
            epidemic_stats_hist[current_time][2] = countOccurrences(dstatus, 'I');
            epidemic_stats_hist[current_time][3] = countOccurrences(dstatus, 'R');

            for (int aidx = 0; aidx < N; aidx++) { if (vaccinated[aidx] == true && age[aidx] >= 65) { epidemic_stats_hist[current_time][5]++; } }
            //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            dstatus = new_dstatus;
            //
            //================================================================
            //[RESULT #?] LOG SAVING PHASE
            //================================================================    
            int c_count = std::count(compliance.begin(), compliance.end(), true);  //for (bool com: compliance){ if(com == true) c_count++; }
            compliance_percentage_hist[current_time] = std::round( 10000.0* (1.0 * c_count) / N) / 10000.0; 
            //----------------------------------------------------------------
            int v_count = std::count(vaccinated.begin(), vaccinated.end(), true);  //for (bool vac: vaccinated){ if(vac == true) v_count++; }
            vaccinated_percentage_hist[current_time] = std::round( 10000.0* (1.0 * v_count) / N) / 10000.0; 
        }


        if (r_idx == 0 && write_log_flag == 1){
            //----------------------------------------------------------------
            ofs.open(outputDirectory+"/Ins0_p00_log_" + argStr + ".json"); 
            ofs << Ins0_p00log; ofs.close();    
            //----------------------------------------------------------------
            ofs.open(outputDirectory+"/Ins0_ps_log_" + argStr + ".json"); 
            ofs << Ins0_pslog; ofs.close();    
            //----------------------------------------------------------------
            ofs.open(outputDirectory+"/Ins0_u00log_" + argStr + ".json"); 
            ofs << Ins0_u00log; ofs.close();    
            //----------------------------------------------------------------
            ofs.open(outputDirectory+"/Ins0_uslog_" + argStr + ".json"); 
            ofs << Ins0_uslog; ofs.close();    
            //----------------------------------------------------------------
            ofs.open(outputDirectory+"/Ins0_clog_" + argStr + ".json"); 
            ofs << Ins0_clog; ofs.close();    
            //----------------------------------------------------------------
            ofs.open(outputDirectory+"/Ins0_vlog_" + argStr + ".json"); 
            ofs << Ins0_vlog; ofs.close();    
            //----------------------------------------------------------------
            ofs.open(outputDirectory+"/Ins0_dlog_" + argStr + ".json"); 
            ofs << Ins0_dlog; ofs.close();    
            //----------------------------------------------------------------
            ofs.open(outputDirectory+"/Ins0_dev_" + argStr + ".json"); 
            ofs << Ins0_dev; ofs.close();    
        }
        //for (double val: vaccinated_percentage_hist){ std::cout << val << " ";} 
        if (write_log_flag < 2){
            JO_clog.push_back( compliance_percentage_hist ); // Compliance percentage logs of instances
            JO_vlog.push_back( vaccinated_percentage_hist ); // Vaccination percentage logs of instances
            JO_dlog.push_back( epidemic_stats_hist ); 
            //----------------------------------------------------------------
            JO_ne10_age.push_back(newly_exposed_per_10perc_age);  // for "newly_exposed_per_10perc_age"
            JO_ne10_rl.push_back(newly_exposed_per_10perc_rl);   // for "newly_exposed_per_10perc_rl"
            JO_ne10_cdeg.push_back(newly_exposed_per_10perc_cdeg); // for "newly_exposed_per_10perc_cdeg"
            //----------------------------------------------------------------
            JO_dev10_age.push_back(deviation_per_10perc_age); // for "deviation_per_10perc_age"
            JO_dev10_rl.push_back(deviation_per_10perc_rl);  // for "deviation_per_10perc_rl"
            JO_dev10_cdeg.push_back(deviation_per_10perc_cdeg);// for "deviation_per_10perc_cdeg"
            //----------------------------------------------------------------
            JO_com10_age.push_back(compliance_per_10perc_age);
            JO_com10_rl.push_back(compliance_per_10perc_rl);
            JO_com10_cdeg.push_back(compliance_per_10perc_cdeg);
            //----------------------------------------------------------------
            JO_vac10_age.push_back(vaccinated_per_10perc_age);
            JO_vac10_rl.push_back(vaccinated_per_10perc_rl);
            JO_vac10_cdeg.push_back(vaccinated_per_10perc_cdeg);
        }
        //====================================================================
        double score_current = 0.0;
        for (int tidx = 0; tidx < 64; tidx++){ // 79 >> 65 (till day 450)
            //std::cout << tidx << std::endl;
            int tval = cal_time[tidx]; 
            //----------------------------------------------------------------
            if (cal_obsd[tidx] > 0){ 
                double tmp_val = 0.0;
                tmp_val = cal_obsd[tidx] - epidemic_stats_hist[tval][4]; 
                score_current += (tmp_val * tmp_val);      
            }            
            //----------------------------------------------------------------
        }
        score += score_current; 
    }
    // "score" -- the cumulative score over n_reps.
    //===========================================================================
    // [SIMUL END] Iteratively run the requested number of simulation instances
    //===========================================================================
    if (write_log_flag < 2){
        ofs.open(outputDirectory+"/clog_" + argStr + ".json"); 
        ofs << JO_clog; ofs.close();    
        ofs.open(outputDirectory+"/vlog_" + argStr + ".json"); 
        ofs << JO_vlog; ofs.close();    
        ofs.open(outputDirectory+"/dlog_" + argStr + ".json"); 
        ofs << JO_dlog; ofs.close();    
        //===========================================================================
        ofs.open(outputDirectory+"/com10_age_" + argStr + ".json"); 
        ofs << JO_com10_age; ofs.close();    
        ofs.open(outputDirectory+"/com10_rl_" + argStr + ".json"); 
        ofs << JO_com10_rl; ofs.close();    
        ofs.open(outputDirectory+"/com10_cdeg_" + argStr + ".json"); 
        ofs << JO_com10_cdeg; ofs.close();    
        //===========================================================================
        ofs.open(outputDirectory+"/vac10_age_" + argStr + ".json"); 
        ofs << JO_vac10_age; ofs.close();    
        ofs.open(outputDirectory+"/vac10_rl_" + argStr + ".json"); 
        ofs << JO_vac10_rl; ofs.close();    
        ofs.open(outputDirectory+"/vac10_cdeg_" + argStr + ".json"); 
        ofs << JO_vac10_cdeg; ofs.close();    
        //===========================================================================
        ofs.open(outputDirectory+"/ne10_age_" + argStr + ".json"); 
        ofs << JO_ne10_age; ofs.close();    
        ofs.open(outputDirectory+"/ne10_rl_" + argStr + ".json"); 
        ofs << JO_ne10_rl; ofs.close();    
        ofs.open(outputDirectory+"/ne10_cdeg_" + argStr + ".json"); 
        ofs << JO_ne10_cdeg; ofs.close();    
        //===========================================================================
        ofs.open(outputDirectory+"/dev10_age_" + argStr + ".json"); 
        ofs << JO_dev10_age; ofs.close();    
        ofs.open(outputDirectory+"/dev10_rl_" + argStr + ".json"); 
        ofs << JO_dev10_rl; ofs.close();    
        ofs.open(outputDirectory+"/dev10_cdeg_" + argStr + ".json"); 
        ofs << JO_dev10_cdeg; ofs.close();    
    }
    std::cout << score << std::endl;
    return 0;
}


