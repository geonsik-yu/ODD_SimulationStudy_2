// ===========================================================================
// <<  FUNCTION DEFINITIONS (PART 2: SIMULATION CONTEXT FUNCTIONS)  >>
// ===========================================================================

#include <vector>     // For std::vector
#include <cmath>      // For std::pow, std::abs
#include <algorithm>  // For std::max_element, std::max
#include <iostream>   // for std::cout
#include "include/globals_1.h" //#include "include/global_2.h"
#include "include/simul-logic.h"

// ---------------------------------------------------------------------------
// 1) Utility function
double uf(double& resource) {
    double result = 0.0;
    double tmp_diff = UPOWV - 1.00;
    if (abs(tmp_diff) < 0.01){
        result = resource;
    } else {
        std::vector<double> v = {C_p, C_q, C_q + C_n};
        double bound = *std::max_element(std::begin(v), std::end(v));
        double tmp = bound + resource;
        result = std::pow(tmp, UPOWV);
    }
    return result;
}

// ---------------------------------------------------------------------------
// 2) Decision function -- Vaccination willingness // not vaccinated yet
bool v_willing(double& age, double& Y, double& p00, double& ps, int& X_p, int& X_q) {
    // Inputs:  age (age), Y (resource level), 
    //          p00 (perceived p00 - infectious with 00), 
    //          ps  (perceived ps - side effect), 
    //          X_p (NPI related policy flag),
    //              >>   X_p == 0:   C_p is not applied to anyone.              (  0  ~ t_1 )
    //              >>   X_p == 1:   C_p is applied to all non-compliant agents ( t_1 ~ t_2 )
    //              >>   X_p == 2:   C_p is applied to all non-compliant agents ( t_2 ~ end )
    //          X_q (Queue congestion flag).
    double alpha = ((age < 65) ? alpha_le64 : alpha_ge65);  // Resource modifier (infectious)
    double xi = alpha;                                      // Resource modifier (side effect)
    double eu00, eu10, eu01, eu11;  // Utility storage variables
    double tmp1, tmp2, tmp3;        // Temporary storage variables for computation
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Social penalty costs of non-compliance with NPI
    double effective_C_p_1; // Effective penalty for 00: no vaccination + no NPI compliance
    double effective_C_p_2; // Effective penalty for 01: get vaccinated + no NPI compliance
    if( X_p == 0 ) {  effective_C_p_1 = 0;    effective_C_p_2 = 0;   }
    if( X_p == 1 ) {  effective_C_p_1 = C_p;  effective_C_p_2 = C_p; }  
    if( X_p == 2 ) {  effective_C_p_1 = C_p;  effective_C_p_2 = 0;   }  
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Behavior Code (0,0): 
    // >> Not vaccinated + No NPI compliance
    tmp1 = (alpha*Y - effective_C_p_1 - AbsoluteLoss); // Case of infection    || Social penalty if X_p = True 
    tmp2 = (      Y - effective_C_p_1               ); // Case of no infection || Social penalty if X_p = True
    eu00 =    p00     * uf(tmp1) 
         + (1-p00)    * uf(tmp2);
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Behavior Code (1,0): 
    // >> Not vaccinated + complying with NPI recommendations
    tmp1 = (alpha*Y - C_n - AbsoluteLoss); // Case of infection    || Cost of compliance 
    tmp2 = (      Y - C_n               ); // Case of no infection || Cost of compliance
    eu10 =    P10     * uf(tmp1)     
         + (1-P10)    * uf(tmp2);
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Behavior Code (0,1): 
    // >> Willing to get vaccinated + No NPI compliance
    tmp1 = (alpha*Y - effective_C_p_2 - X_q * C_q - AbsoluteLoss); // Case of infection    || Social penalty if X_p = True || queue cost
    tmp2 = (   xi*Y - effective_C_p_2 - X_q * C_q - AbsoluteLoss); // Case of side effect  || Social penalty if X_p = True || queue cost 
    tmp3 = (      Y - effective_C_p_2 - X_q * C_q               ); // Case of no infection || Social penalty if X_p = True || queue cost 
    eu01 =    P01     * uf(tmp1)
         +    ps      * uf(tmp2)
         + (1-P01-ps) * uf(tmp3);
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Behavior Code (1,1): 
    // >> Willing to get vaccinated + complying with NPI recommendations
    tmp1 = (alpha*Y - C_n - X_q * C_q - AbsoluteLoss); // Case of infection    || Cost of compliance || queue cost
    tmp2 = (   xi*Y - C_n - X_q * C_q - AbsoluteLoss); // Case of side effect  || Cost of compliance || queue cost 
    tmp3 = (      Y - C_n - X_q * C_q               ); // Case of no infection || Cost of compliance || queue cost 
    eu11 =    P11     * uf(tmp1)
         +    ps      * uf(tmp2)       
         + (1-P11-ps) * uf(tmp3);
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // returns "TRUE" if the agent is willing to get vaccinated.
    // returns "FALSE" otherwise.
    bool willingness = false;
    if (std::max(eu01, eu11) >= std::max(eu00, eu10)){ willingness = true; }
    return willingness;
}

// ---------------------------------------------------------------------------
// 3) Decision function -- Compliance decision (not vaccinated yet)
bool decision_B0(double& age, double& Y, double& p00, int& X_p) {
    // Inputs:  age (age), Y (resource level), 
    //          p00 (perceived p00 - infectious with 00), 
    //          ps  (perceived ps - side effect), 
    //          X_p (NPI related policy flag),
    //              >>   X_p == 0:   C_n is not applied to anyone.              (  0  ~ t_1 )
    //              >>   X_p == 1:   C_n is applied to all non-compliant agents ( t_1 ~ t_2 )
    //              >>   X_p == 2:   C_n is applied to all non-compliant agents ( t_2 ~ end )
    //          X_q (Queue congestion flag).

    double eu00, eu10;
    double alpha = ((age < 65) ? alpha_le64 : alpha_ge65);
    double tmp1, tmp2;
    double effective_C_p = C_p;
    if( X_p == 0 ) { effective_C_p = 0; }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Expected utility of do nothing
    tmp1 = (alpha*Y - effective_C_p - AbsoluteLoss); // Case of infection    || Social penalty if X_p = True
    tmp2 = (      Y - effective_C_p               ); // Case of no infection || Social penalty if X_p = True
    eu00 =     p00  * uf(tmp1) 
          + (1-p00) * uf(tmp2);
    // Expected utility of NPI compliance
    tmp1 = (alpha*Y - C_n - AbsoluteLoss); // Case of infection    || Cost of compliance
    tmp2 = (      Y - C_n               ); // Case of no infection || Cost of compliance
    eu10 =     P10  * uf(tmp1)                 
          + (1-P10) * uf(tmp2); 
    bool compliance = ((eu00 >= eu10) ? false : true);
    return compliance; 
    // returns 0 for not complying with NPIs
    // returns 1 for otherwise

}

// ---------------------------------------------------------------------------
// 4) Decision function -- Compliance decision (vaccinated already)
bool decision_B1(double& age, double& Y, int& X_p) {
    // Inputs:  age (age), Y (resource level), 
    //          p00 (perceived p00 - infectious with 00), 
    //          ps  (perceived ps - side effect), 
    //          X_p (NPI related policy flag),
    //              >>   X_p == 0:   C_n is not applied to anyone.              (  0  ~ t_1 )
    //              >>   X_p == 1:   C_n is applied to all non-compliant agents ( t_1 ~ t_2 )
    //              >>   X_p == 2:   C_n is applied to all non-compliant agents ( t_2 ~ end )
    //          X_q (Queue congestion flag).

    double eu01, eu11;
    double alpha = ((age < 65) ? alpha_le64 : alpha_ge65);
    double tmp1, tmp2;
    double effective_C_p = C_p;
    if( X_p == 0 || X_p == 2 ) { effective_C_p = 0; }
    //------------------------------------------------------------------------
    // Expected utility of do nothing
    tmp1 = (alpha*Y - effective_C_p - AbsoluteLoss); // Case of infection    || Social penalty if X_p = True
    tmp2 = (      Y - effective_C_p               ); // Case of no infection || Social penalty if X_p = True
    eu01 =     P01  * uf(tmp1)          
          + (1-P01) * uf(tmp2);
    // Expected utility of NPI compliance
    tmp1 = (alpha*Y - C_n - AbsoluteLoss); // Case of infection    || Cost of compliance
    tmp2 = (      Y - C_n               ); // Case of no infection || Cost of compliance
    eu11 =     P11  * uf(tmp1) 
          + (1-P11) * uf(tmp2);
    bool compliance = ((eu01 >= eu11) ? false : true);
    return compliance; 
    // returns 0 for not complying with NPIs
    // returns 1 for otherwise
}

