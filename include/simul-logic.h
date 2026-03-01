#ifndef SIMULLOGIC_H
#define SIMULLOGIC_H

// 1) Utility function
double uf(double& resource);

// 2) Decision function -- Vaccination willingness // not vaccinated yet
bool v_willing(double& age, double& Y, double& p00, double& ps, int& X_p, int& X_q);

// 3) Decision function -- Compliance decision (not vaccinated yet)
bool decision_B0(double& age, double& Y, double& p00, int& X_p);

// 4) Decision function -- Compliance decision (vaccinated already)
bool decision_B1(double& age, double& Y, int& X_p);

#endif
