#ifndef GLOBALS_2_H
#define GLOBALS_2_H

#include <vector>

// ===========================================================================
// <<  HARD-CODED CALIBRATION SEQUENCES & VACCINE ELIGIBILITY SCHEDULES  >>
// ===========================================================================
extern std::vector<int>    cal_time;
extern std::vector<double> cal_obsb, cal_obsv, cal_obsd;
extern std::vector<int> vaccine_update_times1;
extern std::vector<int> vaccine_thresholds_age;
extern std::vector<int> vaccine_update_times2;
extern std::vector<int> vaccine_thresholds_qts;


extern std::vector<int> vaccine_update_times_age;
extern std::vector<int> vaccine_thresholds_age;

extern std::vector<int> vaccine_update_times_rl;
extern std::vector<int> vaccine_thresholds_rl;

extern std::vector<int> vaccine_update_times_cdeg;
extern std::vector<int> vaccine_thresholds_cdeg;

#endif
