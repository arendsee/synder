#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

struct {
    // Likelihood of each event (should sum to 1)
    double flip_rate;
    double copy_rate;
    double move_rate;
    double delete_rate;
    double repeat_rate;

    // Length of affected region
    double flip_width;
    double copy_width;
    double move_width;
    double delete_width;
    double repeat_width;
} SimulationParameters;

struct S_Parameters default_parameters();

struct S_Parameters init_s_parameters();

#endif
