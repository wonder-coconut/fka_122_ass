#pragma once

double inst_temp_calc(double, double);
double inst_press_calc(double, double, double);
double alpha_T_calc(double, double, double, double);
double alpha_P_calc(double , double, double, double);
void eqb_temp_1_step(double **, double , double, double, double, double, int);
void eqb_press_1_step(double **, double, double, double, double, double, double, int, double*);