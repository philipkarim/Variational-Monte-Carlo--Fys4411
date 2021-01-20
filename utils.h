#include <iostream>


extern double h_bar, a, beta, alpha, a_ho, r;
extern double m, w_ho, w_z;
extern int interaction;
extern char trap;


double local_energy(int x, int y);
double V_external(double x, double y, double z,char trap);
double V_internal(double r_distance);
double distance_particle(double x1,double y1,double z1, double x2, double y2, double z2);
double g_function(double x, double y, double z);
double f_function(double a, double r_dist);
void declare_trap(char ho_trap, int interacting);
