#include <iostream>
#include "utils.h"

using namespace std;

//#define Globalvariables søk på det
//Constants, most of them is set equal to one, for a dimensionless calculation
double h_bar=1.0;
double a, beta, alpha; //Constants to be defined later
double a_ho, r; //Diameter of bosons
double m=1, w_ho=1, w_z=1;

//Choose if the trap is speherical (S) or elliptical (E)
char trap='S';
//Choose if we have a non-interaction case (0) or an interaction (1)
int interaction=0;


int main() {
    declare_trap(trap, interaction);
    cout << "Hello World! hei\n";
    double x=distance_particle(3,2,0,9,7,0);
    cout<<x;
    cout << "Hello World! hei\n";
    cout<<beta;

    return 0;
}
