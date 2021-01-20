//Includes necessary header files
#include "utils.h"

//Includes necessary packages
#include <iostream>
#include <math.h>

/*
//#define Globalvariables søk på det
//Constants, most of them is set equal to one, for a dimensionless calculation
double h_bar=1.0;
double a, beta, alpha; //Constants to be defined later
double a_ho; //Diameter of bosons
double m=1, w_ho=1, w_z=1;

//Choose if the trap is speherical (S) or elliptical (E)
char trap='S';
//Choose if we have a non-interaction case (0) or an interaction (1)
int interaction=0;

if (trap=='S'){
  beta=1.0;
}
else if(trap=='E'){
  beta=2.82843; //REMEMBER to change beta value!!
}

if (interaction==0){
  a=0.0
  alpha=0.5/(a_ho*a_ho);
}
else if(interaction==1){
  a=1.0//THIS IS WRONG JUST WROTE A NUMBER
  alpha=; //REMEMBER to change beta value!!
}

//r in V_ext, no idea if its the position or the distance between particles?
//6 or 3 position arguments
double r;
*/


//Function that returns the local energy
double local_energy(int x, int y){
  int z;
  z=x+y;

  return z;
}


double V_external(double x, double y, double z){
  double V_ext;
  if (trap=='S'){
    V_ext=0.5*m*w_ho*w_ho*r*r;
  }
  else if(trap=='E'){
    V_ext=0.5*m*(w_ho*w_ho*(x*x+y*y)+w_z*w_z*z*z);
  }

  return V_ext;
}

double V_internal(double r_distance){
  double V_int;
  if (r_distance>a){
    V_int=0;
  }
  else{
    V_int=std::numeric_limits<double>::infinity();
;
  }

  return V_int;
  }


//Function that calculates the distance between two particles
double distance_particle(double x1,double y1,double z1, double x2, double y2, double z2){
  double distance=sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2));

  return distance;
}

double g_function(double x, double y, double z){
  double g=exp(-alpha*(x*x+y*y+beta*z*z));
  return g;
}

double f_function(double a, double r_dist){
  double f;
  if (r_dist>a){
    f=1-a/abs(r_dist);
  }
  else{
    f=0;
}
  return f;
}

//Declaring some constants
void declare_trap(char ho_trap, int interacting){
  if (ho_trap=='S'){
    beta=1.0;
  }
  else if(ho_trap=='E'){
    beta=2.82843; //REMEMBER to change beta value!!
  }

  if (interacting==0){
    a=0.0;
    alpha=0.5/(a_ho*a_ho);
  }
  else if(interacting==1){
    a=1.0;//THIS IS WRONG JUST WROTE A NUMBER
    alpha=1.0; //REMEMBER to change the alpha value!!
  }


}

//Calculating the derivatives:

void derivative(double x, double y, double z){
  //Forward euler??



}
