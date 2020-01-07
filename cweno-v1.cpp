#include "basic.h"

using namespace std;


void output(int, vector<float>, vector<float>); //file output

/*RK-stage step: u, dt, alpha, source-rhs*/
float rkstep(float, float,float, float);

void init_calc(int &imax, float &length,
	       float &dx, vector<float> &xg,
	       vector<float> &u);

void param(int &itermax, float &dt)
{
  cout<<"Enter dt\n";
  cin>>dt;
  cout<<"Enter itermax\n";
  cin>>itermax;
}

int main()
{
  int imax;
  float length,dx;
  vector<float> xg, u;
  
  init_calc(imax, length, dx, xg, u);

  cout<<"dx="<<dx<<endl;
  
  int iter,itermax;
  float dt;
  param(itermax,dt);

  vector<float> un,uc;
  un.resize(imax);
  uc.resize(imax);

  un = u;
  uc = u;

  int i;
  float c = 1.0;

  float utemp, uin, source;

  /*Coefficient for 4th CD*/
  float a0 = 1.0/12.0;
  float a1 = 2.0/3.0;
  
  /*main loop*/
  iter =0;
  do{

    /*Constructing value for i+(1/2) */
    //for(i=1; i<imax-1; ++i){
    //uc[i] = 0.5*(u[i+1] + u[i]);
    //}

    /*Constructing value at higher-order*/
    for(i=2; i<imax-2; ++i){
      uc[i] = (-u[i-2] + 7.0*u[i] + 7.0*u[i+1] -u[i+2])/12.0;
    }
    
    /*Time march for constructured value at i+(1/2) */
    //for(i=1; i<imax-1; ++i){
    //uc[i] = uc[i] -(0.5*c*dt/dx)*(u[i+1]-u[i-1]);
    //}

    /*First step for time marching*/
    for(i=2; i<imax-2; ++i){
      uin = uc[i];
      source = -(c/dx)*(a0*u[i-2] - a1*u[i-1] + a1*u[i+1] - a0*u[i+2]);
      utemp = rkstep(uin,dt,1.0,source);
      uc[i] = utemp;
    }

    /*Reconstructing value for i */
    //for(i=1; i<imax-1; ++i){
    //  u[i] = 0.5*(uc[i-1] + uc[i]);
    //}

    /*Reconstruct for higher-order*/
    for(i=2; i<imax-2; ++i){
      u[i] = (-uc[i-2] + 7.0*uc[i] + 7.0*uc[i+1] -uc[i+2])/12.0;
    }
    
    output(iter, xg, u); //file output
    
    iter +=1;    
    
  }while(iter<itermax);
  
}
