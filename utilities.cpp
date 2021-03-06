/*
This file contains utility functions such as file output,
simple functions etc.
The content will be compiled first and linked later.
 */
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<iomanip>

using namespace std;


/*File output*/
void output(int steps, vector<float> x, vector<float> fx)
{
  ofstream fp;
  stringstream buf;
  string filenumber;

  buf<<setfill('0');
  filenumber = to_string(steps);
  buf<<setw(5)<<filenumber;

  string filename="./data/f"+buf.str()+".csv";
  fp.open(filename, ios::out);
  fp<<"x, fx\n";
  for(int i=0; i<x.size(); ++i){
    fp<<x[i]<<", "
      <<fx[i]<<"\n";
  }
  fp.close();
  buf.str(string()); //clear buffer
}

/*sign function*/
template<typename T>int sign(T val){
  return (T(0) < val) - (val < T(0));
}

/*float sign*/
float sign_f(float val){
 return ( 0.0 < val) - (val < 0.0);
}

/*Runge-Kutta step function*/
float rkstep(float u, float dt,
	     float alpha, float source)
{
  u = u + alpha*dt*source;
  return u;
}

//WENO5 reconstruction procedure.
float weno_recon(float stn[5]) //calculate the right-hand side
{
  float fr;

  float gamma1, gamma2, gamma3;
  gamma1 = 0.1f;
  gamma2 = 0.6f;
  gamma3 = 0.3f;
    
  float p1, p2, p3;

  p1 = (1.0/3.0)*stn[0]+(-7.0/6.0)*stn[1]+(11.0/6.0)*stn[2];
  p2 = (-1.0/6.0)*stn[1]+(5.0/6.0)*stn[2]+(1.0/3.0)*stn[3];
  p3 = (1.0/3.0)*stn[2]+(5.0/6.0)*stn[3]+(-1.0/6.0)*stn[4];

  /*Calculating smoothness indicator */
  float b1 = (13.0/12.0)*(stn[0]-2.0*stn[1]+stn[2])*(stn[0]-2.0*stn[1]+stn[2])
    +0.25*(3.0*stn[0]-4.0*stn[1]+stn[2])*(3.0*stn[0]-4.0*stn[1]+stn[2]);
  
  float b2 = (13.0/12.0)*(stn[1]-2.0*stn[2]+stn[3])*(stn[1]-2.0*stn[2]+stn[3])
    +0.25*(3.0*stn[1]-4.0*stn[2]+stn[3])*(3.0*stn[1]-4.0*stn[2]+stn[3]);
  
  float b3 = (13.0/12.0)*(stn[2]-2.0*stn[3]+stn[4])*(stn[2]-2.0*stn[3]+stn[4])
    +0.25*(3.0*stn[2]-4.0*stn[3]+stn[4])*(3.0*stn[2]-4.0*stn[3]+stn[4]);
  
    
 float wtilde1 = gamma1/( (1.0e-6 + b1)*(1.0e-6 + b1) );
 float wtilde2 = gamma2/( (1.0e-6 + b2)*(1.0e-6 + b2) );
 float wtilde3 = gamma3/( (1.0e-6 + b3)*(1.0e-6 + b3) );
  

 float w1 = wtilde1/(wtilde1+wtilde2+wtilde3);
 float w2 = wtilde2/(wtilde1+wtilde2+wtilde3);
 float w3 = wtilde3/(wtilde1+wtilde2+wtilde3);
    
  fr = w1*p1 + w2*p2 + w3*p3;
  
  return fr;
}

/*Initiate step function*/
void init_calc(int &imax, float &length,
	       float &dx, vector<float> &xg,
	       vector<float> &u)
{

  int i;
  
  cout<<"Enter imax\n";
  cin>>imax;
  cout<<"Enter length\n";
  cin>>length;

  dx = length/(imax-1);

  xg.resize(imax);
  u.resize(imax);
  
  for(i=0; i<imax; ++i){
    xg[i] = i*dx;
    u[i] = 0.0;
    if(xg[i] >=0.4*length && xg[i]<=0.6*length){
      u[i] = 1.0;
    }
  }
  
};
