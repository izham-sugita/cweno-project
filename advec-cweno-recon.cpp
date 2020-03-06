#include<iostream>
#include<algorithm> //for swap()
#include<fstream> //for file I/O
#include<cmath> // math function
#include<vector>

using namespace std;

void output(int steps, vector<float> x, vector<float> fx);

float rkstep(float u, float dt,
	     float alpha, float source);

float cweno4(float stn[5]);

int direction(float val);

int main()
{
  vector<float> uex, u, uinit;
  int imax;
  float dx, dt;

  /*Input*/
  cout<<"Enter imax "<<endl;
  cin>>imax;

  dx = 10.0/(imax-1);
  cout<<"Enter dt (dx= "<<dx<<" ) "<<endl;;
  cin>>dt;

  int itermax;
  cout<<"Enter itermax "<<endl;
  cin>>itermax;
  
  /*Initial condition*/
  uex.resize(imax);
  u.resize(imax);
  uinit.resize(imax);

  vector<float> xg;
  xg.resize(imax);

  for(int i=0; i<imax; ++i){
    xg[i] = i*dx;
  }

  for(int i=0; i<imax; ++i){
    u[i] = 0.0;
    uex[i] = 0.0;
    if( xg[i]>4.0 && xg[i] < 6.0){
      u[i] = 1.0;
    }
  }  
  
  uinit = u;
  /*Creating analytical solution*/
  int iter=0;
  do{

    for(int i=2; i<imax-1; ++i){
      uex[i] = u[i-1];  
    }
    //update
    swap(u,uex);
    iter +=1;

  }while(iter<itermax+1);

  /*starting the WENO scheme*/
  cout<<"Solution marching time for analytical solution"<<endl;
  cout<<itermax*dx<<endl;
  cout<<"Enter new itermax to match"<<endl;
  cin>>itermax;
  
  
  vector<float> u1,u2;
  u1.resize(imax);
  u2.resize(imax);
  

  u1 = uinit;
  u2 = uinit;
  u = uinit;

    float stn[5];
   float fl, fr;
   float fleft, fright;
   float uleft, uright;
   float dfdx;
  
  float source;
  float utemp;
  float uin;
  float alpha;
  int steps = 1;
  iter =0;

  int up;
  float c;
  cout<<"Enter c (+1, or -1)"<<endl;
  cin>>c;

  up = direction(c);  

  do{
   
   for(int i=3; i<imax-3; ++i){
   uin = u[i];
   alpha=0.25f;

   /*------------------------------------------------------*/      
   stn[0] = (u[i-2*up] - u[i-3*up]) / dx ;
   stn[1 ] = up*( u[i-up] - u[i-2*up]) / dx ;
   stn[2 ] = up*( u[i] - u[i-up] ) / dx ; // reconstruction at [i]
   stn[3 ] = up*( u[i+up] - u[i] ) / dx ;
   stn[4 ] = up*( u[i+2*up] - u[i+1*up] ) /dx;
    dfdx = cweno4(stn);
    
   source = -c*dfdx;
   utemp = rkstep(uin,dt,alpha,source);
   u1[i] = utemp;
   } //first step

   for(int i=3; i<imax-3; ++i){
   uin = u[i];
   alpha=0.333333333f;

   /*
   stn[0] = (u1[i-2] - u1[i-3])/dx ;
   stn[1] = (u1[i-1] - u1[i-2])/dx ;
   stn[2] = (u1[i] - u1[i-1])/dx ; // reconstruction at [i]
   stn[3] = (u1[i+1] - u1[i])/dx ;
   stn[4] = (u1[i+2] - u1[i+1])/dx;
    dfdx = cweno4(stn);
    */

   stn[0 ]  =  up*(u1[i-2*up] - u1[i-3*up]) / dx ;
   stn[1 ]      =  up*( u1[i-up] - u1[i-2*up]) / dx ;
   stn[2 ]           =  up*( u1[i] - u1[i-up] ) / dx ; // reconstruction at [i]
   stn[3 ]    =  up*( u1[i+up] - u1[i] ) / dx ;
   stn[4 ] = up*( u1[i+2*up] - u1[i+1*up] ) /dx;
    dfdx = cweno4(stn);
    
   source = -c*dfdx;
   utemp = rkstep(uin,dt,alpha,source);
   u2[i] = utemp;
   } //second step

   for(int i=3; i<imax-3; ++i){
   uin = u[i];
   alpha=0.5f;

   /*
   stn[0] = (u2[i-2] - u2[i-3])/dx ;
   stn[1] = (u2[i-1] - u2[i-2])/dx ;
   stn[2] = (u2[i] - u2[i-1])/dx ; // reconstruction at [i]
   stn[3] = (u2[i+1] - u2[i])/dx ;
   stn[4] = (u2[i+2] - u2[i+1])/dx;
    dfdx = cweno4(stn);
   */

   stn[0 ]  =  up*(u2[i-2*up] - u2[i-3*up]) / dx ;
   stn[1 ]      =  up*( u2[i-up] - u2[i-2*up]) / dx ;
   stn[2 ]           =  up*( u2[i] - u2[i-up] ) / dx ; // reconstruction at [i]
   stn[3 ]    =  up*( u2[i+up] - u2[i] ) / dx ;
   stn[4 ] = up*( u2[i+2*up] - u2[i+1*up] ) /dx;
   dfdx = cweno4(stn);
	  
   source = -c*dfdx; 
   utemp = rkstep(uin,dt,alpha,source);
   u1[i] = utemp;
   } //third step

   for(int i=3; i<imax-3; ++i){
   uin = u[i];
   alpha=1.0f;

   /*
   stn[0] = (u1[i-2] - u1[i-3])/dx ;
   stn[1] = (u1[i-1] - u1[i-2])/dx ;
   stn[2] = (u1[i] - u1[i-1])/dx ; // reconstruction at [i]
   stn[3] = (u1[i+1] - u1[i])/dx ;
   stn[4] = (u1[i+2] - u1[i+1])/dx;
    dfdx = cweno4(stn);
   */

   stn[0 ]  =  up*(u1[i-2*up] - u1[i-3*up]) / dx ;
   stn[1 ]      =  up*( u1[i-up] - u1[i-2*up]) / dx ;
   stn[2 ]           =  up*( u1[i] - u1[i-up] ) / dx ; // reconstruction at [i]
   stn[3 ]    =  up*( u1[i+up] - u1[i] ) / dx ;
   stn[4 ] = up*( u1[i+2*up] - u1[i+1*up] ) /dx;
   dfdx = cweno4(stn);
         
   source = -c*dfdx;
   utemp = rkstep(uin,dt,alpha,source);
   u2[i] = utemp;
   } //forth step

   swap(u,u2);

   if(iter%steps == 0){ output( iter, xg, u );}
   
     iter +=1;
     
 }while(iter<itermax+1);

 
}
