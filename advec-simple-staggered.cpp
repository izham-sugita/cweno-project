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

  float dudxL, dudxR;
  float dudx_ave;
  float ucenter;
  
  do{
   
   for(int i=3; i<imax-3; ++i){
     //source = - ( c*u[i+1] - c*u[i] )/dx ;
     dudxL = 0.5*(c*u[i+1] - c*u[i-1] ) /dx;
     dudxR = 0.5*(c*u[i+2] - c*u[i] ) /dx;
     dudx_ave = 0.5*dudxL + 0.5*dudxR;
     source = -dudx_ave;
     
     uin  = 0.5*(u[i] + u[i+1]) ;
     u1[i] = rkstep(uin, dt, 0.5, source); //update uhalf
   }

   /*Update u[i] from u1[i], which is a uhalf*/
   for(int i=2; i<imax-2; ++i){
     u2[i] = 0.5* (u1[i] + u1[i - 1]); //reconstruct from uhalf
   }
   /*End first RK*/
   
   
   /*Second RK integration loop*/
    for(int i=3; i<imax-3; ++i){
      //      source = - ( c*u1[i+1] - c*u1[i] )/dx ;
     dudxL = 0.5*(c*u2[i+1] - c*u2[i-1] ) /dx;
     dudxR = 0.5*(c*u2[i+2] - c*u2[i] ) /dx;
     dudx_ave = 0.5*dudxL + 0.5*dudxR;
     source = -dudx_ave;
      
     uin  = 0.5*(u[i] + u[i+1]) ;
     u1[i] = rkstep(uin, dt, 1.0, source); //update uhalf

   }
    
   /*Update u[i] from u1[i], which is a uhalf*/
   for(int i=2; i<imax-2; ++i){
     u2[i] = 0.5* (u1[i] + u1[i - 1]); //reconstruct from uhalf
     
   }
   
   swap(u,u2);
     
     if(iter%steps == 0){ output( iter, xg, u );}
     iter +=1;
     
 }while(iter<itermax+1);

 
}
