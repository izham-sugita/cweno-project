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
        
     fl = 0.5*( c*u[i] + c*u[i-1] ) - 0.5*abs(c)*( u[i] - u[i-1] );
     fr = 0.5*( c*u[i] + c*u[i+1] ) - 0.5*abs(c)*( u[i+1] - u[i] );

     //u1[i] = u[i] - c*(dt/dx)*(up)*( u[i] - u[i-up] );
     
     u1[i] = u[i] - (dt/dx)*( fr - fl );

   }
   swap(u,u1);
     
     if(iter%steps == 0){ output( iter, xg, u );}
     iter +=1;
     
 }while(iter<itermax+1);

 
}
