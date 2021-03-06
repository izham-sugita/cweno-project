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

  float dudx_i0, dudx_ip1;
  float dudx_iphalf;
  float ucenter;
  
  do{
   
   for(int i=3; i<imax-3; ++i){

     /* for derivative reconstruction at i */
     stn[0] = (c*u[i - 2] - c*u[i-3]) / dx ;
     stn[1] = (c*u[i - 1] - c*u[i - 2]) / dx ;
     stn[2] = (c*u[i] - c*u[i-1]) / dx;
     stn[3] = (c*u[i + 1] - c*u[i]) / dx;
     stn[4] = (c*u[i + 2] - c*u[i + 1]) /dx;
     dudx_i0 = cweno4( stn );

      /* for derivative reconstruction at i+1 */
     stn[0] = (c*u[i - 1] - c*u[i - 2]) / dx ;
     stn[1] = (c*u[i] - c*u[i - 1]) / dx ;
     stn[2] = (c*u[i+1] - c*u[i]) / dx;
     stn[3] = (c*u[i + 2] - c*u[i + 1]) / dx;
     stn[4] = (c*u[i + 3] - c*u[i + 2]) /dx;
     dudx_ip1 = cweno4( stn );

     dudx_iphalf = 0.5*(dudx_i0 + dudx_ip1);

     
     stn[0] = u[i - 2];
     stn[1] = u[i - 1];
     stn[2] = u[i];
     stn[3] = u[i + 1]; 
     stn[4] = u[i + 2];
     ucenter = cweno4(stn);
     
     source = - dudx_iphalf ;
     // uin  = 0.5*(u[i] + u[i+1]) ; //uhalf -> quite diffusive, only 2nd-order approx.
     uin  = ucenter; // uhalf from cweno4 ->less diffusive
     u1[i] = rkstep(uin, dt, 0.25, source); //update uhalf

   }

   /*Update u[i] from u1[i], which is a uhalf*/
   for(int i=3; i<imax-3; ++i){

     // u2[i] = 0.5* (u1[i] + u1[i - 1]); //reconstruct from uhalf

     u2[i+1] = 0.5* (u1[i] + u1[i + 1]); //reconstruct from uhalf

   }
   /*End first RK*/
   
   
   /*Second RK integration loop*/
    for(int i=3; i<imax-3; ++i){

       /* for derivative reconstruction at i */
     stn[0] = (c*u2[i - 2] - c*u2[i-3]) / dx ;
     stn[1] = (c*u2[i - 1] - c*u2[i - 2]) / dx ;
     stn[2] = (c*u2[i] - c*u2[i-1]) / dx;
     stn[3] = (c*u2[i + 1] - c*u2[i]) / dx;
     stn[4] = (c*u2[i + 2] - c*u2[i + 1]) /dx;
     dudx_i0 = cweno4( stn );

      /* for derivative reconstruction at i+1 */
     stn[0] = (c*u2[i - 1] - c*u2[i - 2]) / dx ;
     stn[1] = (c*u2[i] - c*u2[i - 1]) / dx ;
     stn[2] = (c*u2[i+1] - c*u2[i]) / dx;
     stn[3] = (c*u2[i + 2] - c*u2[i + 1]) / dx;
     stn[4] = (c*u2[i + 3] - c*u2[i + 2]) /dx;
     dudx_ip1 = cweno4( stn );

     dudx_iphalf = 0.5*(dudx_i0 + dudx_ip1);     

     stn[0] = u[i - 2];
     stn[1] = u[i - 1];
     stn[2] = u[i];
     stn[3] = u[i + 1]; 
     stn[4] = u[i + 2];
     ucenter = cweno4(stn);
     
     source = - dudx_iphalf ;
     uin  = ucenter ; // uhalf from cweno
     u1[i] = rkstep(uin, dt, 0.33333333, source); //update uhalf

   }
    
   /*Update u[i] from u1[i], which is a uhalf*/
   for(int i=3; i<imax-3; ++i){
     //     u2[i] = 0.5* (u1[i] + u1[i - 1]); //reconstruct from uhalf
     u2[i+1] = 0.5* (u1[i] + u1[i + 1]); //reconstruct from uhalf
   }

   
   /*Third RK integration loop*/
    for(int i=3; i<imax-3; ++i){

       /* for derivative reconstruction at i */
     stn[0] = (c*u2[i - 2] - c*u2[i-3]) / dx ;
     stn[1] = (c*u2[i - 1] - c*u2[i - 2]) / dx ;
     stn[2] = (c*u2[i] - c*u2[i-1]) / dx;
     stn[3] = (c*u2[i + 1] - c*u2[i]) / dx;
     stn[4] = (c*u2[i + 2] - c*u2[i + 1]) /dx;
     dudx_i0 = cweno4( stn );

      /* for derivative reconstruction at i+1 */
     stn[0] = (c*u2[i - 1] - c*u2[i - 2]) / dx ;
     stn[1] = (c*u2[i] - c*u2[i - 1]) / dx ;
     stn[2] = (c*u2[i+1] - c*u2[i]) / dx;
     stn[3] = (c*u2[i + 2] - c*u2[i + 1]) / dx;
     stn[4] = (c*u2[i + 3] - c*u2[i + 2]) /dx;
     dudx_ip1 = cweno4( stn );

     dudx_iphalf = 0.5*(dudx_i0 + dudx_ip1);     

     stn[0] = u[i - 2];
     stn[1] = u[i - 1];
     stn[2] = u[i];
     stn[3] = u[i + 1]; 
     stn[4] = u[i + 2];
     ucenter = cweno4(stn);
     
     source = - dudx_iphalf ;
     uin  = ucenter ; // uhalf from cweno
     u1[i] = rkstep(uin, dt, 0.5, source); //update uhalf

   }
    
   /*Update u[i] from u1[i], which is a uhalf*/
   for(int i=3; i<imax-3; ++i){
     //     u2[i] = 0.5* (u1[i] + u1[i - 1]); //reconstruct from uhalf
     u2[i+1] = 0.5* (u1[i] + u1[i + 1]); //reconstruct from uhalf
   }

   
    /*4th RK integration loop*/
    for(int i=3; i<imax-3; ++i){

       /* for derivative reconstruction at i */
     stn[0] = (c*u2[i - 2] - c*u2[i-3]) / dx ;
     stn[1] = (c*u2[i - 1] - c*u2[i - 2]) / dx ;
     stn[2] = (c*u2[i] - c*u2[i-1]) / dx;
     stn[3] = (c*u2[i + 1] - c*u2[i]) / dx;
     stn[4] = (c*u2[i + 2] - c*u2[i + 1]) /dx;
     dudx_i0 = cweno4( stn );

      /* for derivative reconstruction at i+1 */
     stn[0] = (c*u2[i - 1] - c*u2[i - 2]) / dx ;
     stn[1] = (c*u2[i] - c*u2[i - 1]) / dx ;
     stn[2] = (c*u2[i+1] - c*u2[i]) / dx;
     stn[3] = (c*u2[i + 2] - c*u2[i + 1]) / dx;
     stn[4] = (c*u2[i + 3] - c*u2[i + 2]) /dx;
     dudx_ip1 = cweno4( stn );

     dudx_iphalf = 0.5*(dudx_i0 + dudx_ip1);     

     stn[0] = u[i - 2];
     stn[1] = u[i - 1];
     stn[2] = u[i];
     stn[3] = u[i + 1]; 
     stn[4] = u[i + 2];
     ucenter = cweno4(stn);
     
     source = - dudx_iphalf ;
     uin  = ucenter ; // uhalf from cweno
     u1[i] = rkstep(uin, dt, 1.0, source); //update uhalf

   }
    
   /*Update u[i] from u1[i], which is a uhalf*/
   for(int i=3; i<imax-3; ++i){
     //     u2[i] = 0.5* (u1[i] + u1[i - 1]); //reconstruct from uhalf
     u2[i+1] = 0.5* (u1[i] + u1[i + 1]); //reconstruct from uhalf
   }

   
   swap(u,u2);
     
     if(iter%steps == 0){ output( iter, xg, u );}
     iter +=1;
     
 }while(iter<itermax+1);

 
}
