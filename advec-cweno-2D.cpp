#include "basic.h"

using namespace std; // basic.h provides basic C++ lib. for my project

float cweno4(float stn[5]);
/*sign function*/
int direction(float val);

float rkstep(float u, float dt,
	     float alpha, float source);

void output2D(int steps, int imax, int jmax, 
              vector< vector<float> > x, 
              vector< vector<float> > y, 
              vector< vector<float> > z, 
              vector< vector<float> > fx);


int main()
{

  int imax =201, jmax = 201;
  int i, j;

  float dx = 1.0/(imax-1);
  float dy = 1.0/(jmax-1);
  float dt = 0.0001;
  int itermax = 1000;
  
  cout<<"dx="<<dx<<endl;
  cout<<"Enter dt"<<endl;
  cin>>dt;
  cout<<"Enter maximum iteration"<<endl;
  cin>>itermax;
  
  float xc = 0.5;
  float yc = 0.5;

  float omega = 2.0*pi; // rad/secs
  
  vector< vector<float> > u, v, f, f1, f2, fn;
  vector< vector<float> > xg, yg;
  
  u.resize(imax);
  v.resize(imax);
  f.resize(imax);
  f1.resize(imax);
  f2.resize(imax);
  fn.resize(imax);

  xg.resize(imax);
  yg.resize(imax);

  for(i =0; i<imax; ++i){
    u[i].resize(jmax);
    v[i].resize(jmax);
    f[i].resize(jmax);
    f1[i].resize(jmax);
    f2[i].resize(jmax);
    fn[i].resize(jmax);

    xg[i].resize(jmax);
    yg[i].resize(jmax);
  }

  for(i=0; i<imax; ++i){
    for(j=0; j<jmax; ++j){
      xg[i][j] = i*dx;
      yg[i][j] = j*dy;

      /*Rotational speed*/
      
      u[i][j] = -(yg[i][j] - yc)*omega;
      v[i][j] = (xg[i][j] - xc)*omega;
      
      if( xg[i][j] >=0.2 && xg[i][j] <=0.4){
	if(yg[i][j] >=0.2 && yg[i][j]<=0.4){
	  f[i][j] = 1.0;
	}
      }else{
	f[i][j] = 0.0;
      }
    }
  }

  /*Initiate value*/
  f1 = f;
  f2 = f;
  
  int iter = 0;
  
  float dfdx, dfdy;
  int shift;
  float stn[5];
  float source;
  float fin, ftemp;

  int istart=3;
  int iend = imax-4;
  int jstart=3;
  int jend = jmax-4;
  int steps=10;
  
  float fmhalf, fphalf;

  float fleft, fright, fl, fr;
  float uright, uleft;  
  do{

    for(i=istart; i<iend; ++i){
      for(j=jstart; j<jend; ++j){

	// x-axis for dfdx
      /* for flux at i-(1/2) */
     stn[0] = u[i - 3][j] * f[i - 3][j];
     stn[1] = u[i - 2][j] * f[i - 2][j];
     stn[2] = u[i - 1][j] * f[i - 1][j];
     stn[3] = u[i][j] * f[i][j];
     stn[4] = u[i + 1][j]*f[i + 1][j];
     fleft = cweno4( stn ); //flux

     stn[0] = f[i - 3][j];
     stn[1] = f[i - 2][j];
     stn[2] = f[i - 1][j];
     stn[3] = f[i][j];
     stn[4] = f[i + 1][j];
     uleft = cweno4( stn ); //flux

     stn[0] = u[i + 2][j] * f[i + 2][j];
     stn[1] = u[i + 1][j] * f[i + 1][j];
     stn[2] = u[i][j] * f[i][j];
     stn[3] = u[i - 1][j]*f[i - 1][j];
     stn[4] = u[i - 2][j]*f[i - 2][j];
     fright = cweno4( stn );

     stn[0] = f[i + 2][j];
     stn[1] = f[i + 1][j];
     stn[2] = f[i][j];
     stn[3] = f[i - 1][j];
     stn[4] = f[i - 2][j];
     uright = cweno4( stn );
     
     fl = 0.5*( fleft + fright ) - 0.5* max( abs( u[i][j] ), abs( u[i-1][j] ) ) * ( uright - uleft );

     /* for flux at i+(1/2) */
     stn[0] = u[i - 2][j] * f[i - 2][j];
     stn[1] = u[i - 1][j] * f[i - 1][j];
     stn[2] = u[i][j] *  f[i][j];
     stn[3] = u[i + 1][j] *  f[i + 1][j];
     stn[4] = u[i + 2][j] *  f[i + 2][j];
     fleft = cweno4( stn );

     stn[0] = f[i - 2][j];
     stn[1] = f[i - 1][j];
     stn[2] = f[i][j];
     stn[3] = f[i + 1][j];
     stn[4] = f[i + 2][j];
     uleft = cweno4( stn );

     stn[0] = u[i+3][j] * f[i+3][j];
     stn[1] = u[i+2][j] * f[i+2][j];
     stn[2] = u[i+1][j] * f[i+1][j];
     stn[3] = u[i][j] * f[i][j];
     stn[4] = u[i-1][j] * f[i-1][j];
     fright = cweno4( stn );

     stn[0] = f[i+3][j];
     stn[1] = f[i+2][j];
     stn[2] = f[i+1][j];
     stn[3] = f[i][j];
     stn[4] = f[i-1][j];
     uright = cweno4( stn );

     fr = 0.5*( fleft + fright ) - 0.5* max(abs( u[i][j] ), abs( u[i+1][j] )) *( uright - uleft );

     dfdx = (fr - fl)/dx; 
     /* --------------------------------------------------------------------- */

      // y-axis for dfdy
      /* for flux at j-(1/2) */
     stn[0] = v[i][j - 3] * f[i][j - 3];
     stn[1] = v[i][j - 2] * f[i][j - 2];
     stn[2] = v[i][j - 1] * f[i][j - 1];
     stn[3] = v[i][j] * f[i][j];
     stn[4] = v[i][j+1] * f[i][j + 1];
     fleft = cweno4( stn );

     stn[0] = f[i][j - 3];
     stn[1] = f[i][j - 2];
     stn[2] = f[i][j - 1];
     stn[3] = f[i][j];
     stn[4] = f[i][j + 1];
     uleft = cweno4( stn );

     stn[0] = v[i][j + 2] * f[i][j + 2];
     stn[1] = v[i][j + 1] * f[i][j + 1];
     stn[2] = v[i][j] * f[i][j];
     stn[3] = v[i][j - 1] * f[i][j - 1];
     stn[4] = v[i][j - 2] * f[i][j - 2];
     fright = cweno4( stn );

     stn[0] = f[i][j + 2];
     stn[1] = f[i][j + 1];
     stn[2] = f[i][j];
     stn[3] = f[i][j - 1];
     stn[4] = f[i][j - 2];
     uright = cweno4( stn );
     
     fl = 0.5*( fleft + fright ) - 0.5*max( abs( v[i][j] ), abs( v[i][j-1] ) )*( uright - uleft );

     /* for flux at j+(1/2) */
     stn[0] = v[i][j - 2] * f[i][j - 2];
     stn[1] = v[i][j - 1] * f[i][j - 1];
     stn[2] = v[i][j] * f[i][j];
     stn[3] = v[i][j + 1] * f[i][j + 1];
     stn[4] = v[i][j + 2] * f[i][j + 2];
     fleft = cweno4( stn );

     stn[0] = f[i][j - 2];
     stn[1] = f[i][j - 1];
     stn[2] = f[i][j];
     stn[3] = f[i][j + 1];
     stn[4] = f[i][j + 2];
     uleft = cweno4( stn );

     stn[0] = v[i][j + 3] * f[i][j + 3];
     stn[1] = v[i][j + 2] * f[i][j + 2];
     stn[2] = v[i][j + 1] * f[i][j + 1];
     stn[3] = v[i][j] * f[i][j];
     stn[4] = v[i][j - 1] * f[i][j - 1];
     fright = cweno4( stn );

     stn[0] = f[i][j + 3];
     stn[1] = f[i][j + 2];
     stn[2] = f[i][j + 1];
     stn[3] = f[i][j];
     stn[4] = f[i][j - 1];
     uright = cweno4( stn );
     
     fr = 0.5*( fleft + fright ) - 0.5*max( abs( v[i][j] ), abs( v[i][j+1] ) )*( uright - uleft );

     dfdy = (fr - fl)/dy; 
     /* --------------------------------------------------------------------- */

	/*Time integration*/
	 source = - dfdx -  dfdy ; 
	 fin = f[i][j];
	 ftemp = rkstep(fin, dt, 0.5, source);
	 f1[i][j] = ftemp;
	 	 
	}
      } //end of i-j loop

    /* ---- 2nd-step RK-------- */

      for(i=istart; i<iend; ++i){
      for(j=jstart; j<jend; ++j){

	// x-axis for dfdx
      /* for flux at i-(1/2) */
     stn[0] = u[i - 3][j] * f1[i - 3][j];
     stn[1] = u[i - 2][j] * f1[i - 2][j];
     stn[2] = u[i - 1][j] * f1[i - 1][j];
     stn[3] = u[i][j] * f1[i][j];
     stn[4] = u[i + 1][j]*f1[i + 1][j];
     fleft = cweno4( stn ); //flux

     stn[0] = f1[i - 3][j];
     stn[1] = f1[i - 2][j];
     stn[2] = f1[i - 1][j];
     stn[3] = f1[i][j];
     stn[4] = f1[i + 1][j];
     uleft = cweno4( stn ); //flux

     stn[0] = u[i + 2][j] * f1[i + 2][j];
     stn[1] = u[i + 1][j] * f1[i + 1][j];
     stn[2] = u[i][j] * f1[i][j];
     stn[3] = u[i - 1][j]*f1[i - 1][j];
     stn[4] = u[i - 2][j]*f1[i - 2][j];
     fright = cweno4( stn );

     stn[0] = f1[i + 2][j];
     stn[1] = f1[i + 1][j];
     stn[2] = f1[i][j];
     stn[3] = f1[i - 1][j];
     stn[4] = f1[i - 2][j];
     uright = cweno4( stn );
     
     fl = 0.5*( fleft + fright ) - 0.5* max( abs( u[i][j] ), abs( u[i-1][j] ) ) * ( uright - uleft );

     /* for flux at i+(1/2) */
     stn[0] = u[i - 2][j] * f1[i - 2][j];
     stn[1] = u[i - 1][j] * f1[i - 1][j];
     stn[2] = u[i][j] *  f1[i][j];
     stn[3] = u[i + 1][j] *  f1[i + 1][j];
     stn[4] = u[i + 2][j] *  f1[i + 2][j];
     fleft = cweno4( stn );

     stn[0] = f1[i - 2][j];
     stn[1] = f1[i - 1][j];
     stn[2] = f1[i][j];
     stn[3] = f1[i + 1][j];
     stn[4] = f1[i + 2][j];
     uleft = cweno4( stn );

     stn[0] = u[i+3][j] * f1[i+3][j];
     stn[1] = u[i+2][j] * f1[i+2][j];
     stn[2] = u[i+1][j] * f1[i+1][j];
     stn[3] = u[i][j] * f1[i][j];
     stn[4] = u[i-1][j] * f1[i-1][j];
     fright = cweno4( stn );

     stn[0] = f1[i+3][j];
     stn[1] = f1[i+2][j];
     stn[2] = f1[i+1][j];
     stn[3] = f1[i][j];
     stn[4] = f1[i-1][j];
     uright = cweno4( stn );

     fr = 0.5*( fleft + fright ) - 0.5* max(abs( u[i][j] ), abs( u[i+1][j] )) *( uright - uleft );

     dfdx = (fr - fl)/dx; 
     /* --------------------------------------------------------------------- */

      // y-axis for dfdy
      /* for flux at j-(1/2) */
     stn[0] = v[i][j - 3] * f1[i][j - 3];
     stn[1] = v[i][j - 2] * f1[i][j - 2];
     stn[2] = v[i][j - 1] * f1[i][j - 1];
     stn[3] = v[i][j] * f1[i][j];
     stn[4] = v[i][j+1] * f1[i][j + 1];
     fleft = cweno4( stn );

     stn[0] = f1[i][j - 3];
     stn[1] = f1[i][j - 2];
     stn[2] = f1[i][j - 1];
     stn[3] = f1[i][j];
     stn[4] = f1[i][j + 1];
     uleft = cweno4( stn );

     stn[0] = v[i][j + 2] * f1[i][j + 2];
     stn[1] = v[i][j + 1] * f1[i][j + 1];
     stn[2] = v[i][j] * f1[i][j];
     stn[3] = v[i][j - 1] * f1[i][j - 1];
     stn[4] = v[i][j - 2] * f1[i][j - 2];
     fright = cweno4( stn );

     stn[0] = f1[i][j + 2];
     stn[1] = f1[i][j + 1];
     stn[2] = f1[i][j];
     stn[3] = f1[i][j - 1];
     stn[4] = f1[i][j - 2];
     uright = cweno4( stn );
     
     fl = 0.5*( fleft + fright ) - 0.5*max( abs( v[i][j] ), abs( v[i][j-1] ) )*( uright - uleft );

     /* for flux at j+(1/2) */
     stn[0] = v[i][j - 2] * f1[i][j - 2];
     stn[1] = v[i][j - 1] * f1[i][j - 1];
     stn[2] = v[i][j] * f1[i][j];
     stn[3] = v[i][j + 1] * f1[i][j + 1];
     stn[4] = v[i][j + 2] * f1[i][j + 2];
     fleft = cweno4( stn );

     stn[0] = f1[i][j - 2];
     stn[1] = f1[i][j - 1];
     stn[2] = f1[i][j];
     stn[3] = f1[i][j + 1];
     stn[4] = f1[i][j + 2];
     uleft = cweno4( stn );

     stn[0] = v[i][j + 3] * f1[i][j + 3];
     stn[1] = v[i][j + 2] * f1[i][j + 2];
     stn[2] = v[i][j + 1] * f1[i][j + 1];
     stn[3] = v[i][j] * f1[i][j];
     stn[4] = v[i][j - 1] * f1[i][j - 1];
     fright = cweno4( stn );

     stn[0] = f1[i][j + 3];
     stn[1] = f1[i][j + 2];
     stn[2] = f1[i][j + 1];
     stn[3] = f1[i][j];
     stn[4] = f1[i][j - 1];
     uright = cweno4( stn );
     
     fr = 0.5*( fleft + fright ) - 0.5*max( abs( v[i][j] ), abs( v[i][j+1] ) )*( uright - uleft );

     dfdy = (fr - fl)/dy; 
     /* --------------------------------------------------------------------- */

	/*Time integration*/
	 source = - dfdx -  dfdy ; 
	 fin = f[i][j];
	 ftemp = rkstep(fin, dt, 1.0, source);
	 f2[i][j] = ftemp;
	 	 
	}
      } //end of i-j loop
    
    
    /*Update*/
    swap(f,f2);

    if(iter%steps == 0){
    output2D(iter, imax, jmax, 
              xg, 
              yg, 
              f, 
              f);
    }
    
    iter +=1;
  }while(iter<itermax);
  
    
  
}
