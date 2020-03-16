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
  
  float xc0 = 0.5;
  float yc0 = 0.5;

  float omega = 2.0*pi; // rad/secs
  
  vector< vector<float> > u, v, f, f1, f2, fn;
  vector< vector<float> > xg, yg;

  vector< vector<float> > fc;
  vector< vector<float> > xc, yc;
  
  u.resize(imax);
  v.resize(imax);
  f.resize(imax);
  f1.resize(imax);
  f2.resize(imax);
  fn.resize(imax);
  
  fc.resize(imax);
  xc.resize(imax);
  yc.resize(imax);
  
  xg.resize(imax);
  yg.resize(imax);

  for(i =0; i<imax; ++i){
    u[i].resize(jmax);
    v[i].resize(jmax);
    f[i].resize(jmax);
    f1[i].resize(jmax);
    f2[i].resize(jmax);
    fn[i].resize(jmax);

    fc[i].resize(jmax);
    xc[i].resize(jmax);
    yc[i].resize(jmax);

    xg[i].resize(jmax);
    yg[i].resize(jmax);
  }

  for(i=0; i<imax; ++i){
    for(j=0; j<jmax; ++j){
      xg[i][j] = i*dx;
      yg[i][j] = j*dy;

      /*Rotational speed*/
      //u[i][j] = -(yg[i][j] - yc0)*omega;
      //v[i][j] = (xg[i][j] - xc0)*omega;

      u[i][j] = 1.0;
      v[i][j] = 1.0;
      
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
  fc = f;

  xc = xg;
  yc = yg;
  
  int iter = 0;
  
  
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

  vector< vector<float> > dfdx, dfdy;
  dfdx.resize(imax);
  dfdy.resize(imax);
  for(i=0; i<imax; ++i){
    dfdx[i].resize(jmax);
    dfdy[i].resize(jmax);
  }
  dfdx = f;
  dfdy = f;

  //calculating position of staggered grid
  for(i=istart; i<iend; ++i){
    for(j=jstart; j<jend; ++j){
      xc[i][j] = xc[i][j] + 0.5*dx;
      yc[i][j] = yc[i][j] + 0.5*dy;
    }
  }

  //time loop will start here 
  do{
  
  //Interpolating direction-by-direction i+1/2, j+1/2
    for(i=istart; i<iend; ++i){
      for(j=jstart; j<jend; ++j){
	//interpolation for i+1/2, j
	stn[0] = f[i-2][j];
	stn[1] = f[i-1][j];
	stn[2] = f[i][j];
	stn[3] = f[i+1][j];
	stn[4] = f[i+2][j];
	f1[i][j] = cweno4(stn);
	}
      } 
    for(i=istart; i<iend; ++i){
      for(j=jstart; j<jend; ++j){
	//interpolation for i+1/2, j+1/2
	stn[0] = f1[i-2][j];
	stn[1] = f1[i-1][j];
	stn[2] = f1[i][j];
	stn[3] = f1[i+1][j];
	stn[4] = f1[i+2][j];
	fc[i][j] = cweno4(stn);
      }
    } 
    //end interpolation

    //reconstruction of derivatives d(uf)dx at (i,j)
       for(i=istart; i<iend; ++i){
      for(j=jstart; j<jend; ++j){
	stn[0] = ( u[i-2][j]*f[i-2][j]  -  u[i-3][j]*f[i-3][j] ) / dx ;
	stn[1] = ( u[i-1][j]*f[i-1][j] -   u[i-2][j]*f[i-2][j] ) / dx;
	stn[2] = ( u[i][j]*f[i][j]   -  u[i-1][j]*f[i-1][j] ) / dx ;
	stn[3] = ( u[i+1][j]*f[i+1][j] - u[i][j]*f[i][j] ) / dx;
	stn[4] = ( u[i+2][j]*f[i+2][j] - u[i+1][j]*f[i+1][j] ) / dx;
	dfdx[i][j] = cweno4(stn);
	}
      }
       //reconstruction of derivatives d(vf)dy at (i,j)
       for(i=istart; i<iend; ++i){
      for(j=jstart; j<jend; ++j){
	stn[0] = ( v[i][j-2]*f[i][j-2] - v[i][j-3]*f[i][j-3] ) / dy ;
	stn[1] = ( v[i][j-1]*f[i][j-1] -  v[i][j-2]*f[i][j-2] ) / dy;
	stn[2] = ( v[i][j]*f[i][j]    -  v[i][j-1]*f[i][j-1] ) / dy ;
	stn[3] = ( v[i][j+1]*f[i][j+1] - v[i][j]*f[i][j] ) / dy;
	stn[4] = ( v[i][j+2]*f[i][j+2] - v[i][j+1]*f[i][j+1] ) / dy;
	dfdy[i][j] = cweno4(stn);
	}
      } 
       //advection step
       for(i = istart; i < iend; ++i){
	 for(j = jstart; j< jend; ++j){
	   source = -0.5*(dfdx[i][j]+dfdx[i+1][j]) -0.5*(dfdy[i][j]+dfdy[i][j+1]);
	   fin = fc[i][j];
	   f2[i][j] = rkstep(fin, dt,
	     0.5, source);
	 }
       }

       //reconstruction
       
       for(i = istart; i<iend; ++i){
	 for(j= jstart; j<jend; ++j){
	   f1[i][j] = 0.25*( f2[i][j] + f2[i-1][j] + f2[i-1][j-1] + f2[i][j-1]  ); 
	 }
       }
       
       
     //Second RK-step
     //reconstructing dfdx, dfdy at (i,j) //reconstruction of derivatives dfdx at (i,j)
       for(i=istart; i<iend; ++i){
      for(j=jstart; j<jend; ++j){
	stn[0] = ( u[i-2][j]*f1[i-2][j]  -  u[i-3][j]*f1[i-3][j] ) / dx ;
	stn[1] = ( u[i-1][j]*f1[i-1][j] -   u[i-2][j]*f1[i-2][j] ) / dx;
	stn[2] = ( u[i][j]*f1[i][j]   -  u[i-1][j]*f1[i-1][j] ) / dx ;
	stn[3] = ( u[i+1][j]*f1[i+1][j] - u[i][j]*f1[i][j] ) / dx;
	stn[4] = ( u[i+2][j]*f1[i+2][j] - u[i+1][j]*f1[i+1][j] ) / dx;
	dfdx[i][j] = cweno4(stn);
	}
      }
       //reconstruction of derivatives dfdy at (i,j)
       for(i=istart; i<iend; ++i){
      for(j=jstart; j<jend; ++j){
	stn[0] = ( v[i][j-2]*f1[i][j-2] - v[i][j-3]*f1[i][j-3] ) / dy ;
	stn[1] = ( v[i][j-1]*f1[i][j-1] -  v[i][j-2]*f1[i][j-2] ) / dy;
	stn[2] = ( v[i][j]*f1[i][j]    -  v[i][j-1]*f1[i][j-1] ) / dy ;
	stn[3] = ( v[i][j+1]*f1[i][j+1] - v[i][j]*f1[i][j] ) / dy;
	stn[4] = ( v[i][j+2]*f1[i][j+2] - v[i][j+1]*f1[i][j+1] ) / dy;
	dfdy[i][j] = cweno4(stn);
	}
      }
        //advection step
       for(i = istart; i < iend; ++i){
	 for(j = jstart; j< jend; ++j){
	   source = -0.5*(dfdx[i][j]+dfdx[i+1][j]) -0.5*(dfdy[i][j]+dfdy[i][j+1]);
	   fin = fc[i][j];
	   f2[i][j] = rkstep(fin, dt, 1.0, source);
	   	 }
       }
          //reconstruction
       
       for(i = istart; i<iend; ++i){
	 for(j= jstart; j<jend; ++j){
	   f1[i][j] = 0.25*( f2[i][j] + f2[i-1][j] + f2[i-1][j-1] + f2[i][j-1]  ); 
	 }
       }
       
       
     // end of second RK step
    // ----------------------------------------------------------------- //    

     swap(f, f1);
     iter = iter+1;

     if(iter%steps == 0) output2D(iter, imax, jmax, xg, yg, f, f);
     
}while(iter < itermax);

     vector< vector<float> > zg;
     zg.resize(imax);
     for(i=0; i<imax; ++i){
       zg[i].resize(jmax);
     }
     for(i=0; i<imax; ++i){
       for(j=0; j<jmax; ++j){
	 zg[i][j] = 0.0;
       }
     }
     
     //output2D(1, imax, jmax, xc, yc, fc, fc);
     //output2D(2, imax, jmax, xg, yg, f, f);
     //output2D(3, imax, jmax, xg, yg, zg, dfdx);

    


  
    
  
}
