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
      /*
      u[i][j] = -(yg[i][j] - yc)*omega;
      v[i][j] = (xg[i][j] - xc)*omega;
      */

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
  
  int iter = 0;
  
  float dfdx, dfdy;
  int shift;
  float stn[5];
  float source;
  float fin, ftemp;

  int istart=3;
  int iend = imax-2;
  int jstart=3;
  int jend = jmax-2;
  int steps=10;
  
  float fmhalf, fphalf;
  
  do{

    for(i=istart; i<iend; ++i){
      for(j=jstart; j<jend; ++j){
	//check for u-direction
	if(u[i][j] >= 0.0){
	  stn[0] = f[i-3][j];
	  stn[1] = f[i-2][j];
	  stn[2] = f[i-1][j];
	  stn[3] = f[i][j];
	  stn[4] = f[i+1][j];
	  fmhalf = cweno4(stn);
	  dfdx = 2.0*(f[i][j] - fmhalf)/dx;

	}else{

	  stn[0] = f[i-2][j];
	  stn[1] = f[i-1][j];
	  stn[2] = f[i][j];
	  stn[3] = f[i+1][j];
	  stn[4] = f[i+2][j];
	  fphalf = cweno4(stn);
	  dfdx = 2.0*(fphalf - f[i][j])/dx;

	}

	//check for v-direction
	if(v[i][j] >= 0.0){
	  stn[0] = f[i][j-3];
	  stn[1] = f[i][j-2];
	  stn[2] = f[i][j-1];
	  stn[3] = f[i][j];
	  stn[4] = f[i][j+1];
	  fmhalf = cweno4(stn);
	  dfdy = 2.0*(f[i][j] - fmhalf)/dy;

	}else{

	  stn[0] = f[i][j-2];
	  stn[1] = f[i][j-1];
	  stn[2] = f[i][j];
	  stn[3] = f[i][j+1];
	  stn[4] = f[i][j+2];
	  fphalf = cweno4(stn);
	  dfdy = 2.0*(fphalf - f[i][j])/dy;

	}

	/*Time integration*/
	 source = -u[i][j]*dfdx - v[i][j]*dfdy; 
	 fin = f[i][j];
	 ftemp = rkstep(fin, dt, 1.0, source);
	 f1[i][j] = ftemp;
	 	 
	}
      } //end of i-j loop

    /*Second block, second step RK*/

           //start of i-j loop
    /*
    for(i=istart; i<iend; ++i){
	for(j=jstart; j<jend; ++j){
	
	
	dfdx = 0.5*(udfdx_phalf + udfdx_mhalf);

	
	dfdy = 0.5*(udfdx_phalf + udfdx_mhalf);

	 source = -dfdx - dfdy; //minus is included in reconstruction step
	 fin = f[i][j];
	 f2[i][j] = rkstep(fin, dt,1.0, source);
	 
	}
      } //end of i-j loop
      
    */

    /*Update*/
    swap(f,f1);

    if(iter%steps == 0){
    output2D(iter, imax, jmax, 
              xg, 
              yg, 
              f, 
              f);
    }
    
    iter +=1;
  }while(iter<itermax);
  
  
  /*file output*/
  /*
  ofstream fp;
  fp.open("init-2D.csv", ios::out);
  fp<<"x,y,z,u,v,f\n";
  for(i=0; i<imax; ++i){
    for(j=0; j<jmax; ++j){
      fp<<xg[i][j]<<","
	<<yg[i][j]<<","
	<<0.0<<","
	<<u[i][j]<<","
	<<v[i][j]<<","
	<<f[i][j]<<"\n";
    }
  }
  fp.close();
  */
  
}
