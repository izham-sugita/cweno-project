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

float cweno4(float stencil[5])
{
  float um2, um1, u0, up1, up2;
  float eps = 1.0e-6;
  float c0, c1, c2;
  float b0, b1, b2;
  float a0, a1, a2;
  float w0, w1, w2;
  float uhalf;

  um2 = stencil[0];
  um1 = stencil[1];
  u0 = stencil[2];
  up1 = stencil[3];
  up2 = stencil[4];

  eps = 1.0e-6;
  c0 = 1.0/6.0;
  c1 = 2.0/3.0;
  c2 = 1.0/6.0;

  b0 = 0.25*(um2-4.0*um1+3.0*u0)*(um2-4.0*um1+3.0*u0)
    + (13.0/12.0)*(um2-2.0*um1+u0)*(um2-2.0*um1+u0);
  
  b1 = 0.25*(up1-um1)*(up1-um1)
    + (13.0/12.0)*(um1-2.0*u0+up1)*(um1-2.0*u0+up1);
  
  b2 = 0.25*(3.0*u0-4.0*up1+up2)*(3.0*u0-4.0*up1+up2)
    + (13.0/12.0)*(u0-2.0*up1+up2)*(u0-2.0*up1+up2);

  a0 = c0/( ( eps + b0)*( eps + b0)  );
  a1 = c1/( ( eps + b1)*( eps + b1)  );
  a2 = c2/( ( eps + b2)*( eps + b2)  );

  w0 = a0/( a0+a1+a2  );
  w1 = a1/( a0+a1+a2  );
  w2 = a2/( a0+a1+a2  );

  uhalf = 0.5*( -w0*um1 + (3.0*w0+w1)*u0 +(3.0*w2+w1)*up1 -w2*up2  );
    
  return uhalf;
  
}

float cweno4_v2(float stencil[5], float weight[3])
{
  float um2, um1, u0, up1, up2;
  float eps = 1.0e-6;
  float c0, c1, c2;
  float b0, b1, b2;
  float a0, a1, a2;
  float w0, w1, w2;
  float uhalf;

  um2 = stencil[0];
  um1 = stencil[1];
  u0 = stencil[2];
  up1 = stencil[3];
  up2 = stencil[4];

  eps = 1.0e-6;
  c0 = 1.0/6.0;
  c1 = 2.0/3.0;
  c2 = 1.0/6.0;

  b0 = 0.25*(um2-4.0*um1+3.0*u0)*(um2-4.0*um1+3.0*u0)
    + (13.0/12.0)*(um2-2.0*um1+u0)*(um2-2.0*um1+u0);
  
  b1 = 0.25*(up1-um1)*(up1-um1)
    + (13.0/12.0)*(um1-2.0*u0+up1)*(um1-2.0*u0+up1);
  
  b2 = 0.25*(3.0*u0-4.0*up1+up2)*(3.0*u0-4.0*up1+up2)
    + (13.0/12.0)*(u0-2.0*up1+up2)*(u0-2.0*up1+up2);

  a0 = c0/( ( eps + b0)*( eps + b0)  );
  a1 = c1/( ( eps + b1)*( eps + b1)  );
  a2 = c2/( ( eps + b2)*( eps + b2)  );

  w0 = a0/( a0+a1+a2  );
  w1 = a1/( a0+a1+a2  );
  w2 = a2/( a0+a1+a2  );

  weight[0] = w0;
  weight[1] = w1;
  weight[2] = w2;

  uhalf = 0.5*( -w0*um1 + (3.0*w0+w1)*u0 +(3.0*w2+w1)*up1 -w2*up2  );
    
  return uhalf;
  
}

float cweno4_v3(float stencil[5], float weight[3] )
{
  float um2, um1, u0, up1, up2;
  float eps = 1.0e-6;
  float c0, c1, c2;
  float b0, b1, b2;
  float a0, a1, a2;
  float w0, w1, w2;
  float uhalf;

  um2 = stencil[0];
  um1 = stencil[1];
  u0 = stencil[2];
  up1 = stencil[3];
  up2 = stencil[4];

  w0 = weight[0];
  w1 = weight[1];
  w2 = weight[2];
  
  uhalf = 0.5*( -w0*um1 + (3.0*w0+w1)*u0 +(3.0*w2+w1)*up1 -w2*up2  );
    
  return uhalf;
  
}


float sign_f(float);

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

  vector<float> u1,u2, u3;
  u1.resize(imax);
  u2.resize(imax);
  u3.resize(imax);

  vector<float> w[3];
  w[0].resize(imax);
  w[1].resize(imax);
  w[2].resize(imax);

  float weight[3];
  
  un = u;
  uc = u;
  u1 = u;
  u2 = u;
  u3 = u;
 
 
  int i;
  float c = 1.0;

  float alpha = 1.0;
  float utemp, uin, source;

  float uphalf, umhalf;
  float stn[5];
  float u_1;
  float dudx;

  float vec;
  float pi = 3.142;
  int up;
  
  /*main loop*/
  iter =0;
  do{

    vec = 0.5*pi*c*iter*dt;
    vec = sin(vec);

    for(i=3; i<imax-4; ++i){
     stn[0] = u[i-3];
	 stn[1] = u[i-2];
	 stn[2] = u[i-1] ;
	 stn[3] = u[i] ;
	 stn[4] = u[i+1] ;
	 umhalf = cweno4(stn);

 stn[0] = u[i-2];
	 stn[1] = u[i-1];
	 stn[2] = u[i] ;
	 stn[3] = u[i+1] ;
	 stn[4] = u[i+2] ;

	 uphalf = cweno4_v2(stn, &weight[0]);

	 w[0][i] = weight[0];
	 w[1][i] = weight[1];
	 w[2][i] = weight[2];
	 
	 uc[i] = 0.5*(umhalf + uphalf);
	 
    }
    
    /*mid-point value interpolation*/
    for(i =3; i<imax-4; ++i){
            
	 stn[0] = (u[i-2]-u[i-3])/dx;
      stn[1] = (u[i-1]-u[i-2])/dx;
      stn[2] = (u[i]-u[i-1])/dx ;
      stn[3] = (u[i+1]-u[i])/dx ;
      stn[4] = (u[i+2]-u[i+1])/dx ;

       weight[0] = w[0][i];
       weight[1] = w[1][i];
       weight[2] = w[2][i];

       dudx = cweno4_v3(stn, weight);

	 source = -sign_f(vec)*c*dudx;

	 //getting the average right is important!
	 utemp = uc[i];
	 
      u1[i] = rkstep(utemp, dt,0.33333333, source);
      
    }

       /*mid-point value interpolation*/
    for(i =3; i<imax-4; ++i){

      stn[0] = (u1[i-2]-u1[i-3])/dx;
      stn[1] = (u1[i-1]-u1[i-2])/dx;
      stn[2] = (u1[i]-u1[i-1])/dx ;
      stn[3] = (u1[i+1]-u1[i])/dx ;
      stn[4] = (u1[i+2]-u1[i+1])/dx ;

         weight[0] = w[0][i];
	 weight[1] = w[1][i];
	 weight[2] = w[2][i];
	 dudx = cweno4_v3(stn, weight);
	 source = -sign_f(vec)*c*dudx;
	 utemp = uc[i];
	 
      u2[i] = rkstep(utemp, dt,0.5, source);
      
    }

     for(i =3; i<imax-4; ++i){

       stn[0] = (u2[i-2]-u2[i-3])/dx;
       stn[1] = (u2[i-1]-u2[i-2])/dx;
       stn[2] = (u2[i]-u2[i-1])/dx ;
       stn[3] = (u2[i+1]-u2[i])/dx ;
       stn[4] = (u2[i+2]-u2[i+1])/dx ;
         weight[0] = w[0][i];
	 weight[1] = w[1][i];
	 weight[2] = w[2][i];
	 dudx = cweno4_v3(stn, weight);
	 source = -sign_f(vec)*c*dudx;
	 utemp = uc[i];
	 
      u3[i] = rkstep(utemp, dt,1.0, source);
      
    }

    swap(u,u3);
     
    output(iter, xg, u); //file output
    
    iter +=1;    
    
  }while(iter<itermax);
  
  
}
