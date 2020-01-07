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
  
  /*main loop*/
  iter =0;
  do{

    /*Time march*/
    for(i=1; i<imax-1; ++i){
    un[i] = u[i] -(c*dt/dx)*(u[i]-u[i-1]);
    }

    /*Update*/
    u = un; // I love vector!

    output(iter, xg, u); //file output
    
    iter +=1;    
    
  }while(iter<itermax);
  
}
