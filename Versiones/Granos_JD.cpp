//Metodos de Integracion 
#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
using namespace std;

const double g=9.8;
const double K=1e4, Gamma=50, Kcundall=10, MU=0.4; 
const double Lx=100,Ly=100;
const int Nx=1,Ny=1,N=Nx*Ny;

const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;

class Cuerpo;
class Colisionador;

//----------------------Clase Cuerpo---------------
class Cuerpo{
private:
  vector3D r,V,F,omega,tau; double m,R,theta,I;
public:
  void Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,
	      double theta0,double omega0,
	      double m0,double R0);
  void BorreFuerzayTorque(void);
  void AgregueFuerza(vector3D F0);
  void AgregueTorque(vector3D tau0);
  void Mueva_r(double dt,double Constante);
  void Mueva_V(double dt,double Constante);
  void Dibujese(void);
  double Getx(void){return r.x();}; //Función Inline
  double Gety(void){return r.y();}; //Función Inline
  double GetV(void){return norma(V);}; //Función Inline
  double GetVx(void){return V.x();}; //Función Inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,
	      double theta0,double omega0,
	      double m0,double R0){
  r.cargue(x0,y0,z0); V.cargue(Vx0,Vy0,Vz0); omega.cargue(0,0,omega0);
  theta=theta0;  m=m0; R=R0; I=2.0/5*m*R*R;
}
void Cuerpo::BorreFuerzayTorque(void){
  F.cargue(0,0,0);   tau.cargue(0,0,0); 
}
void Cuerpo::AgregueFuerza(vector3D F0){
  F+=F0;
}
void Cuerpo::AgregueTorque(vector3D tau0){
  tau+=tau0;
}
void Cuerpo::Mueva_r(double dt,double Constante){
  r+=V*(Constante*dt); theta+=omega.z()*Constante*dt;
}
void Cuerpo::Mueva_V(double dt,double Constante){
  V+=F*(Constante*dt/m); omega+=tau*(Constante*dt/I);
}
void Cuerpo::Dibujese(void){
  cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
      <<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t"; 
}
//----------------------Clase Colisionador---------------
class Colisionador{
private:
  vector3D ele[N+4][N+4]; bool EstabaEnContacto[N+4][N+4];
public:
  void Inicie(void);
  void CalculeTodasLasFuerzas(Cuerpo* Grano,double dt);
  void CalculeLaFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2,
			    vector3D & ele, bool & EstabaEnContacto,double dt);
};
void Colisionador::Inicie(void){
  int i,j;
  for(i=0;i<N;i++)
    for(j=i+1;j<N+4;j++){
      ele[i][j].cargue(0,0,0); EstabaEnContacto[i][j]=false;
    }
}
void Colisionador::CalculeTodasLasFuerzas(Cuerpo* Grano,double dt){
  int i,j;
  vector3D g_vector; g_vector.cargue(0,-g,0);
  //Borrar todas las fuerzas y torques
  for(i=0;i<N+4;i++) Grano[i].BorreFuerzayTorque();
  //Agregue la fuerza de gravedad
  for(i=0;i<N;i++) Grano[i].AgregueFuerza(Grano[i].m*g_vector);
  //Calcular todas las fuerzas entre parejas de granos
  for(i=0;i<N;i++)
    for(j=i+1;j<N+4;j++)
      CalculeLaFuerzaEntre(Grano[i],Grano[j],ele[i][j],EstabaEnContacto[i][j],dt);
}
void Colisionador::CalculeLaFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2,
				       vector3D & ele, bool & EstabaEnContacto, double dt){
  vector3D r21,n,Vc,Vcn,Vct,t,Fn,Ft,F2; 
  double R1,R2,d21,s,m1,m2,m12,componenteVcn,componenteFn,normaVct,Ftmax,normaFt;
  double ERFF=1e-8;
  r21=Grano2.r-Grano1.r; d21=norma(r21);  s=(Grano1.R+Grano2.R)-d21;
  if(s>0){ //Si se chocan,
    //Geometría y dinámica del contacto
    m1=Grano1.m;   m2=Grano2.m;   m12=(m1*m2)/(m1+m2);
    R1=Grano1.R;   R2=Grano2.R;
    n=r21/d21;
   //Calcular velocidad de contacto y el vector tangente
    Vc=(Grano2.V-Grano1.V)-(Grano2.omega^n)*R2-(Grano1.omega^n)*R1;
    componenteVcn=Vc*n; Vcn=n*componenteVcn; Vct=Vc-Vcn;  normaVct=norma(Vct);
    if(normaVct<ERFF) t.cargue(0,0,0); else t=Vct/normaVct;

    //FUERZAS NORMALES
    //Fuerza de Hertz
    componenteFn=K*pow(s,1.5); 
    //Disipacion plástica
    componenteFn-=m12*sqrt(s)*Gamma*componenteVcn; if(componenteFn<0) componenteFn=0;
    Fn=n*componenteFn;

    //FUERZAS TANGENCIALES
    //fuerza estática
    ele+=(Vct*dt);
    Ft=ele*(-Kcundall);
    //fuerza cinética
    Ftmax=MU*componenteFn; normaFt=norma(Ft);
    if(normaFt>Ftmax) Ft=ele*(-Ftmax/norma(ele));

    //Construir la fuerza total
    F2=Fn+Ft;
    Grano2.AgregueFuerza(F2);      Grano2.AgregueTorque((n*(-R2))^Ft);      
    Grano1.AgregueFuerza(F2*(-1)); Grano1.AgregueTorque((n*R1)^(Ft*(-1))); 

    EstabaEnContacto=true;
  }
  else if(EstabaEnContacto==true){
    ele.cargue(0,0,0); EstabaEnContacto=false;
  }
    
}
//----------------------Funciones Globales---------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'MiGrano.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-10:110]"<<endl;
  cout<<"set yrange [-10:110]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<100/7<<"*t,0";
    cout<<" , "<<100/7<<"*t,100";
    cout<<" , 0,"<<100/7<<"*t";
    cout<<" , 100,"<<100/7<<"*t";
}
void TermineCuadro(void){
    cout<<endl;
}


//--------------------Programa Principal------------------
int main(void){
  double t, dt=1e-3;
  double tdibujo; int Ndibujos;
  Cuerpo Grano[N+4]; int i,j;
  Colisionador Newton;
  Crandom ran64(1); double theta;

  double m0=1,R0=6,V=10,omega0=10;
  double Rpared=10000,Mpared=1000; 

  double T=Lx/V, teq=10*T, tmax=2*T;

  double dx=Lx/(Nx+1), dy=Ly/(Ny+1);

  Ndibujos=2000;
  InicieAnimacion(); 
  
  //PAREDES
  //---------------(x0  ,       y0,z0,Vx0,Vy0,Vz0,theta0,omega0,    m0,     R0);
  //Pared arriba
  Grano[N  ].Inicie(Lx/2,Ly+Rpared, 0,  0,  0,  0,     0,     0,Mpared, Rpared);
  //Pared abajo
  Grano[N+1].Inicie(Lx/2,  -Rpared, 0,  0,  0,  0,     0,     0,Mpared, Rpared);
  //Pared derecha
  Grano[N+2].Inicie(Lx+Rpared,Ly/2, 0,  0,  0,  0,     0,     0,Mpared, Rpared);
  //Pared izquierda
  Grano[N+3].Inicie(  -Rpared,Ly/2, 0,  0,  0,  0,     0,     0,Mpared, Rpared);

  //GRANOS
  for(i=0;i<Nx;i++)
    for(j=0;j<Ny;j++){
      //------------------(x0      ,y0      ,z0,         Vx0,         Vy0,Vz0,theta0,omega0,m0,R0);
      theta=2*M_PI*ran64.r();
      Grano[i+Nx*j].Inicie((i+1)*dx,(j+1)*dy, 0,           0,           0,  0,     0,omega0,m0,R0);
    }
  
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    if(tdibujo>tmax/Ndibujos){
      //    if(t>teq)
      //      for(i=0;i<N;i++) cout<<Grano[i].GetVx()<<endl;
      
      InicieCuadro();
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();
      
      tdibujo=0;
    }
    
    //Muevase con Omelyan PEFRL
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,Zeta);
    Newton.CalculeTodasLasFuerzas(Grano,dt); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,(1-2*Lambda)/2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,Xi);
    Newton.CalculeTodasLasFuerzas(Grano,dt); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,Lambda);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,1-2*(Xi+Zeta));
    Newton.CalculeTodasLasFuerzas(Grano,dt); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,Lambda);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,Xi);
    Newton.CalculeTodasLasFuerzas(Grano,dt); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,(1-2*Lambda)/2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,Zeta);
  }
  
  return 0;
}
