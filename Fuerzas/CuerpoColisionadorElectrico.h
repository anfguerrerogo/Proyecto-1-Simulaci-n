#include <iostream>
#include <cmath>
using namespace std;

//-----------------------Constantes globales--------------
const int g=9.8;
const double K=1e4, Gamma=50, Kcundall=10, MU=0.4;
const int Nx=5,Ny=3,N=Nx*Ny;   //Número partículas
//-------------------------Clase Cuerpo-------------------

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
  double Getx(void){return r.x();};
  double Gety(void){return r.y();}; 
  double GetV(void){return norma(V);};
  double GetVx(void){return V.x();};
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0, double y0, double z0,
	      double Vx0, double Vy0, double Vz0,
	      double theta0, double omega0, double m0, double R0){
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

//-----------------------Clase Colisionador-------------------------------------
class Colisionador{
private:
  vector3D ele[N+5][N+5]; bool EstabaEnContacto[N+5][N+5];
public:
  void Inicie(void);
  void CalculeTodasLasFuerzas(Cuerpo* Particula,double dt);
  void CalculeLaFuerzaEntre(Cuerpo & Particula1,Cuerpo & Particula2,
			    vector3D & ele, bool & EstabaEnContacto,double dt);
};
void Colisionador::Inicie(void){
  int i,j;
  for(i=0;i<N;i++)
    for(j=i+1;j<N+5;j++){
      ele[i][j].cargue(0,0,0); EstabaEnContacto[i][j]=false;
    }
}
void Colisionador::CalculeTodasLasFuerzas(Cuerpo* Particula,double dt){
  int i,j;
  vector3D g_vector; g_vector.cargue(0,-g,0);
  //Borrar todas las fuerzas y torques
  for(i=0;i<N+5;i++) Particula[i].BorreFuerzayTorque();
  //Agregue la fuerza de gravedad
  for(i=0;i<N;i++) Particula[i].AgregueFuerza(Particula[i].m*g_vector);
  //Calcular todas las fuerzas entre parejas de granos
  for(i=0;i<N;i++)
    for(j=i+1;j<N+5;j++)
      CalculeLaFuerzaEntre(Particula[i],Particula[j],ele[i][j],EstabaEnContacto[i][j],dt);
}
void Colisionador::CalculeLaFuerzaEntre(Cuerpo & Particula1,Cuerpo & Particula2,
				       vector3D & ele, bool & EstabaEnContacto, double dt){
  vector3D r21,n,Vc,Vcn,Vct,t,Fn,Ft,F2; 
  double R1,R2,d21,s,m1,m2,m12,componenteVcn,componenteFn,normaVct,Ftmax,normaFt;
  double ERFF=1e-8;
  r21=Particula2.r-Particula1.r; d21=norma(r21);  s=(Particula1.R+Particula2.R)-d21;
  if(s>0){ //Si se chocan,
    //Geometría y dinámica del contacto
    m1=Particula1.m;   m2=Particula2.m;   m12=(m1*m2)/(m1+m2);
    R1=Particula1.R;   R2=Particula2.R;
    n=r21/d21;
   //Calcular velocidad de contacto y el vector tangente
    Vc=(Particula2.V-Particula1.V)-(Particula2.omega^n)*R2-(Particula1.omega^n)*R1;
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
    Particula2.AgregueFuerza(F2);      Particula2.AgregueTorque((n*(-R2))^Ft);      
    Particula1.AgregueFuerza(F2*(-1)); Particula1.AgregueTorque((n*R1)^(Ft*(-1))); 

    EstabaEnContacto=true;
  }
  else if(EstabaEnContacto==true){
    ele.cargue(0,0,0); EstabaEnContacto=false;
  }
    
}
