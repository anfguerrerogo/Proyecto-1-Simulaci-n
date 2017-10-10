//Simulación Reloj de Arena
#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
#include "CuerpoColisionador.h"
using namespace std;
//---------------------Constantes globales--------
const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;

const double Lx=12,Ly=19;
const double Rpared=10000,Mpared=1000, Rparedcircular=0.48*Lx;//Propiedades paredes
const double DeltaLy=Lx*0.8;//Fracción de la caja donde comienzan las partículas
const double m0=1,R0=0.1;//Propiedades partículas

//Se requiere Nx<Lx/(2*R0), Ny<DeltaLy/(2*R0) para que no se toquen al inicio.
const double Espacioparticulasx=(Lx-2*Nx*R0)/(Nx+1), Espacioparticulasy=(DeltaLy-2*Ny*R0)/(Ny+1);//Espacio entre bordes
const double dx=2*R0 + Espacioparticulasx, dy=2*R0 + Espacioparticulasy;//Espacio entre centros

const double tmax=50;
//----------------------Funciones animación--------------
void InicieAnimacion(void);
void InicieCuadro(void);
void TermineCuadro(void);
void Dibujar(Cuerpo* Particula, double tdibujo, int Ndibujos);

//----------------------Funciones programa---------------
void MuevaParticulas(Cuerpo* Particula, Colisionador Newton, double dt);
void Flujo(Cuerpo* Particula, double & Numpart,  double & Flowrate,  double & Numpartabajo, double t);
void CondicionesIniciales(Cuerpo* Particula);

//----------------------Programa Principal---------------
int main(void){
  double t, dt=1e-4;
  Cuerpo Particula[N+5];
  Colisionador Newton;
  Crandom ran64(1); 
  double Numpartabajo=0, Flowrate=0, NumpartabajoOld=0, Numpart=0;
  int Tint=0;
  double tdibujo=0;
  double Pospreviasy[N];//Guarda posiciones previas en y
  int Ndibujos=1000;
  
  CondicionesIniciales(Particula);
   
  //InicieAnimacion();
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    //Dibujar(Particula, tdibujo, Ndibujos);
    for(int i=0;i<N;i++)Pospreviasy[i]=Particula[i].Gety();
    MuevaParticulas(Particula, Newton, dt);
    for(int i=0;i<N;i++){if(Particula[i].Gety()*Pospreviasy[i]<0) Numpart++ ;}
    Numpartabajo+=Numpart;
    Tint=(int) (t/dt);
    //cout<<t<<" "<<Numpart<<" "<<Numpartabajo<<" "<<Tint<<endl;
    if(Tint % 100 == 0) Flujo(Particula,  Numpart, Flowrate, Numpartabajo, t);
    for(int i=0;i<N;i++){
      if(Particula[i].Gety()<-1)Particula[i].Inicie((Lx-2*R0)*ran64.r(), Ly-2*R0, 0, 0, 0, 0, 0, 0, m0, R0);
    }
    if(Numpartabajo==2*N) break;
  }
  return 0;
}

//----------------------Funciones Globales---------------
void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl; 
  //cout<<"set output 'MiParticula.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-3:15]"<<endl;
  cout<<"set yrange [-1:21]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    //cout<<" , "<<100/7<<"*t,0";//Pared abajo
    cout<<" , "<<Lx/7<<"*t,19";//Pared arriba
    cout<<" , 0,"<<Ly/7<<"*t"; //Pared izquierda
    cout<<" , 12,"<<Ly/7<<"*t";//Pared derecha
    cout<<", "<<Rparedcircular<<"*cos(t/2),"<<Rparedcircular<<"*sin(t/2)";//Dibuja pared circular izquierda
    cout<<", "<<Lx<<"+"<<Rparedcircular<<"*cos(t/2),"<<Rparedcircular<<"*sin(t/2)";//Dibuja pared circular derecha
}
void TermineCuadro(void){
    cout<<endl;
}

void Dibujar(Cuerpo* Particula, double tdibujo, int Ndibujos){
  if(tdibujo>tmax/Ndibujos){
      InicieCuadro();
      for(int i=0;i<N;i++) Particula[i].Dibujese();
      TermineCuadro();     
      tdibujo=0;
      }}

void MuevaParticulas(Cuerpo* Particula,Colisionador Newton, double dt){
  int i;
  //Muevase con Omelyan PEFRL
  for(i=0;i<N;i++) Particula[i].Mueva_r(dt,Zeta);
  Newton.CalculeTodasLasFuerzas(Particula,dt); for(i=0;i<N;i++) Particula[i].Mueva_V(dt,(1-2*Lambda)/2);
  for(i=0;i<N;i++) Particula[i].Mueva_r(dt,Xi);
  Newton.CalculeTodasLasFuerzas(Particula,dt); for(i=0;i<N;i++) Particula[i].Mueva_V(dt,Lambda);
  for(i=0;i<N;i++) Particula[i].Mueva_r(dt,1-2*(Xi+Zeta));
  Newton.CalculeTodasLasFuerzas(Particula,dt); for(i=0;i<N;i++) Particula[i].Mueva_V(dt,Lambda);
  for(i=0;i<N;i++) Particula[i].Mueva_r(dt,Xi);
  Newton.CalculeTodasLasFuerzas(Particula,dt); for(i=0;i<N;i++) Particula[i].Mueva_V(dt,(1-2*Lambda)/2);
  for(i=0;i<N;i++) Particula[i].Mueva_r(dt,Zeta);
}
void Flujo(Cuerpo* Particula, double & Numpart, double & Flowrate, double & Numpartabajo, double t){
//Método 2 para sacar el Flow Rate
    Flowrate=Numpart;
    cout<<t<<" "<<Flowrate<<" "<<Numpartabajo<<endl;
    Numpart=0;
    Flowrate=0;
}

void CondicionesIniciales(Cuerpo* Particula){
  int i=0, j=0;
  //PAREDES
  //---------------(x0, y0, z0, Vx0,Vy0, Vz0, theta0, omega0, m0, R0);
  Particula[N  ].Inicie(Lx/2,Ly+Rpared, 0,  0,  0,  0,     0,     0, Mpared, Rpared,0);  //Pared arriba
  Particula[N+1].Inicie(Lx+Rpared,Ly/2, 0,  0,  0,  0,     0,     0, Mpared, Rpared,0);  //Pared derecha
  Particula[N+2].Inicie(  -Rpared,Ly/2, 0,  0,  0,  0,     0,     0, Mpared, Rpared,0);  //Pared izquierda
  Particula[N+3].Inicie(0,  0,  0,   0, 0,  0,     0,     0, Mpared, Rparedcircular,0);  //Pared circular izquierda
  Particula[N+4].Inicie(Lx, 0,  0,   0, 0,  0,     0,     0, Mpared, Rparedcircular,0);  //Pared circular derecha

  double carga=1;
  //GRANOS
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      //------------(x0, y0, z0, Vx0, Vy0, Vz0, theta0, omega0, m0, R0);
      Particula[i+Nx*j].Inicie(Espacioparticulasx+(i)*dx, Ly-Espacioparticulasy-(j)*dy, 0, 0, 0, 0, 0, 0, m0, R0, carga);
    }
  }
}

/*
//Método 2 para sacar flow rate
void Flujo(Cuerpo* Particula, double & Numpartabajo, double & NumpartabajoOld, double t, double & Flowrate){
  int i;
  NumpartabajoOld=Numpartabajo;
  Numpartabajo=0;
  for(i=0;i<N;i++) {if(Particula[i].Gety()<0) Numpartabajo++;}
  Flowrate=Numpartabajo-NumpartabajoOld;
  cout<<t<<" "<<Flowrate<<" "<<Numpartabajo<<" "<<NumpartabajoOld<<endl;
}*/
