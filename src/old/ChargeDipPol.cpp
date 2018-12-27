#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime> 

#include <iomanip>
#include <fstream>
#include <iostream>

using namespace std;
#include "ChargeDipPol.h"
#include "vecdefs.h"
#include "lapackblas.h"
#include "constants.h"

//
//  this is for setting up a molecule with interacting atomic polarizability
//  InvA is the matrix needed to compute the induced dipoles from the "external" E field
//

//////////////////////////////////////////////////////////////
//
//   InvA = (alpha^-1 - T)^-1 
//
//   this is the matrix that can be multiplied with the electric field at
//   the polarizable sites to yield the induced dipoles:
//
//   mu = InvA * E
//
//   dim(A) = 3*nSiyes x 3*nSites
//
void InvChargeDipolePol(int nSites, const double *R, double *InvA, int nFullerenes, const int *NoAtomsArray)
{

  BuildPolarizationMatrix(nSites, R, InvA, nFullerenes, NoAtomsArray);  
  InvertCDMatrix(4*nSites+nFullerenes, InvA);

 int n = nSites;
int n4plus1 = 4*nSites+nFullerenes; 
double alpha_xx = 0.0;
double alpha_yy = 0.0;
double alpha_zz = 0.0;
double alpha_xy = 0.0;
double alpha_yx = 0.0;
double alpha_xz = 0.0;
double alpha_zx = 0.0;
double alpha_yz = 0.0;
double alpha_zy = 0.0;
double temp=0.0;
/*

double Q_xx = 0.0;
double Q_yy = 0.0;
double Q_zz = 0.0;
double Q_xy = 0.0;
double Q_xz = 0.0;
double Q_yz = 0.0;

  for (int i = 0; i < n; ++i){
    double Ri2 = R[3*i+0]*R[3*i+0] + R[3*i+1]*R[3*i+1] + R[3*i+2]*R[3*i+2];
    Q_xx += 3.0*InvA[n4plus1*i]*R[3*i+0]*R[3*i+0] - Ri2;
    Q_yy += 3.0*InvA[n4plus1*i]*R[3*i+1]*R[3*i+1] - Ri2;
    Q_zz += 3.0*InvA[n4plus1*i]*R[3*i+2]*R[3*i+2] - Ri2;
    Q_xy += 3.0*InvA[n4plus1*i]*R[3*i+0]*R[3*i+1];
    Q_xz += 3.0*InvA[n4plus1*i]*R[3*i+0]*R[3*i+2];
    Q_yz += 3.0*InvA[n4plus1*i]*R[3*i+1]*R[3*i+2];
  }
*/

  for (int j = 0; j < n; ++j){
    for (int i = 0; i < n; ++i){

//For the charge-charge term

       alpha_xx += InvA[n4plus1*i+j]*R[3*i+0]*R[3*j+0];
       alpha_xy += InvA[n4plus1*i+j]*R[3*i+0]*R[3*j+1];
       alpha_xz += InvA[n4plus1*i+j]*R[3*i+0]*R[3*j+2];
       alpha_yx += InvA[n4plus1*i+j]*R[3*i+1]*R[3*j+0];
       alpha_yy += InvA[n4plus1*i+j]*R[3*i+1]*R[3*j+1];
       alpha_yz += InvA[n4plus1*i+j]*R[3*i+1]*R[3*j+2];
       alpha_zx += InvA[n4plus1*i+j]*R[3*i+2]*R[3*j+0];
       alpha_zy += InvA[n4plus1*i+j]*R[3*i+2]*R[3*j+1];
       alpha_zz += InvA[n4plus1*i+j]*R[3*i+2]*R[3*j+2];


// For the dipole-chrage term
       alpha_xx -= 2*InvA[n*n4plus1+j*3*n4plus1                  +i]*R[3*i+0]; 
       alpha_xy -= 2*InvA[n*n4plus1+j*3*n4plus1                  +i]*R[3*i+1]; 
       alpha_xz -= 2*InvA[n*n4plus1+j*3*n4plus1                  +i]*R[3*i+2]; 
       alpha_yx -= 2*InvA[n*n4plus1+j*3*n4plus1+ n4plus1         +i]*R[3*i+0]; 
       alpha_yy -= 2*InvA[n*n4plus1+j*3*n4plus1+ n4plus1         +i]*R[3*i+1]; 
       alpha_yz -= 2*InvA[n*n4plus1+j*3*n4plus1+ n4plus1         +i]*R[3*i+2]; 
       alpha_zx -= 2*InvA[n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i]*R[3*i+0]; 
       alpha_zy -= 2*InvA[n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i]*R[3*i+1]; 
       alpha_zz -= 2*InvA[n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i]*R[3*i+2]; 


//      cout<<"R["<<3*i+0<<"]="<<R[3*i+0]<<endl;
//      cout<<"R["<<3*i+1<<"]="<<R[3*i+1]<<endl;
//      cout<<"R["<<3*i+2<<"]="<<R[3*i+2]<<endl;
//      cout<<"InvA[i*n4plus1+n+3*j+0] =" <<InvA[i*n4plus1+n+3*j+0] <<endl;
//      cout<<"muti =" <<InvA[i*n4plus1+n+3*j+0] * R[3*i+0] <<endl;
//      cout<<"InvA[i*n4plus1+n+3*j+1] =" <<InvA[i*n4plus1+n+3*j+1] <<endl;
//      cout<<"InvA[i*n4plus1+n+3*j+2] =" <<InvA[i*n4plus1+n+3*j+2] <<endl;

//        cout<<"alpha_xx = "<<alpha_xx<<endl;


// For the dipole-dipole term
       alpha_xx += InvA[n*n4plus1+3*i*n4plus1 +                  n+j*3+0];
       alpha_xy += InvA[n*n4plus1+3*i*n4plus1 +                  n+j*3+1]; 
       alpha_xz += InvA[n*n4plus1+3*i*n4plus1 +                  n+j*3+2]; 
       alpha_yx += InvA[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+0]; 
       alpha_yy += InvA[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+1];
       alpha_yz += InvA[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+2]; 
       alpha_zx += InvA[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+0]; 
       alpha_zy += InvA[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+1]; 
       alpha_zz += InvA[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+2];

//       cout<<"2 alpha_xx "<<InvA[n*n4plus1+3*i*n4plus1 +                  n+j*3+0]<<endl;
 
/*
       temp += InvA[i*n4plus1+n+3*j+0];
       cout<<"temp = "<<temp<<" "<< InvA[i*n4plus1+n+3*j+0]<<endl;
       alpha_xx += InvA[i*n4plus1+n+3*j+0];
       alpha_yy += InvA[i*n4plus1+n+3*j+1];
       alpha_zz += InvA[i*n4plus1+n+3*j+2];
*/

    }
  } 
 cout<<"---Polarizability  -------"<<endl;
 cout<<"alpha_xx = "<<alpha_xx<<endl;
 cout<<"alpha_yy = "<<alpha_yy<<endl;
 cout<<"alpha_zz = "<<alpha_zz<<endl;
 cout<<"alpha_xy = "<<alpha_xy<<endl;
 cout<<"alpha_yx = "<<alpha_yx<<endl;
 cout<<"alpha_xz = "<<alpha_xz<<endl;
 cout<<"alpha_zx = "<<alpha_zx<<endl;
 cout<<"alpha_yz = "<<alpha_yz<<endl;
 cout<<"alpha_zy = "<<alpha_zy<<endl;

/*
 cout<<"---Quadrupole moment  -------"<<endl;
 cout<<"Q_xx = "<<Q_xx<<endl;
 cout<<"Q_yy = "<<Q_yy<<endl;
 cout<<"Q_zz = "<<Q_zz<<endl;
 cout<<"Q_xy = "<<Q_xy<<endl;
 cout<<"Q_xz = "<<Q_xz<<endl;
 cout<<"Q_yz = "<<Q_yz<<endl;
*/

/*
  for (int i = 0; i < n4plus1*n4plus1; ++i){
     cout<<" InvA["<<i<<"]= "<< InvA[i]<<endl;;
  }
*/

}

//////////////////////////////////////////////////////////////
//
//  build matrix A = (alpha^-1 - T)
//
//  this is a super matrix consisting of:
//  
//  T(qq) T(qd) T(1)
//  T(qd) T(dd) T(1)
//  T(1)  T(1)  0 
//
void BuildPolarizationMatrix2(int nSites, const double *R, double *A, int nFullerenes, const int *NoAtomsArray)
{

  int n = nSites;
  int n3 = n*3 ;

  for (int k = 0; k < n3*n3 ; ++k)
    A[k] = 0.0;
  
//  double R_const = 0.62*1.889725989 ; 
  double R_const = 0.68620399*1.889725989 ; 
  double chi = sqrt(2.0/PI)/R_const  ;  
  double Rqq = R_const*sqrt(2.0)  ;
  //cout << "R_const , R_qq" << R_const << " "<<Rqq << endl ; 
  //T-lagrange
  int current_atom = 0 ; 
  int old_atom     = 0 ;    
  
  double alpha = 3.0*sqrt(PI/2)*R_const*R_const*R_const ;
  //Tpq and Tpp
  for (int j = 0; j < n; ++j){
    for (int i = 0; i < n; ++i){
      if(i==j){
        A[3*i*n3 +        j*3+0] = 1.0/alpha ;
        A[3*i*n3 +n3+     j*3+1] = 1.0/alpha ;
        A[3*i*n3 +n3+ n3+ j*3+2] = 1.0/alpha ;
      }
      else {
        double rij[3] ;
        rij[0] = R[3*i+0] - R[3*j+0];
        rij[1] = R[3*i+1] - R[3*j+1];
        rij[2] = R[3*i+2] - R[3*j+2];
        double Rij = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
        double fac, fac1, fac2, fac_pq ; 
        double Txy, Txz, Tyz ;
        fac = erf(Rij/Rqq)-2.0*Rij*exp(-1.0*(Rij/Rqq)*(Rij/Rqq))/(sqrt(PI)*Rqq) ;
        fac_pq = fac/(Rij*Rij*Rij) ; 
        fac1 = fac_pq/(Rij*Rij) ; 
        fac2 = 4.0*exp(-1.0*(Rij/Rqq)*(Rij/Rqq))/(sqrt(PI)*Rqq*Rqq*Rqq*Rij*Rij) ;
        Txy = 3.0*rij[0]*rij[1]*fac1 - rij[0]*rij[1]*fac2 ;
        Txz = 3.0*rij[0]*rij[2]*fac1 - rij[0]*rij[2]*fac2 ;
        Tyz = 3.0*rij[1]*rij[2]*fac1 - rij[1]*rij[2]*fac2 ;
        A[3*i*n3 +        j*3+0] = -1.0*((3.0*rij[0]*rij[0]-Rij*Rij)*fac1 - rij[0]*rij[0]*fac2) ; 
        A[3*i*n3 +n3+     j*3+1] = -1.0*((3.0*rij[1]*rij[1]-Rij*Rij)*fac1 - rij[1]*rij[1]*fac2) ; 
        A[3*i*n3 +n3+ n3+ j*3+2] = -1.0*((3.0*rij[2]*rij[2]-Rij*Rij)*fac1 - rij[2]*rij[2]*fac2) ; 
        A[3*i*n3 +        j*3+1] = -1.0*Txy ; 
        A[3*i*n3 +n3+     j*3+0] = -1.0*Txy ; 
        A[3*i*n3 +        j*3+2] = -1.0*Txz ; 
        A[3*i*n3 +n3+ n3+ j*3+0] = -1.0*Txz ; 
        A[3*i*n3 +n3+     j*3+2] = -1.0*Tyz ; 
        A[3*i*n3 +n3+ n3+ j*3+1] = -1.0*Tyz ;
         
        //Tpq

      }
     
    }
  }


}


//////////////////////////////////////////////////////////////
//
//  build matrix A = (alpha^-1 - T)
//
//  this is a super matrix consisting of:
//  
//  T(qq) T(qd) T(1)
//  T(qd) T(dd) T(1)
//  T(1)  T(1)  0 
//
void BuildPolarizationMatrix(int nSites, const double *R, double *A, int nFullerenes, const int *NoAtomsArray)
{

  int n = nSites;
  int n4 = n*4 ;
  int n4plus1 = 4*n + nFullerenes ; // dimension of super-matrix A and leading dimension of any block B
  int nplus1 = n + nFullerenes ; 

   cout<<"n4plus1 = "<<n4plus1<<endl;

  for (int k = 0; k < n4plus1*n4plus1 ; ++k)
    A[k] = 0.0;
  
  double* Matrix= new double[nplus1*nplus1];
  double ChargePen=-1.5;
// double R_const = 0.62*1.889725989 ; 
  double R_const = 0.68620399*1.889725989 ; 
 // double R_const = 2; 
  double chi = sqrt(2.0/PI)/R_const  ;  
  double Rqq = R_const*sqrt(2.0)  ;
  //cout << "R_const , R_qq" << R_const << " "<<Rqq << endl ; 
  //Tqq
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      if (i == j) {
        A[n4plus1*i+j] = chi; 
        Matrix[nplus1*i+j] = chi;
       // A[n4plus1*i+j] = 0.0; 
      //  cout << "A" << n4plus1*i+j << "= " << A[n4plus1*i+j] << endl; 
      }
      else {
        double rij[3] ;
        rij[0] = R[3*i+0] - R[3*j+0];
        rij[1] = R[3*i+1] - R[3*j+1];
        rij[2] = R[3*i+2] - R[3*j+2];
        double Rij = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
        double CT=4.0*pow(10,4)*exp(-1.5*Rij);
/*

        if ((i < 60 && j >= 60) || (i >= 60 && j < 60) )  {
          // A[n4plus1*i+j] = erf(Rij/Rqq)/(Rij*CT) ; 
          if ( CT <= 1.0  )   A[n4plus1*i+j] = erf(Rij/Rqq)/(Rij*CT) ; 
          else A[n4plus1*i+j] = erf(Rij/Rqq)/Rij ;
        //   cout<<"i and j and CT Rij :  "<<i<<", "<<j<<" "<<CT<<" " <<Rij<<endl;
       //   A[n4plus1*i+j] = 0.0 ;
        }
        else
*/

/*
        double C = 1.00;
        double B = 0.43; 
        double w = Rij*B;
        if ((i < 60 && j >= 60) || (i >= 60 && j < 60) )  {
         A[n4plus1*i+j] = erf(Rij/Rqq)/Rij *  C*exp(-w)*(1+w+1.0/3.0*w*w) ; 
        }
        else  A[n4plus1*i+j] = erf(Rij/Rqq)/Rij ; 
*/
         A[n4plus1*i+j] = erf(Rij/Rqq)/Rij ; 
         Matrix[nplus1*i+j] = erf(Rij/Rqq)/Rij ;
         // A[n4plus1*i+j] = erf(Rij/Rqq)/(Rij*CT) ; 
         // cout<<" erf(Rij/Rqq)/(Rij*CT) =  "<< erf(Rij/Rqq)/(Rij*CT) <<endl; 


       // A[n4plus1*i+j] = 0.0 ; 
         
       // cout<<"Rij ="<<Rij<<" Rqq = "<<Rqq<<" erf(Rij/Rqq)= "<<erf(Rij/Rqq)<<endl;
      //  cout << "A" << n4plus1*i+j << "= " << A[n4plus1*i+j] << endl; 
      }
    }
  }
/*
  //Tlagrangei
  cout <<"  warning, no restricted charge flow " << endl ; 
  for (int i = 0; i < n; ++i){
    A[i*n4plus1+n4plus1-1] = 1.0 ;
    A[i+n4plus1*(n4plus1-1)] = 1.0 ;
  //    cout << "A" << i*n4plus1+n4plus1-1 << "= " << A[i*n4plus1+n4plus1-1] << endl;
  //    cout << "A" << i+n4plus1*(n4plus1-1) << "= " << A[i+n4plus1*(n4plus1-1)] << endl; 
  }
  A[n4plus1*n4plus1 -1]= 0.0 ;
*/     

  //T-lagrange
  int current_atom = 0 ; 
  int old_atom     = 0 ;    
  //double Aconst=-341.; // For 8.0 Angstrom of dimer
  // double Aconst=-11.85; // For 9.0 Angstrom of dimer
  // double Aconst=-4.02; // For 9.25 Angstrom of dimer
  // double Aconst=-2.08; // For 9.50 Angstrom of dimer
  // double Aconst=-1.03; // For 9.75 Angstrom of dimer
 //  double Aconst=-0.525; // For 10.05 Angstgrom of dimer
  // double Aconst=-0.345; // For 10.25 Angstgrom of dimer
  //double Aconst=-0.21; // For 10.50 Angstgrom of dimer
 // double Aconst=-0.128; // For 10.75 Angstgrom of dimer
 // double Aconst=-0.079; // For 11.0 Angstrom of dimer 
  //double Aconst=-0.0091; // For 12.0 Angstrom of dimer 
  //double Aconst=-0.0003; // For 13.0 Angstrom of dimer 
  double Aconst=-0.00000001; // For 14.0 Angstrom of dimer 
 // double Aconst= 0.0000000;
  //double Bconst= -11.85; 
  double Bconst= 0.000000; 
  for ( int iFuller = 0 ; iFuller < nFullerenes ; ++iFuller){
    old_atom      = current_atom ;
    current_atom += NoAtomsArray[iFuller] ;
    for (int i = old_atom; i < current_atom ; ++i) {
      A[i*n4plus1+n4 + iFuller]               = 1.0 ; // corresponds to columns of lagrange vector
      A[n4plus1*n4 + n4plus1*iFuller + i]     = 1.0 ; // corresponds to rows of lagrange rows
      Matrix[i*nplus1+n + iFuller]               = 1.0 ; // corresponds to columns of lagrange vector
      Matrix[nplus1*n + nplus1*iFuller + i]     = 1.0 ; // corresponds to rows of lagrange rows


/*
     if (iFuller ==0) {
      A[i*n4plus1+n4 + iFuller]               = 1.0 ; // corresponds to columns of lagrange vector
      A[i*n4plus1+n4 + iFuller+1]               = 1.0 ; // corresponds to columns of lagrange vector
      A[n4plus1*n4 + n4plus1*iFuller + i]     = 1.0 ; // corresponds to rows of lagrange rows
      A[n4plus1*n4 + n4plus1*iFuller + i + n4plus1]     = 1.0 ; // corresponds to rows of lagrange rows
      Matrix[i*nplus1+n + iFuller]               = 1.0 ; // corresponds to columns of lagrange vector
      Matrix[i*nplus1+n + iFuller+1]               = 1.0 ; // corresponds to columns of lagrange vector
      Matrix[nplus1*n + nplus1*iFuller + i]     = 1.0 ; // corresponds to rows of lagrange rows
      Matrix[nplus1*n + nplus1*iFuller + i + nplus1]     = 1.0 ; // corresponds to rows of lagrange rows
    }
      if (iFuller ==1) { 
      A[i*n4plus1+n4 + iFuller]               = -1.0 ; // corresponds to columns of lagrange vector
      A[i*n4plus1+n4 + iFuller-1]               = 1.0 ; // corresponds to columns of lagrange vector
      A[n4plus1*n4 + n4plus1*iFuller + i]     = -1.0 ; // corresponds to rows of lagrange rows
      A[n4plus1*n4 + n4plus1*iFuller + i - n4plus1]     = 1.0 ; // corresponds to rows of lagrange rows
      Matrix[i*nplus1+n + iFuller]               = -1.0 ; // corresponds to columns of lagrange vector
      Matrix[i*nplus1+n + iFuller-1]               = 1.0 ; // corresponds to columns of lagrange vector
      Matrix[nplus1*n + nplus1*iFuller + i]     = -1.0 ; // corresponds to rows of lagrange rows
      Matrix[nplus1*n + nplus1*iFuller + i - nplus1]     = 1.0 ; // corresponds to rows of lagrange rows
     }
*/
     // cout<< "A[ "<<i*n4plus1+n4 + iFuller<<"]= "<<A[i*n4plus1+n4 + iFuller]<<endl; 
    //  cout<< "A[ "<<i*n4plus1+n4 + iFuller-1<<"]= "<<A[i*n4plus1+n4 + iFuller-1]<<endl; 
     // cout<< "A[ "<<n4plus1*n4 + n4plus1*iFuller + i<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + i]<<endl; 
     // A[i*n4plus1+n4 + iFuller]               = 0.0 ; // corresponds to columns of lagrange vector
     //  A[n4plus1*n4 + n4plus1*iFuller + i]     = 0.0 ; // corresponds to rows of lagrange rows
    }
    if (iFuller == 0 ) {
//       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-3]     = Aconst; // corresponds to rows of lagrange rows
//       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]     = -Aconst; // corresponds to rows of lagrange rows
//       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]     = 0.0 ; // corresponds to rows of lagrange rows

       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]     = Aconst; // corresponds to rows of lagrange rows
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]     = -Aconst; // corresponds to rows of lagrange rows
    cout<< "0 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-3<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]<<endl; 
    cout<< "0 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-2<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]<<endl; 
    cout<< "0 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-1<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]<<endl; 
    }
    else if (iFuller == 1 ) {
      // A[n4plus1*n4 + n4plus1*iFuller + n4plus1-3]     =  -Aconst; // corresponds to rows of lagrange rows
     //  A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]     = Aconst+Bconst; // corresponds to rows of lagrange rows
     //  A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]     = -Bconst ; // corresponds to rows of lagrange rows
     
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]     = -Aconst; // corresponds to rows of lagrange rows
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]     = Aconst ; // corresponds to rows of lagrange rows
   //    A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]     = 4.0*Aconst ; // corresponds to rows of lagrange rows
    //   Matrix[nplus1*n + nplus1*iFuller + nplus1-1]     = 4.0*Aconst ; // corresponds to rows of lagrange rows

       cout<< " 1 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-3<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-3]<<endl; 
       cout<< " 1 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-2<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]<<endl; 
       cout<< " 1 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-1<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]<<endl; 
    }
    else if (iFuller == 2 ) {
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-3]     =  0.0; // corresponds to rows of lagrange rows
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]     = -Bconst; // corresponds to rows of lagrange rows
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]     = Bconst ; // corresponds to rows of lagrange rows
       cout<< " 2 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-3<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-3]<<endl; 
       cout<< " 2 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-2<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]<<endl; 
       cout<< " 2 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-1<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]<<endl; 
    }

  }
  
  double alpha = 3.0*sqrt(PI/2)*R_const*R_const*R_const ;
  cout <<"alpha : "<<alpha<<endl;
  //Tpq and Tpp
  for (int j = 0; j < n; ++j){
    for (int i = 0; i < n; ++i){
      if(i==j){
        A[n*n4plus1+3*i*n4plus1 +                  n+j*3+0] = 1.0/alpha ;
        A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+1] = 1.0/alpha ;
        A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+2] = 1.0/alpha ;
       // cout<<" i=j 1 : "<<n*n4plus1+3*i*n4plus1 +                  n+j*3+0 <<endl;
      //  cout<<" i=j 2 : "<<n*n4plus1+3*i*n4plus1 + n4plus1+         n+j*3+1 <<endl;
      //  cout<<" i=j 3 : "<<n*n4plus1+3*i*n4plus1 + n4plus1+ n4plus1+n+j*3+2 <<endl;
      }
      else {
        double rij[3] ;
        rij[0] = R[3*i+0] - R[3*j+0];
        rij[1] = R[3*i+1] - R[3*j+1];
        rij[2] = R[3*i+2] - R[3*j+2];
        double Rij = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
        double fac, fac1, fac2, fac_pq ; 
        double Txy, Txz, Tyz ;
        fac = erf(Rij/Rqq)-2.0*Rij*exp(-1.0*(Rij/Rqq)*(Rij/Rqq))/(sqrt(PI)*Rqq) ;
        fac_pq = fac/(Rij*Rij*Rij) ; 
       // cout<<"Rij = "<<Rij<<endl;
       // cout<<"Rqq = "<<Rqq<<endl;
       // cout<<"fac = "<<fac<<endl;
        fac1 = fac_pq/(Rij*Rij) ; 
        fac2 = 4.0*exp(-1.0*(Rij/Rqq)*(Rij/Rqq))/(sqrt(PI)*Rqq*Rqq*Rqq*Rij*Rij) ;
 
        Txy = 3.0*rij[0]*rij[1]*fac1 - rij[0]*rij[1]*fac2 ;
        Txz = 3.0*rij[0]*rij[2]*fac1 - rij[0]*rij[2]*fac2 ;
        Tyz = 3.0*rij[1]*rij[2]*fac1 - rij[1]*rij[2]*fac2 ;






        A[n*n4plus1+3*i*n4plus1 +                  n+j*3+0] = -1.0*((3.0*rij[0]*rij[0]-Rij*Rij)*fac1 - rij[0]*rij[0]*fac2) ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+1] = -1.0*((3.0*rij[1]*rij[1]-Rij*Rij)*fac1 - rij[1]*rij[1]*fac2) ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+2] = -1.0*((3.0*rij[2]*rij[2]-Rij*Rij)*fac1 - rij[2]*rij[2]*fac2) ; 


        A[n*n4plus1+3*i*n4plus1 +                  n+j*3+1] = -1.0*Txy ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+0] = -1.0*Txy ; 
        A[n*n4plus1+3*i*n4plus1 +                  n+j*3+2] = -1.0*Txz ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+0] = -1.0*Txz ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+2] = -1.0*Tyz ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+1] = -1.0*Tyz ;

         
/*
        cout<<" Txy: "<<n*n4plus1+3*i*n4plus1 +                  n+j*3+1<<endl; 
        cout<<" Txy: "<<n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+0<<endl; 
        cout<<" Txz: "<<n*n4plus1+3*i*n4plus1 +                  n+j*3+2<<endl; 
        cout<<" Txz: "<<n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+0<<endl; 
        cout<<" Tyz: "<<n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+2<<endl; 
        cout<<" Tyz: "<<n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+1<<endl;
*/



        //Tpq
       
/*
        
        A[i*n4plus1+n+3*j+0] = 0;
        A[i*n4plus1+n+3*j+1] = 0;
        A[i*n4plus1+n+3*j+2] = 0;


        A[n*n4plus1+i*3*n4plus1                  +j] = 0;
        A[n*n4plus1+i*3*n4plus1+ n4plus1         +j] = 0;
        A[n*n4plus1+i*3*n4plus1+ n4plus1+ n4plus1+j] = 0;
*/

       


//        if ((i < 60 && j >= 60) || (i >= 60 && j < 60) )  {
//        A[i*n4plus1+n+3*j+0] = -1.0*rij[0]*fac_pq*(1.0 - exp(ChargePen*Rij));
//        A[i*n4plus1+n+3*j+1] = -1.0*rij[1]*fac_pq*(1.0 - exp(ChargePen*Rij));
//        A[i*n4plus1+n+3*j+2] = -1.0*rij[2]*fac_pq*(1.0 - exp(ChargePen*Rij));
//        }
//        else {
        A[i*n4plus1+n+3*j+0] = -1.0*rij[0]*fac_pq;
        A[i*n4plus1+n+3*j+1] = -1.0*rij[1]*fac_pq;
        A[i*n4plus1+n+3*j+2] = -1.0*rij[2]*fac_pq;
//        }


/*

        A[n*n4plus1+i*3*n4plus1                  +j] = rij[0]*fac_pq*-1.0;
        A[n*n4plus1+i*3*n4plus1+ n4plus1         +j] = rij[1]*fac_pq*-1.0;
        A[n*n4plus1+i*3*n4plus1+ n4plus1+ n4plus1+j] = rij[2]*fac_pq*-1.0;


*/


//       if ((i < 60 && j >= 60) || (i >= 60 && j < 60) )  {
//        A[n*n4plus1+j*3*n4plus1                  +i] = -1.0*rij[0]*fac_pq*(1.0 - exp(ChargePen*Rij));
//        A[n*n4plus1+j*3*n4plus1+ n4plus1         +i] = -1.0*rij[1]*fac_pq*(1.0 - exp(ChargePen*Rij));
//        A[n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i] = -1.0*rij[2]*fac_pq*(1.0 - exp(ChargePen*Rij));
//        }
//        else {
        A[n*n4plus1+j*3*n4plus1                  +i] = -1.0*rij[0]*fac_pq;
        A[n*n4plus1+j*3*n4plus1+ n4plus1         +i] = -1.0*rij[1]*fac_pq;
        A[n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i] = -1.0*rij[2]*fac_pq;
//        }



/*

        cout<<"first A["<<i*n4plus1+n+3*j+0<<"] = "<<rij[0]*fac_pq*-1.0<<endl;
        cout<<"first A["<<i*n4plus1+n+3*j+1<<"] = "<<rij[1]*fac_pq*-1.0<<endl;
        cout<<"first A["<<i*n4plus1+n+3*j+2<<"] = "<<rij[2]*fac_pq*-1.0<<endl;
        cout<<"second A["<<n*n4plus1+j*3*n4plus1                  +i<<"]="<<rij[0]*fac_pq*-1.0<<endl;
        cout<<"second A["<<n*n4plus1+j*3*n4plus1+ n4plus1         +i<<"]="<<rij[1]*fac_pq*-1.0<<endl;
        cout<<"second A["<<n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i<<"]="<<rij[2]*fac_pq*-1.0<<endl;

        cout<<"second A["<<n*n4plus1+i*3*n4plus1                  +j<<"]="<<rij[0]*fac_pq*-1.0<<endl;
        cout<<"second A["<<n*n4plus1+i*3*n4plus1+ n4plus1         +j<<"]="<<rij[1]*fac_pq*-1.0<<endl;
        cout<<"second A["<<n*n4plus1+i*3*n4plus1+ n4plus1+ n4plus1+j<<"]="<<rij[2]*fac_pq*-1.0<<endl;
*/
/*
        cout<<"Tqx("<<i<<","<<j<<": "<<i*n4plus1+n+3*j+0 <<" :  "<<  rij[0]*fac_pq*-1.0 <<endl;
        cout<<"Tqy : "<<i*n4plus1+n+3*j+1 <<" :  "<<  rij[1]*fac_pq*-1.0 <<endl;
        cout<<"Tqz : "<<i*n4plus1+n+3*j+2 <<" :  "<<  rij[2]*fac_pq*-1.0 <<endl;
        cout<<"Txq : "<<n*n4plus1+j*3*n4plus1                  +i<<" :  "<<  rij[0]*fac_pq*-1.0 <<endl;
        cout<<"Tyq : "<<n*n4plus1+j*3*n4plus1+ n4plus1         +i<<" :  "<<  rij[1]*fac_pq*-1.0 <<endl;
        cout<<"Tzq : "<<n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i<<" :  "<<  rij[2]*fac_pq*-1.0 <<endl;
*/ 
      }
     
    }
  }

/*
  for (int i = 0; i < n4plus1*n4plus1; ++i){
     cout<<" A["<<i<<"]= "<< A[i]<<endl;;
  }
*/

// Callculate the Hessian matrix diagonalization
  char jobz = 'V';
  char uplo = 'L';
  double* eigenvectors = new double[nplus1*nplus1];
  double* eigenvalues = new double[nplus1];

/*
  for (int i = 0; i < nplus1; ++i){
    cout<<i<<" : "; 
    for (int j = 0; j < nplus1; ++j){
     cout<< Matrix[i*nplus1+j]<<" ";
    }
    cout<<endl;
  }
*/
/*
  for (int i = 0; i < nplus1; ++i)
    memcpy(eigenvectors+i*(n+1), Matrix+(i+1)*(2*n-i)/2-(n-i),(n-i)*sizeof(double));
  int lwork = 12 * n;
  double* work = new double[lwork];
  int info=0;
  dsyev(&jobz, &uplo, &n, eigenvectors, &n, eigenvalues, work, &lwork, &info);

  for (int i = 0; i < nplus1; ++i){
     cout<<" eigenvalues["<<i<<"]= "<< eigenvalues[i]<<endl;;
  }

  if (info != 0) {
    cout << "DMDCBase::CallDSYEV: DSYEV returned with info = " << info << endl;
    exit(1);
  }
  delete[] work;
*/
  delete[] eigenvalues;
  delete[] eigenvectors;
  delete[] Matrix;






/*
  for (int i = 0; i < n4plus1*n4plus1; ++i){
     A[i]=0.0;
   }
*/

//  int icount = 0;
//   for (int i = 0; i < n4plus1; ++i)
//    for (int j = 0; j < n4plus1; ++j) {
//     cout<<"symmetry check ["<<icount<<"], ("<<i<<","<<j<<")=" << " : "<< n4plus1*i+j<<","<<i+n4plus1*j<<" : "<<A[n4plus1*i+j]- A[i+n4plus1*j]<<endl;
//      icount++;
//    }
}

/////////////////////////////////////////////////
//
//  symmetric matrix A is inverted in place
//
void InvertCDMatrix(int n, double *A)
{

  int N = n;
  char uplo = 'U';
  static iVec ipiv; ipiv.resize(n);    // this is called repeatedly to invert 9x9 matrices, but only during setup, so
  static dVec work; work.resize(n*n);  // static vectors should be fine here as long as nobody calls this in mutiple threads
  int lwork = n*n;
  int info;

  dgetrf(&N, &N, A, &n, &ipiv[0], &info);
  if (info != 0) {
    cout << "InverseMatrix, Error in dgetrf = " << info << "\n";
    exit(1);
  }

  dgetri(&N, A, &n, &ipiv[0], &work[0],  &lwork, &info);
  if (info != 0) {
    cout << "InverseMatrix, Error in dgetrfi = " << info << "\n";
    exit(1);
  }

//
//  for (int i = 0; i < n; ++i)
//    for (int j = i+1; j < n; ++j)
//     if (abs(A[n*i+j]- A[i+n*j]) >= 0.00001 ) cout<<"symmetry check =" << A[n*i+j]- A[i+n*j]<<endl;

}


