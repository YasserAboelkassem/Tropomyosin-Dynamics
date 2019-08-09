//--------------------------------------------------------------------------------------//
//                     University of California San Diego                               //
//                         Dept. of Bioengineering                                      //
//                                                                                      //
//                                                                                      //
//--------------------------------------------------------------------------------------//
// Author: Yasser Aboelkassem, PhD                                                      // 
// Assistant Project Scientist at UCSD                                                  //
// yaboelkassem@ucsd.edu                                                                //
//--------------------------------------------------------------------------------------//
// Year  :  2017-2019                                                                   //
//-----------------------------                                                         //
//                          This code uses models the                                   //
//      Tropomyosin (B-C-M) motion as described by a multiwell energy potential         //
//--------------------------------------------------------------------------------------//
//                    |                                       |                         //
//                    |                                       |                         //
//                    |---------New Myofilament Model---------|                         //
//                    |                                       |                         //
//                    |                                       |                         //
//             |Stochastic ODE  to model thin filament activation|                      //
//--------------------------------------------------------------------------------------//
//                                                                                      //
// Copyright:                                                                           // 
//------------                                                                          //
// 1- This code is a part of a larger sofware that I am is underdevelopment to study    // 
//    Tm-Actin-Tn interaction dynamics during thin filament activation.                 //
// 2- You can use this code only to generate the results in our publication.            //
// 3- A complete sofware based on this code will be published seperately and will be    //
//    abailable for the community.
//--------------------------------------------------------------------------------------//
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <time.h>

const long int N_steps  = 3e2;           //  Steps for both 1st and 2dn parts of the sol.
const double   pi       = 4*atan(1);     //  Constant

double A, K_stiff;
double phi_B, phi_C, phi_M, phi_1, phi_2;
double phi_1_min, phi_1_max, phi_2_min, phi_2_max;
double T1, T2, T3, T4;
double dt,lamda, eta, KbT;
double W_noise;
double pCa, pCa50, Ca, Ca50, nH;
double Alpha1, Alpha2, Alpha3, Alpha4;

int main()
{
    vector<double> phi_RU02 = vector<double>(N_steps);
    vector<double> phi_RU03 = vector<double>(N_steps);
    vector<double> phi_RU04 = vector<double>(N_steps);
    vector<double> phi_RU05 = vector<double>(N_steps);
    vector<double> phi_RU06 = vector<double>(N_steps);
    vector<double> phi_RU07 = vector<double>(N_steps);
    vector<double> phi_RU08 = vector<double>(N_steps);
    vector<double> phi_RU09 = vector<double>(N_steps);
    vector<double> phi_RU10 = vector<double>(N_steps);
    vector<double> phi_RU11 = vector<double>(N_steps);
    vector<double> phi_RU12 = vector<double>(N_steps);
    vector<double> phi_RU13 = vector<double>(N_steps);
    vector<double> phi_RU14 = vector<double>(N_steps);
    vector<double> phi_RU15 = vector<double>(N_steps);
    vector<double> phi_RU16 = vector<double>(N_steps);
    vector<double> phi_RU17 = vector<double>(N_steps);
    vector<double> phi_RU18 = vector<double>(N_steps);
    vector<double> phi_RU19 = vector<double>(N_steps);
    vector<double> phi_RU20 = vector<double>(N_steps);
    vector<double> phi_RU21 = vector<double>(N_steps);
    vector<double> phi_RU22 = vector<double>(N_steps);
    vector<double> phi_RU23 = vector<double>(N_steps);
    vector<double> phi_RU24 = vector<double>(N_steps);
    vector<double> phi_RU25 = vector<double>(N_steps);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    normal_distribution<double> distribution (0.0,10.0);

    //--------------------------------
    A       = 1;                        // potential amplitude [Pn.nm]
    K_stiff = 1;                        // Stiff [pn.nm/rad]
    phi_B = 00.0 * pi/180;              // Blocked state "B" angle [rad]
    phi_C = 25.0 * pi/180;              // Closed  state "C" angle [rad]
    phi_M = 35.0 * pi/180;              // Open   state "M" angle  [rad]
    phi_1 = 27.0 * pi/180;              // This will be varied with Ca++ concentration [rad]

    //----------------------------
    Alpha4 =  4*1100;
    Alpha3 = -3*10100;
    Alpha2 =  2*310;
    Alpha1 = -1*10;
    //--------------------------
    // Time-Integration Setup  |
    //-------------------------
    dt             = 1e-4;          // Time step
    eta            = 27*pi;         // Damping eta = 6p*pi*R_rm*(Rtm+Ra)^2*mu  [ Pn.nm.ns]
    KbT            = 4.1 ;          // Thermal fluctuation energy  [pn.nm]
    pCa   = 3.0;
    pCa50 = 5.5;
    nH    = 1.0;
    Ca50  = pow(10,-(pCa50-6));
    Ca    = pow(10,-(pCa-6));                                     
    phi_2 = phi_2_max*(phi_2_min-phi_2_max)*(1./(1+pow(Ca50/Ca,nH)));

    //-----------------------------------------
    // Initial conditions
    //-----------------------------------------
    phi_RU02[0] = 0;
    phi_RU03[0] = 0;
    phi_RU04[0] = 0;
    phi_RU05[0] = 0;
    phi_RU06[0] = 0;
    phi_RU07[0] = 0;
    phi_RU08[0] = 0;
    phi_RU09[0] = 0;
    phi_RU10[0] = 0;
    phi_RU11[0] = 0;
    phi_RU12[0] = 0;
    phi_RU13[0] = 0;
    phi_RU14[0] = 0;
    phi_RU15[0] = 0;
    phi_RU16[0] = 0;
    phi_RU17[0] = 0;
    phi_RU18[0] = 0;
    phi_RU19[0] = 0;
    phi_RU20[0] = 0;
    phi_RU21[0] = 0;
    phi_RU22[0] = 0;
    phi_RU23[0] = 0;
    phi_RU24[0] = 0;
    phi_RU25[0] = 0;
   //------------------------------------------

    //-----------------------
    // Time-Marching part 1
    //-----------------------
    for (int i = 1; i < N_steps; i++)   
        {

        

        phi_RU02[i] = phi_RU02[i]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU02[i-1],3) +
                                         (Alpha3)* pow(phi_RU02[i-1],2) +
                                         (Alpha2)* pow(phi_RU02[i-1],1) +
                                         (Alpha1)
                                         )
                                         
       
       
        phi_RU03[i] = phi_RU03[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU03[i-1],3) +
                                         (Alpha3)* pow(phi_RU03[i-1],2) +
                                         (Alpha2)* pow(phi_RU03[i-1],1) +
                                         (Alpha1)
                                         )
                                        
       
        phi_RU04[i] = phi_RU04[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU04[i-1],3) +
                                         (Alpha3)* pow(phi_RU04[i-1],2) +
                                         (Alpha2)* pow(phi_RU04[i-1],1) +
                                         (Alpha1)
                                         )
                                    
       
        phi_RU05[i] = phi_RU05[i]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU05[i-1],3) +
                                         (Alpha1)* pow(phi_RU05[i-1],2) +
                                         (Alpha1)* pow(phi_RU05[i-1],1) +
                                         (Alpha1)
                                         )
                                        
      

        phi_RU06[i] = phi_RU06[i]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU06[i-1],3) +
                                         (Alpha3)* pow(phi_RU06[i-1],2) +
                                         (Alpha2)* pow(phi_RU06[i-1],1) +
                                         (Alpha1)
                                         )
                                        
      
        phi_RU07[i] = phi_RU07[i]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU07[i-1],3) +
                                         (Alpha3)* pow(phi_RU07[i-1],2) +
                                         (Alpha2)* pow(phi_RU07[i-1],1) +
                                         (Alpha1)
                                         )
                                         
       

        phi_RU08[i] = phi_RU08[i]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU08[i-1],3) +
                                         (Alpha1)* pow(phi_RU08[i-1],2) +
                                         (Alpha1)* pow(phi_RU08[i-1],1) +
                                         (Alpha1)
                                         )
                                     
      
        phi_RU09[i] = phi_RU09[i]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU09[i-1],3) +
                                         (Alpha3)* pow(phi_RU09[i-1],2) +
                                         (Alpha2)* pow(phi_RU09[i-1],1) +
                                         (Alpha1)
                                         )
                                         
  
        phi_RU10[i] = phi_RU10[i]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU10[i-1],3) +
                                         (Alpha3)* pow(phi_RU10[i-1],2) +
                                         (Alpha2)* pow(phi_RU10[i-1],2) +
                                         (Alpha1)
                                         )
                   

        phi_RU11[i] = phi_RU11[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU11[i-1],3) +
                                         (Alpha3)* pow(phi_RU11[i-1],2) +
                                         (Alpha2)* pow(phi_RU11[i-1],1) +
                                         (Alpha1)
                                         )
                                      
       
      

        phi_RU12[i] = phi_RU12[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU12[i-1],3) +
                                         (Alpha3)* pow(phi_RU12[i-1],2) +
                                         (Alpha2)* pow(phi_RU12[i-1],1) +
                                         (Alpha1)
                                         )
                                      
        phi_RU13[i] = phi_RU13[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU13[i-1],3) +
                                         (Alpha3)* pow(phi_RU13[i-1],2) +
                                         (Alpha2)* pow(phi_RU13[i-1],1) +
                                         (Alpha1)
                                         )
                                     
      

        phi_RU14[i] = phi_RU14[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU14[i-1],3) +
                                         (Alpha3)* pow(phi_RU14[i-1],2) +
                                         (Alpha2)* pow(phi_RU14[i-1],1) +
                                         (Alpha1)
                                         )
                         

        phi_RU15[i] = phi_RU15[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU15[i-1],3) +
                                         (Alpha3)* pow(phi_RU15[i-1],2) +
                                         (Alpha2)* pow(phi_RU15[i-1],1) +
                                         (Alpha1)
                                         )
                          
       

        phi_RU16[i] = phi_RU16[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU16[i-1],3) +
                                         (Alpha3)* pow(phi_RU16[i-1],2) +
                                         (Alpha2)* pow(phi_RU16[i-1],1) +
                                         (Alpha1)
                                         )
                     

        phi_RU17[i] = phi_RU17[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU17[i-1],3) +
                                         (Alpha3)* pow(phi_RU17[i-1],2) +
                                         (Alpha2)* pow(phi_RU17[i-1],1) +
                                         (Alpha1)
                                         )
           
       

        phi_RU18[i] = phi_RU18[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU18[i-1],3) +
                                         (Alpha3)* pow(phi_RU18[i-1],2) +
                                         (Alpha2)* pow(phi_RU18[i-1],1) +
                                         (Alpha1)
                                         )
                    

        phi_RU19[i] = phi_RU19[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU19[i-1],3) +
                                         (Alpha3)* pow(phi_RU19[i-1],2) +
                                         (Alpha2)* pow(phi_RU19[i-1],1) +
                                         (Alpha1)
                                         )
                       
      
        phi_RU20[i] = phi_RU20[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU20[i-1],3) +
                                         (Alpha3)* pow(phi_RU20[i-1],2) +
                                         (Alpha2)* pow(phi_RU20[i-1],1) +
                                         (Alpha1)
                                         )
             
       
        phi_RU21[i] = phi_RU21[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU21[i-1],3) +
                                         (Alpha3)* pow(phi_RU21[i-1],2) +
                                         (Alpha2)* pow(phi_RU21[i-1],1) +
                                         (Alpha1)
                                         )
                     
       
        phi_RU22[i] = phi_RU22[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU22[i-1],3) +
                                         (Alpha3)* pow(phi_RU22[i-1],2) +
                                         (Alpha2)* pow(phi_RU22[i-1],1) +
                                         (Alpha1)
                                         )
          
       

        phi_RU23[i] = phi_RU23[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU23[i-1],3) +
                                         (Alpha3)* pow(phi_RU23[i-1],2) +
                                         (Alpha2)* pow(phi_RU23[i-1],1) +
                                         (Alpha1)
                                         )
                     
       

        phi_RU24[i] = phi_RU24[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU24[i-1],3) +
                                         (Alpha3)* pow(phi_RU24[i-1],2) +
                                         (Alpha2)* pow(phi_RU24[i-1],1) +
                                         (Alpha1)
                                         )
         
      


       

    //-----------------------
    // Time-Marching part 21
    //-----------------------
    for (int i = 1; i < N_steps; i++)   
        {

       
        phi_RU02[i] = phi_RU02[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU02[i-1],3) +
                                         (Alpha3)* pow(phi_RU02[i-1],2) +
                                         (Alpha2)* pow(phi_RU02[i-1],1) +
                                         (Alpha1)
                                         )
                                        

        phi_RU01[i] = phi_RU03[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU03[i-1],3) +
                                         (Alpha3)* pow(phi_RU03[i-1],2) +
                                         (Alpha2)* pow(phi_RU03[i-1],1) +
                                         (Alpha1)
                                         )
                                        
       
        phi_RU04[i] = phi_RU04[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU04[i-1],3) +
                                         (Alpha3)* pow(phi_RU04[i-1],2) +
                                         (Alpha2)* pow(phi_RU04[i-1],1) +
                                         (Alpha1)
                                         )
                                       
       
        phi_RU05[i] = phi_RU05[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU05[i-1],3) +
                                         (Alpha3)* pow(phi_RU05[i-1],2) +
                                         (Alpha2)* pow(phi_RU05[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU05[i-1]-phi_RU04[i-1]-phi_RU06[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       

        phi_RU06[i] = phi_RU06[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU06[i-1],3) +
                                         (Alpha3)* pow(phi_RU06[i-1],2) +
                                         (Alpha2)* pow(phi_RU06[i-1],1) +
                                         (Alpha1)
                                         )
                                        
        phi_RU07[i] = phi_RU07[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU07[i-1],3) +
                                         (Alpha3)* pow(phi_RU07[i-1],2) +
                                         (Alpha2)* pow(phi_RU07[i-1],1) +
                                         (Alpha1)
                                         )
                                    

        phi_RU08[i] = phi_RU08[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU08[i-1],3) +
                                         (Alpha3)* pow(phi_RU08[i-1],2) +
                                         (Alpha2)* pow(phi_RU08[i-1],1) +
                                         (Alpha1)
                                         )
                                        
       
        phi_RU09[i] = phi_RU09[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU09[i-1],3) +
                                         (Alpha3)* pow(phi_RU09[i-1],2) +
                                         (Alpha2)* pow(phi_RU09[i-1],1) +
                                         (Alpha1)
                                         )
                                       
       

        phi_RU10[i] = phi_RU10[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU10[i-1],3) +
                                         (Alpha3)* pow(phi_RU10[i-1],2) +
                                         (Alpha2)* pow(phi_RU10[i-1],1) +
                                         (Alpha1)
                                         )
                                      
       
        phi_RU11[i] = phi_RU11[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU11[i-1],3) +
                                         (Alpha3)* pow(phi_RU11[i-1],2) +
                                         (Alpha2)* pow(phi_RU11[i-1],1) +
                                         (Alpha1)
                                         )
                                      
       

        phi_RU12[i] = phi_RU12[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU12[i-1],3) +
                                         (Alpha3)* pow(phi_RU12[i-1],2) +
                                         (Alpha2)* pow(phi_RU12[i-1],1) +
                                         (Alpha1)
                                         )
                            
       
        phi_RU13[i] = phi_RU13[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU13[i-1],3) +
                                         (Alpha3)* pow(phi_RU13[i-1],2) +
                                         (Alpha2)* pow(phi_RU13[i-1],1) +
                                         (Alpha1)
                                         )
                                      
      
        phi_RU14[i] = phi_RU14[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU14[i-1],3) +
                                         (Alpha3)* pow(phi_RU14[i-1],2) +
                                         (Alpha2)* pow(phi_RU14[i-1],1) +
                                         (Alpha1)
                                         )
                                    
       
        phi_RU15[i] = phi_RU15[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU15[i-1],3) +
                                         (Alpha3)* pow(phi_RU15[i-1],2) +
                                         (Alpha2)* pow(phi_RU15[i-1],1) +
                                         (Alpha1)
                                         )
                                        
        phi_RU16[i] = phi_RU16[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU16[i-1],3) +
                                         (Alpha3)* pow(phi_RU16[i-1],2) +
                                         (Alpha2)* pow(phi_RU16[i-1],1) +
                                         (Alpha1)
                                         )
                                         
        phi_RU17[i] = phi_RU17[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU17[i-1],3) +
                                         (Alpha3)* pow(phi_RU17[i-1],2) +
                                         (Alpha2)* pow(phi_RU17[i-1],1) +
                                         (Alpha1)
                                         )
                                        

        phi_RU18[i] = phi_RU18[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU18[i-1],3) +
                                         (Alpha3)* pow(phi_RU18[i-1],2) +
                                         (Alpha2)* pow(phi_RU18[i-1],1) +
                                         (Alpha1)
                                         )
                                        

        phi_RU19[i] = phi_RU19[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU19[i-1],3) +
                                         (Alpha3)* pow(phi_RU19[i-1],2) +
                                         (Alpha2)* pow(phi_RU19[i-1],1) +
                                         (Alpha1)
                                         )
                                      
       

        phi_RU20[i] = phi_RU20[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU20[i-1],3) +
                                         (Alpha3)* pow(phi_RU20[i-1],2) +
                                         (Alpha2)* pow(phi_RU20[i-1],1) +
                                         (Alpha1)
                                         )
                                        

        phi_RU21[i] = phi_RU21[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU21[i-1],3) +
                                         (Alpha3)* pow(phi_RU21[i-1],2) +
                                         (Alpha2)* pow(phi_RU21[i-1],1) +
                                         (Alpha1)
                                         )
                             

        phi_RU22[i] = phi_RU22[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU22[i-1],3) +
                                         (Alpha3)* pow(phi_RU22[i-1],2) +
                                         (Alpha2)* pow(phi_RU22[i-1],1) +
                                         (Alpha1)
                                         )
                                      
       
        phi_RU23[i] = phi_RU23[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU23[i-1],3) +
                                         (Alpha3)* pow(phi_RU23[i-1],2) +
                                         (Alpha2)* pow(phi_RU23[i-1],1) +
                                         (Alpha1)
                                         )
                            

        phi_RU24[i] = phi_RU24[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU24[i-1],3) +
                                         (Alpha3)* pow(phi_RU24[i-1],2) +
                                         (Alpha2)* pow(phi_RU24[i-1],1) +
                                         (Alpha1)
                                         )
                                       
       
        phi_RU25[i] = phi_RU25[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU25[i-1],3) +
                                         (Alpha3)* pow(phi_RU25[i-1],2) +
                                         (Alpha2)* pow(phi_RU25[i-1],1) +
                                         (Alpha1)
                                         )
      



    
    return 0;
}







