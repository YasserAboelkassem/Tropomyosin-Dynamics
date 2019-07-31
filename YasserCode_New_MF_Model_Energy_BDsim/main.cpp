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
#include <chrono>
#include <random>

using namespace std;
const long int N_steps  = 3e6;           //  Steps for both 1st and 2dn parts of the sol.
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
    normal_distribution<double> distribution (0.0,1.0);

    //--------------------------------
    A       = 1;                        // potential amplitude [Pn.nm]
    K_stiff = 1;                        // Stiff [pn.nm/rad]
    phi_B = 00.0 * pi/180;              // Blocked state "B" angle [rad]
    phi_C = 25.0 * pi/180;              // Closed  state "C" angle [rad]
    phi_M = 35.0 * pi/180;              // Open   state "M" angle  [rad]
    phi_1 = 27.0 * pi/180;              // This will be varied with Ca++ concentration [rad]
    //phi_2 = 2.0  * pi/180;            // This will be varied with Ca++ concentration [rad]

    //phi_1_min = 26.0 * pi/180;
    //phi_1_max = 34.0 * pi/180;
    // phi_2_min = 1.00 * pi/180;
    // phi_2_max = 24.0 * pi/180;

    //-------------------------------
    //T1 = phi_1 + phi_2 +phi_B + phi_M;
    //T2 = phi_1*phi_2 + (phi_1+phi_2)*(phi_B+phi_M) + phi_B*phi_M;
    //T3 = (phi_1*phi_2)*(phi_B+phi_M) + phi_B*phi_M * (phi_1+phi_2);
    //T4 = phi_1*phi_2*phi_B*phi_M;
    //----------------------------
    Alpha4 =  4*1100;
    Alpha3 = -3*1100;
    Alpha2 =  2*310;
    Alpha1 = -1*20;
    //--------------------------
    // Time-Integration Setup  |
    //-------------------------
    dt             = 1e-6;          // Time step
    eta            = 27*pi;         // Damping eta = 6p*pi*R_rm*(Rtm+Ra)^2*mu  [ Pn.nm.ns]
    KbT            = 4.1 ;          // Thermal fluctuation energy  [pn.nm]
    pCa   = 3.0;
    pCa50 = 5.5;
    nH    = 1.0;
    Ca50  = pow(10,-(pCa50-6));
    Ca    = pow(10,-(pCa-6));                                     
    phi_2 = phi_2_max+(phi_2_min-phi_2_max)*(1./(1+pow(Ca50/Ca,nH)));

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

        //-------------------------------------------------------------------------------------------------
        // RU2 |
        //-------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU02[i] = phi_RU02[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU02[i-1],3) +
                                         (Alpha3)* pow(phi_RU02[i-1],2) +
                                         (Alpha2)* pow(phi_RU02[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU02[i-1]-phi_RU03[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU3  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU03[i] = phi_RU03[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU03[i-1],3) +
                                         (Alpha3)* pow(phi_RU03[i-1],2) +
                                         (Alpha2)* pow(phi_RU03[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU03[i-1]-phi_RU02[i-1]-phi_RU04[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU4  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU04[i] = phi_RU04[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU04[i-1],3) +
                                         (Alpha3)* pow(phi_RU04[i-1],2) +
                                         (Alpha2)* pow(phi_RU04[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU04[i-1]-phi_RU03[i-1]-phi_RU05[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU5  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU05[i] = phi_RU05[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU05[i-1],3) +
                                         (Alpha3)* pow(phi_RU05[i-1],2) +
                                         (Alpha2)* pow(phi_RU05[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU05[i-1]-phi_RU04[i-1]-phi_RU06[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU6  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU06[i] = phi_RU06[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU06[i-1],3) +
                                         (Alpha3)* pow(phi_RU06[i-1],2) +
                                         (Alpha2)* pow(phi_RU06[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU06[i-1]-phi_RU05[i-1]-phi_RU07[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU7  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU07[i] = phi_RU07[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU07[i-1],3) +
                                         (Alpha3)* pow(phi_RU07[i-1],2) +
                                         (Alpha2)* pow(phi_RU07[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU07[i-1]-phi_RU06[i-1]-phi_RU08[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU8  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU08[i] = phi_RU08[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU08[i-1],3) +
                                         (Alpha3)* pow(phi_RU08[i-1],2) +
                                         (Alpha2)* pow(phi_RU08[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU08[i-1]-phi_RU07[i-1]-phi_RU09[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU9  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU09[i] = phi_RU09[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU09[i-1],3) +
                                         (Alpha3)* pow(phi_RU09[i-1],2) +
                                         (Alpha2)* pow(phi_RU09[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU09[i-1]-phi_RU08[i-1]-phi_RU10[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU10  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU10[i] = phi_RU10[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU10[i-1],3) +
                                         (Alpha3)* pow(phi_RU10[i-1],2) +
                                         (Alpha2)* pow(phi_RU10[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU10[i-1]-phi_RU09[i-1]-phi_RU11[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------


       //--------------------------------------------------------------------------------------------------
       // RU11  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU11[i] = phi_RU11[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU11[i-1],3) +
                                         (Alpha3)* pow(phi_RU11[i-1],2) +
                                         (Alpha2)* pow(phi_RU11[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU11[i-1]-phi_RU10[i-1]-phi_RU12[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU12  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU12[i] = phi_RU12[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU12[i-1],3) +
                                         (Alpha3)* pow(phi_RU12[i-1],2) +
                                         (Alpha2)* pow(phi_RU12[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU12[i-1]-phi_RU11[i-1]-phi_RU13[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU13  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU13[i] = phi_RU13[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU13[i-1],3) +
                                         (Alpha3)* pow(phi_RU13[i-1],2) +
                                         (Alpha2)* pow(phi_RU13[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU13[i-1]-phi_RU12[i-1]-phi_RU14[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU14  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU14[i] = phi_RU14[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU14[i-1],3) +
                                         (Alpha3)* pow(phi_RU14[i-1],2) +
                                         (Alpha2)* pow(phi_RU14[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU14[i-1]-phi_RU13[i-1]-phi_RU15[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU15  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU15[i] = phi_RU15[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU15[i-1],3) +
                                         (Alpha3)* pow(phi_RU15[i-1],2) +
                                         (Alpha2)* pow(phi_RU15[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU15[i-1]-phi_RU14[i-1]-phi_RU16[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU16  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU16[i] = phi_RU16[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU16[i-1],3) +
                                         (Alpha3)* pow(phi_RU16[i-1],2) +
                                         (Alpha2)* pow(phi_RU16[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU16[i-1]-phi_RU15[i-1]-phi_RU17[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU17  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU17[i] = phi_RU17[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU17[i-1],3) +
                                         (Alpha3)* pow(phi_RU17[i-1],2) +
                                         (Alpha2)* pow(phi_RU17[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU17[i-1]-phi_RU16[i-1]-phi_RU18[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU18  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU18[i] = phi_RU18[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU18[i-1],3) +
                                         (Alpha3)* pow(phi_RU18[i-1],2) +
                                         (Alpha2)* pow(phi_RU18[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU18[i-1]-phi_RU17[i-1]-phi_RU19[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU19  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU19[i] = phi_RU19[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU19[i-1],3) +
                                         (Alpha3)* pow(phi_RU19[i-1],2) +
                                         (Alpha2)* pow(phi_RU19[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU19[i-1]-phi_RU18[i-1]-phi_RU20[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU20  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU20[i] = phi_RU20[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU20[i-1],3) +
                                         (Alpha3)* pow(phi_RU20[i-1],2) +
                                         (Alpha2)* pow(phi_RU20[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU20[i-1]-phi_RU19[i-1]-phi_RU21[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU21  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU21[i] = phi_RU21[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU21[i-1],3) +
                                         (Alpha3)* pow(phi_RU21[i-1],2) +
                                         (Alpha2)* pow(phi_RU21[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU21[i-1]-phi_RU20[i-1]-phi_RU22[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU22  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU22[i] = phi_RU22[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU22[i-1],3) +
                                         (Alpha3)* pow(phi_RU22[i-1],2) +
                                         (Alpha2)* pow(phi_RU22[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU22[i-1]-phi_RU21[i-1]-phi_RU23[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU23  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU23[i] = phi_RU23[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU23[i-1],3) +
                                         (Alpha3)* pow(phi_RU23[i-1],2) +
                                         (Alpha2)* pow(phi_RU23[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU23[i-1]-phi_RU22[i-1]-phi_RU24[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU24  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU24[i] = phi_RU24[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU24[i-1],3) +
                                         (Alpha3)* pow(phi_RU24[i-1],2) +
                                         (Alpha2)* pow(phi_RU24[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU24[i-1]-phi_RU23[i-1]-phi_RU25[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU25  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU25[i] = phi_RU25[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU25[i-1],3) +
                                         (Alpha3)* pow(phi_RU25[i-1],2) +
                                         (Alpha2)* pow(phi_RU25[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU25[i-1]-phi_RU24[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------


       

    //-----------------------
    // Time-Marching part 21
    //-----------------------
    for (int i = 1; i < N_steps; i++)   
        {

        //-------------------------------------------------------------------------------------------------
        // RU2 |
        //-------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU02[i] = phi_RU02[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU02[i-1],3) +
                                         (Alpha3)* pow(phi_RU02[i-1],2) +
                                         (Alpha2)* pow(phi_RU02[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU02[i-1]-phi_RU03[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU3  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU03[i] = phi_RU03[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU03[i-1],3) +
                                         (Alpha3)* pow(phi_RU03[i-1],2) +
                                         (Alpha2)* pow(phi_RU03[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU03[i-1]-phi_RU02[i-1]-phi_RU04[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU4  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU04[i] = phi_RU04[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU04[i-1],3) +
                                         (Alpha3)* pow(phi_RU04[i-1],2) +
                                         (Alpha2)* pow(phi_RU04[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU04[i-1]-phi_RU03[i-1]-phi_RU05[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU5  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU05[i] = phi_RU05[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU05[i-1],3) +
                                         (Alpha3)* pow(phi_RU05[i-1],2) +
                                         (Alpha2)* pow(phi_RU05[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU05[i-1]-phi_RU04[i-1]-phi_RU06[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU6  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU06[i] = phi_RU06[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU06[i-1],3) +
                                         (Alpha3)* pow(phi_RU06[i-1],2) +
                                         (Alpha2)* pow(phi_RU06[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU06[i-1]-phi_RU05[i-1]-phi_RU07[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU7  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU07[i] = phi_RU07[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU07[i-1],3) +
                                         (Alpha3)* pow(phi_RU07[i-1],2) +
                                         (Alpha2)* pow(phi_RU07[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU07[i-1]-phi_RU06[i-1]-phi_RU08[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU8  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU08[i] = phi_RU08[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU08[i-1],3) +
                                         (Alpha3)* pow(phi_RU08[i-1],2) +
                                         (Alpha2)* pow(phi_RU08[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU08[i-1]-phi_RU07[i-1]-phi_RU09[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU9  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU09[i] = phi_RU09[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU09[i-1],3) +
                                         (Alpha3)* pow(phi_RU09[i-1],2) +
                                         (Alpha2)* pow(phi_RU09[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU09[i-1]-phi_RU08[i-1]-phi_RU10[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU10  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU10[i] = phi_RU10[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU10[i-1],3) +
                                         (Alpha3)* pow(phi_RU10[i-1],2) +
                                         (Alpha2)* pow(phi_RU10[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU10[i-1]-phi_RU09[i-1]-phi_RU11[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------


       //--------------------------------------------------------------------------------------------------
       // RU11  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU11[i] = phi_RU11[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU11[i-1],3) +
                                         (Alpha3)* pow(phi_RU11[i-1],2) +
                                         (Alpha2)* pow(phi_RU11[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU11[i-1]-phi_RU10[i-1]-phi_RU12[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU12  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU12[i] = phi_RU12[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU12[i-1],3) +
                                         (Alpha3)* pow(phi_RU12[i-1],2) +
                                         (Alpha2)* pow(phi_RU12[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU12[i-1]-phi_RU11[i-1]-phi_RU13[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU13  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU13[i] = phi_RU13[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU13[i-1],3) +
                                         (Alpha3)* pow(phi_RU13[i-1],2) +
                                         (Alpha2)* pow(phi_RU13[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU13[i-1]-phi_RU12[i-1]-phi_RU14[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU14  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU14[i] = phi_RU14[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU14[i-1],3) +
                                         (Alpha3)* pow(phi_RU14[i-1],2) +
                                         (Alpha2)* pow(phi_RU14[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU14[i-1]-phi_RU13[i-1]-phi_RU15[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU15  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU15[i] = phi_RU15[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU15[i-1],3) +
                                         (Alpha3)* pow(phi_RU15[i-1],2) +
                                         (Alpha2)* pow(phi_RU15[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU15[i-1]-phi_RU14[i-1]-phi_RU16[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU16  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU16[i] = phi_RU16[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU16[i-1],3) +
                                         (Alpha3)* pow(phi_RU16[i-1],2) +
                                         (Alpha2)* pow(phi_RU16[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU16[i-1]-phi_RU15[i-1]-phi_RU17[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU17  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU17[i] = phi_RU17[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU17[i-1],3) +
                                         (Alpha3)* pow(phi_RU17[i-1],2) +
                                         (Alpha2)* pow(phi_RU17[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU17[i-1]-phi_RU16[i-1]-phi_RU18[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU18  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU18[i] = phi_RU18[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU18[i-1],3) +
                                         (Alpha3)* pow(phi_RU18[i-1],2) +
                                         (Alpha2)* pow(phi_RU18[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU18[i-1]-phi_RU17[i-1]-phi_RU19[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU19  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU19[i] = phi_RU19[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU19[i-1],3) +
                                         (Alpha3)* pow(phi_RU19[i-1],2) +
                                         (Alpha2)* pow(phi_RU19[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU19[i-1]-phi_RU18[i-1]-phi_RU20[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU20  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU20[i] = phi_RU20[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU20[i-1],3) +
                                         (Alpha3)* pow(phi_RU20[i-1],2) +
                                         (Alpha2)* pow(phi_RU20[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU20[i-1]-phi_RU19[i-1]-phi_RU21[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU21  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU21[i] = phi_RU21[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU21[i-1],3) +
                                         (Alpha3)* pow(phi_RU21[i-1],2) +
                                         (Alpha2)* pow(phi_RU21[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU21[i-1]-phi_RU20[i-1]-phi_RU22[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU22  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU22[i] = phi_RU22[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU22[i-1],3) +
                                         (Alpha3)* pow(phi_RU22[i-1],2) +
                                         (Alpha2)* pow(phi_RU22[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU22[i-1]-phi_RU21[i-1]-phi_RU23[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU23  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU23[i] = phi_RU23[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU23[i-1],3) +
                                         (Alpha3)* pow(phi_RU23[i-1],2) +
                                         (Alpha2)* pow(phi_RU23[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU23[i-1]-phi_RU22[i-1]-phi_RU24[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU24  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU24[i] = phi_RU24[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU24[i-1],3) +
                                         (Alpha3)* pow(phi_RU24[i-1],2) +
                                         (Alpha2)* pow(phi_RU24[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU24[i-1]-phi_RU23[i-1]-phi_RU25[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------

       //--------------------------------------------------------------------------------------------------
       // RU25  |
       //--------------------------------------------------------------------------------------------------
        W_noise = distribution(generator);

        phi_RU25[i] = phi_RU25[i-1]-
                            (dt*A/eta) *(
                                         (Alpha4)* pow(phi_RU25[i-1],3) +
                                         (Alpha3)* pow(phi_RU25[i-1],2) +
                                         (Alpha2)* pow(phi_RU25[i-1],1) +
                                         (Alpha1)
                                         )
                                         - (dt*K_stiff/eta)*(2*phi_RU25[i-1]-phi_RU24[i-1])
                                         + pow(2*dt,0.5)* pow(KbT/eta,0.5) * W_noise;
       //--------------------------------------------------------------------------------------------------


        } 

    cout << " Last value of phi_RU02 = " << phi_RU02[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU03 = " << phi_RU03[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU04 = " << phi_RU04[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU05 = " << phi_RU05[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU06 = " << phi_RU06[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU07 = " << phi_RU07[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU08 = " << phi_RU08[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU09 = " << phi_RU09[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU10 = " << phi_RU10[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU11 = " << phi_RU11[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU12 = " << phi_RU12[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU13 = " << phi_RU13[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU14 = " << phi_RU14[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU15 = " << phi_RU15[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU16 = " << phi_RU16[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU17 = " << phi_RU17[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU18 = " << phi_RU18[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU19 = " << phi_RU19[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU20 = " << phi_RU20[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU21 = " << phi_RU21[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU22 = " << phi_RU22[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU23 = " << phi_RU23[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU24 = " << phi_RU24[N_steps-1]*180/pi << endl;
    cout << " Last value of phi_RU25 = " << phi_RU25[N_steps-1]*180/pi << endl;

    //------------------------------------------------
    // Write data into a text file
    //------------------------------------------------
    ofstream RU02_Phi_Time_out2 ("RU02_Phi_Time.csv"); 
    ofstream RU03_Phi_Time_out3 ("RU03_Phi_Time.csv");
    ofstream RU04_Phi_Time_out4 ("RU04_Phi_Time.csv");
    ofstream RU05_Phi_Time_out5 ("RU05_Phi_Time.csv");
    ofstream RU06_Phi_Time_out6 ("RU06_Phi_Time.csv");
    ofstream RU07_Phi_Time_out7 ("RU07_Phi_Time.csv");
    ofstream RU08_Phi_Time_out8 ("RU08_Phi_Time.csv");
    ofstream RU09_Phi_Time_out9 ("RU09_Phi_Time.csv");
    ofstream RU10_Phi_Time_out10("RU10_Phi_Time.csv");
    ofstream RU11_Phi_Time_out11("RU11_Phi_Time.csv");
    ofstream RU12_Phi_Time_out12("RU12_Phi_Time.csv");
    ofstream RU13_Phi_Time_out13("RU13_Phi_Time.csv");
    ofstream RU14_Phi_Time_out14("RU14_Phi_Time.csv");
    ofstream RU15_Phi_Time_out15("RU15_Phi_Time.csv");
    ofstream RU16_Phi_Time_out16("RU16_Phi_Time.csv");
    ofstream RU17_Phi_Time_out17("RU17_Phi_Time.csv");
    ofstream RU18_Phi_Time_out18("RU18_Phi_Time.csv");
    ofstream RU19_Phi_Time_out19("RU19_Phi_Time.csv");
    ofstream RU20_Phi_Time_out20("RU20_Phi_Time.csv");
    ofstream RU21_Phi_Time_out21("RU21_Phi_Time.csv");
    ofstream RU22_Phi_Time_out22("RU22_Phi_Time.csv");
    ofstream RU23_Phi_Time_out23("RU23_Phi_Time.csv");
    ofstream RU24_Phi_Time_out24("RU24_Phi_Time.csv");
    ofstream RU25_Phi_Time_out25("RU25_Phi_Time.csv");

    int increment = 1;   
    for (int i = 0; i < N_steps; i+=increment)  
		    {
                RU02_Phi_Time_out2 << phi_RU02[i]*180/pi << endl; 
                RU03_Phi_Time_out3 << phi_RU03[i]*180/pi << endl;
                RU04_Phi_Time_out4 << phi_RU04[i]*180/pi << endl;
                RU05_Phi_Time_out5 << phi_RU05[i]*180/pi << endl;
                RU06_Phi_Time_out6 << phi_RU06[i]*180/pi << endl;
                RU07_Phi_Time_out7 << phi_RU07[i]*180/pi << endl;

                RU08_Phi_Time_out8  << phi_RU08[i]*180/pi << endl;
                RU09_Phi_Time_out9  << phi_RU09[i]*180/pi << endl;
                RU10_Phi_Time_out10 << phi_RU10[i]*180/pi << endl;
                RU11_Phi_Time_out11 << phi_RU11[i]*180/pi << endl;
                RU12_Phi_Time_out12 << phi_RU12[i]*180/pi << endl;
                RU13_Phi_Time_out13 << phi_RU13[i]*180/pi << endl;

                RU14_Phi_Time_out14 << phi_RU14[i]*180/pi << endl;
                RU15_Phi_Time_out15 << phi_RU15[i]*180/pi << endl;
                RU16_Phi_Time_out16 << phi_RU16[i]*180/pi << endl;
                RU17_Phi_Time_out17 << phi_RU17[i]*180/pi << endl;
                RU18_Phi_Time_out18 << phi_RU18[i]*180/pi << endl;
                RU19_Phi_Time_out19 << phi_RU19[i]*180/pi << endl;

                RU20_Phi_Time_out20 << phi_RU20[i]*180/pi << endl;
                RU21_Phi_Time_out21 << phi_RU21[i]*180/pi << endl;
                RU22_Phi_Time_out22 << phi_RU22[i]*180/pi << endl;
                RU23_Phi_Time_out23 << phi_RU23[i]*180/pi << endl;
                RU24_Phi_Time_out24 << phi_RU24[i]*180/pi << endl;
                RU25_Phi_Time_out25 << phi_RU25[i]*180/pi << endl;

            }

    cout << " data successfully saved " << endl;

    return 0;
}







