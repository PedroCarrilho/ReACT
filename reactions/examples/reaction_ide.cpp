/* Compile with

 g++ -I/usr/local/include -I/Users/pedro/Documents/GitHub/ReACT/reactions/include -L/usr/local/lib -L/Users/pedro/Documents/GitHub/ReACT/reactions/lib reaction_ide.cpp -lgsl -lcopter -lstdc++ -o test

*/


#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <cmath>


#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <cmath>

#include <omp.h>

#include <Copter/HALO.h>
#include <Copter/Cosmology.h>
#include <Copter/LinearPS.h>
#include <Copter/SpecialFunctions.h>
#include <Copter/SPT.h>

//using namespace std;
using std::ifstream;
using std::string;
using std::istringstream;

vector<vector<double> > mytrans;
vector<vector<double> > mytransl;
vector<vector<double> > mypk;


/* Example code to output the reaction and halo spectra for mg + neutrinos */
int main(int argc, char* argv[]) {

// Which gravity or dark energy model?
// 1: GR  2: f(R) 3: DGP 4: quintessence 5: CPL
int mymodel = 4;

// target redshift
double myz = 0.;

double Omega_nu = 0.00;  // neutrino fraction mv = 0.0ev
//double Omega_nu = 0.0053;  // neutrino fraction mv = 0.24

// Modified gravity active? This allows k* and \mathcal{E} to take on non LCDM values.
bool modg = false;

// Is the transfer being fed to ReACT of the target cosmology?
//If false, the transfer should be LCDM at z=0 and ReACT will rescale P_L using internally computed modified growth - see README.
// If true, the transfer function should be that of the real cosmology (with MG or/and massive neutrinos)
// Note that ReACT does not calculate growth factors for massive neutrino cosmologies and so the real transfer function should be supplied.
bool mgcamb = false;

//output file name
const char* output = "wcdm_z0_w09_xi10_test_virial.dat";
//const char* output = "wcdm_z0_lcdm.dat";

const char* cstr = "transfers/wcdm09";
//const char* cstr = "transfers/lcdm_for_ide";

// integration error
real epsrel = 1e-3;

/*Specify params*/

/* Dustgrain */
double Omega_m = 0.308; // total matter fraction

// CPL parameters
double w0 = -0.9;
double wa = 0.0;

// number of mass bins between 5<Log10[M]<20
double massb = 50.;

// store params for passing into React functions
double vars[7];
    vars[0] = 1./(myz+1.); //  scale factor
    vars[1] = Omega_m;
    vars[2] = w0; //  modified gravity param or w0 (see SpecialFunctions.cpp)
    vars[3] = wa;  // extra, in the CPL case it is wa
    vars[4] = 10*0.678; // extra, is xi for IDE
    vars[5] = massb; // number of mass bins between 5<Log10[M]<18
    vars[6] = Omega_nu;


IOW iow;

// Load cosmology classes
Cosmology Cm(cstr);

// Get linear P(k) from input transfer (0. is the redshift for Copter, which I hear must be zero for consistency with other choices)
LinearPS P_l(Cm, 0.);

// Load halo class witth all linear P(k)
HALO halo(Cm, P_l,P_l,P_l,P_l, epsrel);
SPT spt(Cm, P_l, epsrel);


//initialise spherical collapse quantities and reaction quantities
halo.initialise(vars,mgcamb,modg,mymodel);
// initialise halofit parameters
halo.phinit_pseudo(vars,mgcamb);

double p1,p2,p3,p4,p5;
int Nk = 100;
double kmin = 1e-4;
double kmax = 10.;

 /* Open output file */
 FILE* fp = fopen(output, "w");

 for(int i =0; i <Nk;  i ++) {

      real k = kmin* exp(i*log(kmax/kmin)/(Nk-1));

      p1 = P_l(k); // Linear spectrum
      p2 = halo.reaction_nu(k,vars); // halo model reaction
      p3 = halo.PHALO_pseudo(k,mgcamb); // halofit pseudo spectrum

      printf("%e %e %e %e %e  \n", k, p1,p2, p3,p3*p2); // output to terminal
      fprintf(fp,"%e %e %e %e %e \n", k, p1,p2,p3,p3*p2); // output to file : k , P_linear, R, P_pseudo, P_nl = RxP_pseudo

}

	/*close output file*/
    fclose(fp);
    return 0;
}
