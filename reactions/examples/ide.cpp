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

using namespace std;

using std::ifstream;
using std::string;
using std::istringstream;

vector<vector<double> > mypk;

/* Example code to output the halo model powerspectrum for modified gravity */

int main(int argc, char* argv[]) {
  // Which gravity or dark energy model?
  // 1: GR  2: f(R) 3: DGP 4: quintessence with interaction 5: CPL with interaction
  int mymodel = 4;

  // Modified gravity active?
  bool modg = false;
  // Is the transfer being fed to ReACT of the target cosmology? If false, the transfer should be LCDM at z=0.
  bool mgcamb = false;

	 //output file name
    const char* output = "pl_wcdm11_xi10_z0.dat";
    const char* cstr = "transfers/DS/wcdm11";
// 0: scale factor, 1: omega_total, 2-4: mg param (1e-10 ~ GR for default mg functions ), 5: number of points in halo-mass loop in scol_init , 30 works well.
double vars[7];

  // chosen redshift
    double myz = 0.;
  // chosen omega_matter (total)
    double omega0 = 0.308;

    vars[0] = 1./(1.+myz);
    vars[1] =  omega0;

    vars[2] = -1.1;
    vars[3] = 0.; // wa for CPL
    vars[4] = 0.678*10.; // friction term in IDE

    vars[5] = 50.; // number of mass bins
    vars[6] = 0.0; // omega_neutrinos

    /* Open output file */
    FILE* fp = fopen(output, "w");

    // Keep it z=0 to keep Copter's Growth @ 1
    real z = 0;
    // Relative error in magnitude integrations
    real epsrel = 1e-3;

    Cosmology C(cstr);
    LinearPS P_l(C, z);
    HALO halo(C, P_l,P_l,P_l,P_l, epsrel);
    SPT spt(C, P_l, epsrel);

    IOW iow;
    real p1,p2,p3,p4,p5,p6;


halo.initialise(vars, mgcamb, modg, mymodel);
halo.phinit_pseudo(vars,mgcamb);


//int Nk = mypk.size();
int Nk = 500;
double kmin = 0.001;
double kmax = 100.;

 for(int i =0; i < Nk;  i ++) {

  real k =  kmin * exp(i*log(kmax/kmin)/(Nk-1));

        p1 = halo.plinear_cosmosis(k); // Linear spectrum
        p2 = halo.reaction_nu(k,vars); // halo model reaction
        p3 = halo.PHALO_pseudo(k,mgcamb); // halofit pseudo spectrum

     printf("%d %e %e %e %e \n", i, k, p1, p2, p3); // print to terminal
     fprintf(fp,"%e %e %e %e \n", k, p1,p2, p3); // print to file

}

	/*close output file*/
    fclose(fp);
    return 0;
}
