#include "nfit.h"
#include "modelcalculator.h"
#include "utable.h"
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <fstream>

using namespace std;

#define PI 3.1415926535897932384626433832795
#define SMALLNUM 0.00000000001


void testCCDStrFct()
{
  ModelCalculator mc;
	mc.read_in_utable("../dat/utab.dat");
	
  // Do not change the following parameters
	mc.setModelParameter(8e-13, "Kc"); //6e-13
	mc.setModelParameter(1e13, "B"); //2e13
	mc.setModelParameter(62.8, "D");
	mc.setModelParameter(30, "T"); //37
	mc.setModelParameter(5000, "Lr");
	mc.setModelParameter(10, "Mz");
	mc.setModelParameter(0., "edisp"); //0.
	mc.setModelParameter(0., "mosaic");
	mc.setModelParameter(1.175, "wavelength");
	mc.setModelParameter(0.07113, "pixelSize");
	mc.setModelParameter(2.3, "bFWHM"); //0.
	mc.setModelParameter(360, "s");
    mc.setModelParameter(95., "Kt"); //new
    mc.setModelParameter(12., "a"); //new
	

	
	//ofstream outfile("makeCCDStrFct_1_082714.dat");
	//Create a structure factor by looping over all desired qz-values
    //for (double qz = 0.05; qz <= 0.75; qz += 0.1) {
    //  cout << "qz= " << qz << endl;
	//	vector<double> qvec, sfvec;
	//	mc.getCCDStrFct(0, 0.2, qz, qvec, sfvec);
	//	for (vector<double>::size_type i = 0; i < qvec.size(); i++) {	
	//		outfile << qvec[i] << " " << qz << " " << sfvec[i] << endl;
    //	}
    //}
	//outfile.close();

	vector<double> qvec, sfvec;
	mc.getCCDStrFct(0, 0.2, 0.2, qvec, sfvec);

	for (vector<double>::size_type i = 0; i < qvec.size(); i++) {
		cout << qvec[i] << " " << sfvec[i] << endl;
	}
	cout << "===========================" << endl;
	cout << "qz: " << 0.2 << endl;
	cout << "# of points calculated: " << qvec.size() << endl;
	//saveDoubleColumns(qvec, sfvec, "test_CCD_struct_factor.dat");		
}


void testRotatedStrFct()
{
  ModelCalculator mc;
	mc.read_in_utable("../dat/utab.dat");
	
  // Do not change the following parameters
	mc.setModelParameter(6e-13, "Kc");
	mc.setModelParameter(2e13, "B");
	mc.setModelParameter(62.8, "D");
	mc.setModelParameter(37, "T");
	mc.setModelParameter(5000, "Lr");
	mc.setModelParameter(10, "Mz");
	mc.setModelParameter(0.0134, "edisp");
	mc.setModelParameter(0.1, "mosaic");
	mc.setModelParameter(1.175, "wavelength");
	mc.setModelParameter(0.07113, "pixelSize");
	mc.setModelParameter(2.3, "bFWHM");
	mc.setModelParameter(360, "s");
	
	vector<double> qvec, sfvec;
	mc.getRotatedStrFct(0, 0.2, 0.2, qvec, sfvec);
	
	for (vector<double>::size_type i = 0; i < qvec.size(); i++) {
		cout << qvec[i] << " " << sfvec[i] << endl;
	}
	cout << "===========================" << endl;
	cout << "qz: " << 0.2 << endl;
	cout << "# of points calculated: " << qvec.size() << endl;
	saveDoubleColumns(qvec, sfvec, "test_rotated_struct_factor.dat");		
}


void testMosaicStrFct()
{
	ModelCalculator mc;
	mc.read_in_utable("../dat/utab.dat");
	
	// Do not change the following parameters
	mc.setModelParameter(6e-13, "Kc");
	mc.setModelParameter(2e13, "B");
	mc.setModelParameter(62.8, "D");
	mc.setModelParameter(37, "T");
	mc.setModelParameter(5000, "Lr");
	mc.setModelParameter(10, "Mz");
	mc.setModelParameter(0.0134, "edisp");
	mc.setModelParameter(0.1, "mosaic");
	mc.setModelParameter(1.175, "wavelength");
	mc.setModelParameter(0.07113, "pixelSize");
	mc.setModelParameter(2.3, "bFWHM");
	mc.setModelParameter(360, "s");

	vector<double> qvec, sfvec;
	mc.getMosaicStrFct(0, 0.2, 0.2, qvec, sfvec);
	
	for (vector<double>::size_type i = 0; i < qvec.size(); i++) {
		cout << qvec[i] << " " << sfvec[i] << endl;
	}
	cout << "===========================" << endl;
	cout << "qz: " << 0.2 << endl;
	cout << "# of points calculated: " << qvec.size() << endl;
	saveDoubleColumns(qvec, sfvec, "test_mosaic_struct_factor.dat");
}


void testStrFct()
{
	ModelCalculator mc;
	mc.read_in_utable("../dat/utab.dat");
	  
	mc.setModelParameter(8e-13, "Kc");
	mc.setModelParameter(1e13, "B");
	mc.setModelParameter(62.8, "D");
	mc.setModelParameter(30, "T");
	mc.setModelParameter(5000, "Lr");
	mc.setModelParameter(10, "Mz");
	mc.setModelParameter(0., "edisp");
	mc.setModelParameter(0., "mosaic");
	mc.setModelParameter(1.175, "wavelength");
	mc.setModelParameter(0.07113, "pixelSize");
	mc.setModelParameter(0., "bFWHM");
	mc.setModelParameter(360, "s");
    mc.setModelParameter(80, "Kt"); //new
    mc.setModelParameter(16., "at"); //new

    for (double qz = 0.106; qz <= 0.116; qz += 0.001) {
	  vector<double> qvec, sfvec;
	  mc.getStrFct(0., 0.2, qz, qvec, sfvec);
      cout << "qz= " << qz << endl;
	  //for (vector<double>::size_type i = 0; i < qvec.size(); i++) {  
	  //}
    }

	/*vector<double> qvec, sfvec;
	mc.getStrFct(0, 0.2, 0.2, qvec, sfvec);
	
	for (vector<double>::size_type i = 0; i < qvec.size(); i++) {
		cout << qvec[i] << " " << sfvec[i] << endl;
	}
	cout << "===========================" << endl;
	cout << "qz: " << 0.2 << endl;
	cout << "# of points calculated: " << qvec.size() << endl;
	saveDoubleColumns(qvec, sfvec, "test_struct_factor.dat");*/
}


// Run to read from new utable and test it
void testUtable()
{
    ModelCalculator mc;
	mc.read_in_utable("../dat/utab.dat");
cout << "Done reading utab.dat" << endl;
//	mc.read_in_nstable("../dat/nstab.dat");
//cout << "Done reading nstab.dat" << endl;

    mc.setModelParameter(8e-13, "Kc");
	mc.setModelParameter(1e13, "B");
	mc.setModelParameter(62.8, "D");
	mc.setModelParameter(30, "T");
	mc.setModelParameter(5000, "Lr");
	mc.setModelParameter(10, "Mz");
	mc.setModelParameter(0., "edisp");
	mc.setModelParameter(0., "mosaic");
	mc.setModelParameter(1.175, "wavelength");
	mc.setModelParameter(0.07113, "pixelSize");
	mc.setModelParameter(0., "bFWHM");
	mc.setModelParameter(360, "s");
    mc.setModelParameter(80, "Kt"); //new
    mc.setModelParameter(8., "at"); //new

    /*int n = 0;
    double r = 1;
    double temp;
    cout << mc.get_un_value(n, r) << endl;*/
 
    /*double Kc = 6e-13;
    double B = 2e13;
    double Kt = Kc/100;
    double dspacing = 62.8;
    double a = 8.0;
    double T = 37.0;

    double xi = 1e8 * pow(Kc/B, 0.25);
    double xit_2 = Kc/Kt; //REMEMBER to include prefactor to go from cm^2 to A^2
    double ell = 2 * xit_2 / (xi * xi);
    double eta = 2.16872 * (T+273.15) / sqrt(Kc*B) / dspacing /dspacing;
    //double prefactor = eta * dspacing * dspacing / 2 / PI / PI;  */  

    int n = 0;
	double r = 0.064;
    double temp;
	temp = mc.get_un_value(n, r);
	cout << temp << endl;

    /*for(double r = 100; r < 101; r = r + 1) {
      temp = mc.get_un_value(n, r);
      //cout << "r: " << r << " n: " << n << endl;
      cout << temp << endl;
    }*/
}

void makeStrFct()
{
	ModelCalculator mc;
	mc.read_in_utable("../dat/utab.dat");
cout << "Done reading utab.dat" << endl;
//	mc.read_in_nstable("../dat/nstab.dat");
//cout << "Done reading nstab.dat" << endl;
	  
	mc.setModelParameter(7.8e-13, "Kc");
	mc.setModelParameter(8.5e12, "B"); // was 1e13
	mc.setModelParameter(62.995, "D");
	mc.setModelParameter(30, "T");
	mc.setModelParameter(5000, "Lr");
	mc.setModelParameter(10, "Mz");
	mc.setModelParameter(0., "edisp");
	mc.setModelParameter(0., "mosaic");
	mc.setModelParameter(1.177, "wavelength");
	mc.setModelParameter(0.07113, "pixelSize");
	mc.setModelParameter(0., "bFWHM");
	mc.setModelParameter(359.3, "s");
    mc.setModelParameter(96.2, "Kt"); //new
    mc.setModelParameter(14.6, "at"); //new

//cout << "sumn= " << mc.sumn(1) << endl;
//usleep(1000000);

/*	for (double qz = 0.2; qz <= 0.21; qz += 0.001) {
        cout << "qz= " << qz << endl;
		vector<double> qvec, sfvec;
		mc.getStrFct(0., 0.3, qz, qvec, sfvec);
    }*/

	ofstream outfile("makeStrFct_038_Kt96v3_061914.dat");
	//Create a structure factor by looping over all desired qz-values
    for (double qz = 0.05; qz <= 0.75; qz += 0.001) {
        cout << "qz= " << qz << endl;
		vector<double> qvec, sfvec;
		mc.getStrFct(0, 0.2, qz, qvec, sfvec);
		for (vector<double>::size_type i = 0; i < qvec.size(); i++) {	
			outfile << qvec[i] << " " << qz << " " << sfvec[i] << endl;
    	}
    }
	outfile.close();
}

void makeStrFct2()
{
	ModelCalculator mc;
	mc.read_in_utable("../dat/utab.dat");
cout << "Done reading utab.dat" << endl;
//	mc.read_in_nstable("../dat/nstab.dat");
//cout << "Done reading nstab.dat" << endl;
	  
	mc.setModelParameter(8.0e-13, "Kc");
	mc.setModelParameter(1e13, "B"); // was 1e13
	mc.setModelParameter(62.8, "D");
	mc.setModelParameter(30, "T");
	mc.setModelParameter(5000, "Lr");
	mc.setModelParameter(10, "Mz");
	mc.setModelParameter(0.0134, "edisp");
	mc.setModelParameter(0., "mosaic");
	mc.setModelParameter(1.175, "wavelength");
	mc.setModelParameter(0.07113, "pixelSize");
	mc.setModelParameter(2.3, "bFWHM");
	mc.setModelParameter(360.0, "s");
    mc.setModelParameter(95.0, "Kt"); //new
    mc.setModelParameter(12, "at"); //new
	mc.setModelParameter(5.0, "Ls"); //new 2/16/2015
	mc.setModelParameter(0.0001, "divergeX"); //new 2/26/2015
	mc.setModelParameter(0.0001, "divergeZ"); //new 2/26/2015
	mc.setModelParameter(260, "teff"); //new 3/30/2015

	// determine whether S_0(q) or S(q) is calculated
	mc.set_onlyZERO(1); // if 0 is passed S(q) is calculated
//cout << "sumn= " << mc.sumn(1) << endl;
//usleep(1000000);

/*	for (double qz = 0.2; qz <= 0.21; qz += 0.001) {
        cout << "qz= " << qz << endl;
		vector<double> qvec, sfvec;
		mc.getStrFct(0., 0.3, qz, qvec, sfvec);
    }*/

	ofstream outfile("makeCCDStrFct_So_050615.dat");
	//Create a structure factor by looping over all desired qz-values
    for (double qz = 0.05; qz <= 0.91; qz += 0.001) {
        cout << "qz= " << qz << endl;
		vector<double> qvec, sfvec;
//	mc.getStrFct(0, 0.25, qz, qvec, sfvec);
//		mc.getRotatedStrFct(0, 0.25, qz, qvec, sfvec);
		mc.getCCDStrFct(0, 0.21, qz, qvec, sfvec);
		for (vector<double>::size_type i = 0; i < qvec.size(); i++) {	
			outfile << qvec[i] << " " << qz << " " << sfvec[i] << endl;
    	}
    }
	outfile.close();
}

void makeUtable()
{
  ModelCalculator mc;
  mc.make_un_table();
} 

void makeNstable()
{
  ModelCalculator mc;
//  mc.make_ns_table();
} 

int main()
{
  //makeNstable();
  //makeUtable();
  //makeStrFct();
  makeStrFct2();
  //testUtable();
  //testStrFct();
  //testMosaicStrFct();
  //testRotatedStrFct();
  //testCCDStrFct();

  return 0;
}

