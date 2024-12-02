#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
//#include "alglib/src/stdafx.h"
//#include "alglib/src/interpolation.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <errno.h>
#include <stdexcept>
#include <algorithm>

#include "globalVariables.h"
#include "modelcalculator.h"

using std::string; using std::ifstream; using std::getline; 
using std::istringstream;
using std::vector; using std::cout; using std::endl;
using std::cerr;

//using namespace alglib;

#define PI 3.1415926535897932384626433832795
#define SMALLNUM 0.00000000001
#define WORKSPACE_SIZE 10000
#define KEY 6
#define ARB_SCALE 100000

ModelCalculator::ModelCalculator()
{
  //Turn off the error handler in gsl_integration_qag
  gsl_set_error_handler_off();
  //For qag numerical integration routine
  workspace = gsl_integration_workspace_alloc(WORKSPACE_SIZE);

  //Interpolation support functions
  spsumn.init(s_sumnWrapper, this);
  spStrFct.init(s_StrFctWrapper, this);
  spHr.init(s_HrWrapper, this);
  spRotated.init(s_rotatedWrapper, this);
  spMosaic.init(s_convolveMosaicWrapper, this);
  resetTOL();
}


/******************************************************************************
Cloning constructor that initializes necessary structures and copies
variables to provide a duplicate of a ModelCalculator object already
in use.
******************************************************************************/
ModelCalculator::ModelCalculator(const ModelCalculator& mc)
{
  // Copy necessary variables
  Kc = mc.Kc;
  B = mc.B;
  Kt = mc.Kt;
  at = mc.at;

Ls = mc.Ls; //new 2/26/15
divergeX = mc.divergeX; //new 2/26/15
divergeZ = mc.divergeZ; //new 2/26/15
teff = mc.teff; //new 3/30/15

  dspacing = mc.dspacing;
  T = mc.T;
  avgLr = mc.avgLr;
  avgMz = mc.avgMz;
  edisp = mc.edisp;
  mosaic = mc.mosaic;
  wavelength = mc.wavelength;
  rst = mc.rst; //mfe
  pixelSize = mc.pixelSize;
  bFWHM = mc.bFWHM;
  sdistance = mc.sdistance;
  beamSigma = mc.beamSigma;
  cutoff_r = mc.cutoff_r;
  cutoff_n = mc.cutoff_n;
  in = mc.in;
  out = mc.out;

  LrCutOff = mc.LrCutOff;
  MzCutOff = mc.MzCutOff;
  XiTx = mc.XiTx;
  XiTz = mc.XiTz;
  XiL = mc.XiL;

  workspace = gsl_integration_workspace_alloc(WORKSPACE_SIZE);

  utable = mc.utable;

  spsumn.init(s_sumnWrapper, this);
  spStrFct.init(s_StrFctWrapper, this);
  spHr.init(s_HrWrapper, this);
  spRotated.init(s_rotatedWrapper, this);
  spMosaic.init(s_convolveMosaicWrapper, this);
  resetTOL();

	// call the updaters
	XiTxChanged(); //new 2/26/15
	XiTzChanged(); //new 2/26/15
	XiLChanged(); //new 2/26/15
	LrCutOffChanged();
	MzCutOffChanged();
	set_beamSigma();
}


/******************************************************************************
Clean up resouces, should be called before the object is destructed
******************************************************************************/
void ModelCalculator::cleanup()
{
  gsl_integration_workspace_free(workspace);
}


/******************************************************************************
set tolerance for interpolation
******************************************************************************/
void ModelCalculator::resetTOL()
{
  spStrFct.setol(g_spStrFctAbserr, g_spStrFctRelerr, g_spStrFctMindx, g_spStrFctMaxdx);
  spsumn.setol(g_spsumnAbserr, g_spsumnRelerr, g_spsumnMindx, g_spsumnMaxdx);
  spHr.setol(g_spHrAbserr, g_spHrRelerr, g_spHrMindx, g_spHrMaxdx);
  spRotated.setol(g_spRotatedAbserr, g_spRotatedRelerr, g_spRotatedMindx, g_spRotatedMaxdx);
  spMosaic.setol(g_spMosaicAbserr, g_spMosaicRelerr, g_spMosaicMindx, g_spMosaicMaxdx);
}


/******************************************************************************
Initiate building the interpolants for the pure structure factor, mosaic-
convolved structure factor, and rotated structure factor.
******************************************************************************/
void ModelCalculator::QxSlice(double qxlow, double qxhigh)
{
	//cout << "Kc: " << Kc << endl;
	//cout << "Kt: " << Kt << endl;
	//cout << "at: " << at << endl;
	//sleep(1000000);

if (qxlow < 0 || qxhigh < 0)
		throw domain_error("ERROR: qxlow and qxhigh must both take positive values");
  double sig_qx = get_sigx(qxhigh, sdistance, beamSigma, Ls);
//  sig_qx = beamSigma; /* uncomment to revert to older version */
  buildInterpForRotatedStrFct(0, qxhigh + 20*sig_qx);

if (qxhigh + 2*sig_qx > 0.5) {
	cout << "qxhigh = " << qxhigh << endl;
	cout << "sdistance = " << sdistance << endl;
	cout << "beamSigma = " << beamSigma << endl;
	cout << "Ls = " << Ls << endl;
	cout << "sig_qx = " << sig_qx << endl;
	cout << "Interp qxhigh = " << qxhigh + 20*sig_qx << endl;
}

}


/******************************************************************************
In order to build rotated structure factor interpolant, the program needs
mosaic-convolved structure factor. This function decides the necessary
range of qr for mosaic-convolved structure factor. The upper limit on
the range is slightly larger than qx high limit for rotated structure factor
because qr >= qx for any value of qy.
******************************************************************************/
void ModelCalculator::buildInterpForRotatedStrFct(double qxlow, double qxhigh)
{
  double omegaTop = get_omegaC(qz, wavelength);
  double qylowerLimit = get_qyLimits(qxhigh, qz, wavelength, 0);
  double qyupperLimit = get_qyLimits(qxhigh, qz, wavelength, omegaTop);
  double qy = max(fabs(qylowerLimit),fabs(qyupperLimit));

  double qrhigh = sqrt(qxhigh*qxhigh + qy*qy);
  buildInterpForMosaicStrFct(qxlow, qrhigh);
  spRotated.findPoints( log(qxlow+SMALLNUM), log(qxhigh+SMALLNUM) );
}


void ModelCalculator::buildInterpForMosaicStrFct(double qrlow, double qrhigh)
{
  buildInterpForStrFct(0, qrhigh+0.); // temporary? change as of 2/25/2014
  qrUpperLimit = qrhigh + 0.; // temporary? change as of 2/25/2014
  spMosaic.findPoints( log(qrlow+SMALLNUM), log(qrhigh+SMALLNUM) );
}


void ModelCalculator::buildInterpForStrFct(double qrlow, double qrhigh)
{
  spStrFct.findPoints( log(qrlow+SMALLNUM), log(qrhigh+SMALLNUM) );
}


/******************************************************************************
FuncLmdif::modelCalcThread calls this function.
This function returns the interpolated value of the rotated structure factor
at the input qx. The interpolant returns the log of a rotated structure factor, 
so this function takes the exponential of the interpolated value.

This function seems strange. Why should the non-beam convoluted structure factor
be returned? (OK if beam size on detector is very small. Is 0.05 sufficiently small?)
Maybe should always use beamConvolutedStrFct instead of 
exp( spRotated.val(log(qx + SMALLNUM)) . . .
******************************************************************************/
double ModelCalculator::getCCDStrFct(double qx)
{
  // 0.05 hard coded for now. For some test cases, this worked well.
	if (bFWHM > 0.05) {
		return beamConvolutedStrFct(qx);
	} else {
		return exp( spRotated.val(log(qx + SMALLNUM)) );
	}
}


/******************************************************************************
This function takes care of beam FWHM convolution with rotated structure
factor.
******************************************************************************/
double ModelCalculator::beamConvolutedStrFct(double qx)
{
	currQx = qx;
  double result, abserr;
  gsl_function F;
  F.function = &s_convIntegrand;
  F.params = this;

  // set integration limits
  double sig = get_sigx(qx, sdistance, beamSigma, Ls); //new 2/4/2015
//  sig = beamSigma; /* uncomment to revert to older version */
  double lowerLimit = qx - 10*sig; 
  double upperLimit = qx + 10*sig;

  gsl_integration_qag(&F, lowerLimit, upperLimit, g_epsabs, g_epsrel,
                      WORKSPACE_SIZE, KEY, workspace, &result, &abserr);
  return result;
}


/******************************************************************************
Return the integrand for the integration that performs beam convolution
******************************************************************************/
double ModelCalculator::s_convIntegrand(double qx, void *params)
{
	ModelCalculator *p = (ModelCalculator *)params;
	return exp( p->spRotated.val(log(fabs(qx) + SMALLNUM)) )
	       * p->beamProfile(qx - p->currQx);
}

/******************************************************************************
Calculates omega_C as defined in MJ thesis
******************************************************************************/
double get_omegaC(double qz, double wavelength)
{
	double k = 2*PI/wavelength;
	return atan(qz/k);
}

/******************************************************************************
Calculates an approximation to standard deviation of the beam height
******************************************************************************/
double get_sigx(double qx, double sdistance, double beamSigma, double Ls)
{  
  return beamSigma + Ls*qx/2.3548/sdistance;
}

/******************************************************************************
Return the beam profile in q-space. It is assumed to be Gaussian.
******************************************************************************/
double ModelCalculator::beamProfile(double qx)
{
	double sig = get_sigx(qx, sdistance, beamSigma, Ls);
//	sig = beamSigma; /* uncomment to revert to older version */
	return exp(-qx * qx / 2 / sig / sig) / sqrt(2 * PI) / sig;
}


/******************************************************************************
Static function that returns the logarithmic value of the rotated structure 
factor at the input log(qx). (Note that it is base e, not base 10)
It uses the absolute value of qx. This is to deal with qx = 0 point, where
floating point arithmetic could potentially yield a very slightly negative 
value that represents zero in doulbe precision floating point.
******************************************************************************/
double ModelCalculator::s_rotatedWrapper(double logqx, void *ptr)
{
	ModelCalculator *p = (ModelCalculator *)ptr;
	double tmpQx = fabs(exp(logqx) - SMALLNUM);
	return log(p->rotated(tmpQx));
}

/******************************************************************************
Calculate qy limits as defined in MJ thesis
******************************************************************************/
double get_qyLimits(double qx, double qz, double wavelength, double omega)
{
	double k = 2*PI/wavelength;	
	return -(0.5/k)*(qx*qx + qz*qz) + qz*omega;
}

/******************************************************************************
This function performs the integration over the structure factor along qy 
direction at the input qx value
The lower and upper integration limits depend on qx, so their values get 
retrieved by ModelCalculator::get_qyLimits()
*******************************************************************************/
double ModelCalculator::rotated(double qx)
{
	currQx = qx;
	double result, abserr;
	double omegaTop = get_omegaC(qz, wavelength);
	gsl_function F;
	F.function = &s_qyIntegrandWrapper;
	F.params = this;
	double lowerLimit = get_qyLimits(qx, qz, wavelength, 0);
//	double omegaExpMax = 7.0*PI/180.0;
//if (omegaTop > omegaExpMax)	omegaTop = omegaExpMax;
	double upperLimit = get_qyLimits(qx, qz, wavelength, omegaTop);
//cout << "new upperLimit = " << upperLimit2 << endl;
	gsl_integration_qag(&F, lowerLimit, upperLimit, g_epsabs, g_epsrel,
	                    WORKSPACE_SIZE, KEY, workspace, &result, &abserr);
	return result;
}

/******************************************************************************
This function returns the integrand for the qy integration, which is
the interpolated value of the mosaic convolved structure factor at given qx, 
qy, and qz.
The interpolation over the structure factor is done in (qr,qz) tuple.
Calculates qr corresponding to the input qy and currQx.
This wrapper is necessary for GSL qag.
******************************************************************************/
double ModelCalculator::s_qyIntegrandWrapper(double qy, void *ptr)
{
	ModelCalculator *p = (ModelCalculator *)ptr;
	double qr = sqrt(p->currQx*p->currQx + qy*qy);
	double omega_approx = p->get_omegaApprox(p->currQx,qy,p->qz,p->wavelength); //NEW 3/27/15
	double absorp = p->get_absorp(omega_approx,p->qz,p->wavelength); //NEW 3/27/15
// absorp = 1; //reverts to no absorption correction version

//cout << "qz= " << p->qz << " omega= " << omega_approx << " absorp= " << absorp << endl;
	return p->getMosaicStrFct(qr)*absorp;
}

/******************************************************************************
Calculates approximation to incident angle : defined in MJ thesis
******************************************************************************/
double ModelCalculator::get_omegaApprox(double qx, double qy, double qz, double wavelength)
{	
	double temp = wavelength/4/PI;
	return (1/qz)*(qy + temp*(qx*qx + qz*qz));
}

/******************************************************************************
Calculates absorption correction : defined in MJ thesis
teff is the ratio of the X-ray absorption length to the sample thickness
(it is input from the GUI)
******************************************************************************/
double ModelCalculator::get_absorp(double omega, double qz, double wavelength)
{	
	double k = 2*PI/wavelength;
	double g = 1/sin(omega);
	g += 1/sin(atan(qz/k)-omega);

	return (teff)*(1 - exp(-g/teff))/g;
}

/******************************************************************************
Return the structure factor after mosaic spread is treated.
If mosaic > 0.001, it will return an interpolated value of the structure 
factor with mosaic convolution.
If mosaic is less than 0.001, it will return an interpolated value of the 
structure factor without any additional operation done on itself.
******************************************************************************/
double ModelCalculator::getMosaicStrFct(double qr)
{
  // Takes the absolute value of qr to avoid a very small negative value 
  // representing zero
  const double negligible_mosaic = 0.001 * PI / 180;
  if (mosaic > negligible_mosaic) {
    return exp(spMosaic.val(log(fabs(qr)+SMALLNUM)));
  } else {
    return exp(spStrFct.val(log(fabs(qr)+SMALLNUM)));
  }
}


/****************************************************************************** 
Wrapper needed for FunSupport object, spMosaic
******************************************************************************/
double ModelCalculator::s_convolveMosaicWrapper(double logqr, void *ptr)
{
  ModelCalculator *p = (ModelCalculator *)ptr;
  double tmpQr = fabs(exp(logqr) - SMALLNUM);
  return log(p->convolveMosaic(tmpQr));
}


/******************************************************************************
This function performs convolution of the structure factor with the mosaic
distribution.
******************************************************************************/
double ModelCalculator::convolveMosaic(double qr)
{
  currQr = qr;
  double result, abserr;
  gsl_function F;
  F.function = &s_mosaicIntegrandWrapper;
  F.params = this;
  double lowerLimit = 0;
  double upperLimit = qrUpperLimit;
  gsl_integration_qag(&F, lowerLimit, upperLimit, g_epsabs, g_epsrel,
                      WORKSPACE_SIZE, KEY, workspace, &result, &abserr);
  return result;
}


/******************************************************************************
This function returns the integrand of the mosaic convolution. 
******************************************************************************/
double ModelCalculator::s_mosaicIntegrandWrapper(double qr, void *ptr)
{
  ModelCalculator *p = (ModelCalculator *)ptr;
	return exp(p->spStrFct.val(log(fabs(qr)+SMALLNUM))) * 
	       p->mosaicDist(p->currQr-qr);
}


/******************************************************************************
Mosaic spread distribution. Mosaic angle is approximated by qr/qz, which is
good for small angle.
******************************************************************************/
double ModelCalculator::mosaicDist(double qr)
{
  return 2 * mosaic / PI / (4*qr*qr + mosaic*mosaic);
}


/******************************************************************************
Static wrapper for structure factor along qr direction.
The input is log(qr), base e.
Returns log(S(qr,qz)).
This wrapper is necessary for FunSupport object.
******************************************************************************/
double ModelCalculator::s_StrFctWrapper(double logqr, void *ptr)
{
  ModelCalculator *p = (ModelCalculator *)ptr;
  double tmpQr = fabs(exp(logqr) - SMALLNUM);
  double ret = p->StrFct(tmpQr);
  if (ret < 0) {
    cout << "\n\tNegative structure factor was obtained at" << endl;
    cout << "\tqr: " << p->currQr << " qz: " << p->qz << " Value: " << ret << endl;
    cout << "\t" << p->Kc << " " << p->B << " " << p->LrCutOff << " " << 
					p->MzCutOff << " " << p->dspacing << " " << p->T << " " << 
					p->wavelength << endl;
	cout << "\t" << p->XiTx << " " << p->XiTz << " " << p->XiL << endl;

	double hop = 0.000001;
	double tmp1Qr = tmpQr-1*hop;
	double tmp1 = p->StrFct(tmp1Qr);
  	int i = 1;
	while ( (tmp1 < 0 && i < 1000) ) {
	  i++;	
	  tmp1 = p->StrFct(tmpQr-i*hop);
	  tmp1Qr = tmpQr-i*hop;
	}

	double tmp2Qr = tmpQr+1*hop;
	double tmp2 = p->StrFct(tmp2Qr);
	int j = 1;
	while ( (tmp2 < 0 && j < 1000) ) {
	  j++;	
	  tmp2 = p->StrFct(tmpQr+j*hop);
	  tmp2Qr = tmpQr+j*hop;
	}

	if ( (tmp1 > 0) && (tmp2 > 0) ) {
		cout << "\tStructure factor was calculated " << i+j << " extra times." << endl;
		cout << "\ttmp1: " << tmp1 << " tmp2: " << tmp2 << endl;
		cout << "\ttmp1Qr: " << tmp1Qr << " tmp2Qr: " << tmp2Qr << endl;
        ret = tmp1 + (tmp2-tmp1)*((tmpQr-tmp1Qr)/(tmp2Qr-tmp1Qr)); //linear interpolation
		cout << "\tpositive result --> " << ret << endl;
    } else {
		cout << "\tStructure factor was calculated " << i+j << " extra times." << endl;
		cout << "\tPositive structure factor was not found." << endl;
		// The structure factor at sufficiently large tmpQr is only needed to complete
		// the beam convolution. If a positive structure factor can not be determined,
		// just set structure factor to something small. 
		ret = 0.01;
	cout << "\tret -> " << ret << endl;
	}
  }
  return log(ret);
}

/******************************************************************************
This function calculates the structure factor as a function of qr for a fixed 
value of qz specified by qz variable
******************************************************************************/
double ModelCalculator::StrFct(double qr)
{
  currQr = qr;
//cout << "qr= " << currQr << endl;
  double result, abserr;
  gsl_function F;
  F.function = &s_Fr;
  F.params = this;
  double lowerLimit = 0;
  double upperLimit = cutoff_r;
//cout << "cutoff_r:" << cutoff_r << endl;
   gsl_integration_qag(&F, lowerLimit, upperLimit, g_epsabs, g_epsrel,
                      WORKSPACE_SIZE, KEY, workspace, &result, &abserr);

//cout << "Using qag:                                " << result << endl;

//result=-1;
//////////////////////////////////////////////////////////////
// Function calculates integral by summing many integrals
// with bounds defined by the zeroes of the integrand (J0).
// Usefull for troubleshooting error of numerically integrating
// oscillatory function.
//////////////////////////////////////////////////////////////
/*if (result < 0) {
cout << "Using qag:                                " << result << endl;
    double sum = 0, J0zero_pos;
	int nth_zero = 1;
	J0zero_pos = gsl_sf_bessel_zero_J0(nth_zero)/currQr;
    while (J0zero_pos < lowerLimit) {
      nth_zero += 1;
      J0zero_pos = gsl_sf_bessel_zero_J0(nth_zero)/currQr;
    }
    double ul = 0;
    while (ul < upperLimit) {
	  J0zero_pos = gsl_sf_bessel_zero_J0(nth_zero)/currQr;
	  if (upperLimit < J0zero_pos) {
	    ul = upperLimit;
	  } else {
        ul = J0zero_pos;
      }
//cout << " ll: " << lowerLimit << " ul: " << ul << endl;
      gsl_integration_qag(&F, lowerLimit, ul, g_epsabs, g_epsrel,
                      WORKSPACE_SIZE, KEY, workspace, &result, &abserr);
      lowerLimit = ul;
	  nth_zero += 1;
	  sum += result;
//cout << "ul: " << ul << " result: " << result << " sum: " << sum << endl;
//usleep(100000);
    }
    result = sum;
cout << "Summing integrals between bessel zeroes: " << result << endl;
//usleep(1000000);
  }*/

//usleep(1000000);

  return result;
}


/******************************************************************************
This calculates f_2(r,qr) integrand
It is equal to r * H_r(r) * J_0(qr*r) * f_1(r)
Two functions are their corresponding interpolation functions,
H_r(r) and f_1(r) are both one dimensional interpolation using FunSupport class
J_0(qr*r) is the zeroth order Bessel function of the first kind
******************************************************************************/
double ModelCalculator::s_Fr(double r, void *ptr)
{
// Original 2/25/2014
// Negative structure factor issue is not related to using
// bessel_J0 as opposed to gsl_sf_bessel_J0
	ModelCalculator *p = (ModelCalculator *)ptr;
  return r * p->spHr.val(r) * bessel_J0(p->currQr*r) *
         p->spsumn.val(log(r+SMALLNUM)) / ARB_SCALE;

/*  return r * p->spHr.val(r) * bessel_J0(p->currQr*r) *
         p->sumn(r) / ARB_SCALE;*/

/*	ModelCalculator *p = (ModelCalculator *)ptr;
  return r * p->spHr.val(r) * gsl_sf_bessel_J0(p->currQr*r) *
         p->spsumn.val(log(r+SMALLNUM)) / ARB_SCALE;*/
}


/******************************************************************************
static function that wraps ModelCalculator::sumn, so that it can be passed to
FunSupport::init
FunSupport evaluates this function at constant log(r) step. This function
converts log(r) to r. SMALLNUM is used to avoid Nan coming from exp(log(0)).
When compiled under g++, as of 7/3/2013, log(0) yields -inf, but exp(log(0)) 
yeilds Nan, instead of 0.
SMALLNUM is added to r, like log(r+SMALLNUM), when an interpolated value gets
retrieved.
******************************************************************************/
double ModelCalculator::s_sumnWrapper(double logr, void *ptr)
{
  ModelCalculator *p = (ModelCalculator *)ptr;
  double tmpR = fabs(exp(logr) - SMALLNUM);
  return p->sumn(tmpR);
}

/******************************************************************************
Calculate height difference function for testing.
******************************************************************************/
double ModelCalculator::get_un_value(int n, double r)
{
  double un;
  double xi = 1e8 * pow(Kc/B, 0.25);
  double xit_2 = 1e16 * Kc/Kt; //prefactor to convert from cm^2 to angstrom^2
  double ell = 2 * xit_2 / (xi * xi);
  double eta = 2.16872 * (T+273.15) / sqrt(Kc*B) / dspacing /dspacing;
  un = utable.getHeightDiffFunc(n, r, xi, ell,rst, eta, dspacing, at);
  return un;
}

/******************************************************************************
Write u-table file.
******************************************************************************/
void ModelCalculator::make_un_table()
{
  utable.writeUtableFile(rst);
}

/******************************************************************************
Calculate the sum over layer index, n.
It is equal to f_1(r)
******************************************************************************/
double ModelCalculator::sumn(double r)
{
  double sum = 0;
  double un, t1, t2;
  double qz2 = qz * qz;
  double xi = 1e8 * pow(Kc/B, 0.25);
  double xit_2 = 1e16 * Kc/Kt; //prefactor to convert from cm^2 to angstrom^2
  double ell = 2. * xit_2 / (xi * xi);
  double eta = 2.16872 * (T+273.15) / sqrt(Kc*B) / dspacing /dspacing;

  double sig = get_sigz(qz, wavelength, sdistance, Ls); /* new version */
  double sig2 = sig*sig; /* new version */
//sig2 = (qz*edisp) * (qz*edisp); /* uncomment to revert to older version */

  double max_n = cutoff_n;
if ( onlyZERO == 1 ) {
max_n = 1; /* calculates using only n=0 */
}

//max_n = 20 * avgMz; /* uncomment to revert to older version */

  //Calculate n=0 term separately because it has 1/2 weight relatively to n>0 terms
  int n = 0;
  un = utable.getHeightDiffFunc(n, r, xi, ell,rst, eta, dspacing, at);
//cout << "r= " << r << "; at= " << at << "; un= " << un << endl;
//usleep(10000);
  t1 = 1 + un*sig2;
  t2 = n * dspacing;
  sum += Hz(n) * sqrt(1/t1) * exp( -0.5*(qz2*un+t2*t2*sig2)/t1 ) * cos(t2*qz/t1);
  for(n = 1; n < max_n; n++) {
    un = utable.getHeightDiffFunc(n, r, xi, ell,rst, eta, dspacing, at);
    t1 = 1 + un*sig2;
    t2 = n * dspacing;
    sum += 2 * Hz(n) * sqrt(1/t1) * exp( -0.5*(qz2*un+t2*t2*sig2)/t1 ) * cos(t2*qz/t1);
  }
//cout << "r=" << r << "; sum= " << sum << endl;
  return sum;
}

/******************************************************************************
Calculates an approximation to standard deviation of the beam height
******************************************************************************/
double get_sigz(double qz, double wavelength, double sdistance, double Ls)
{
  double k = 2*PI/wavelength;
  double wC = get_omegaC(qz,wavelength);
  // bFWHM / 2.3548 is equal to beam Gaussian sigma
  
  return (Ls/sdistance/2.3548)*(qz - k*wC/2);
}

/******************************************************************************
Static function that wrapps ModelCalculator::Hr
******************************************************************************/
double ModelCalculator::s_HrWrapper(double r, void *ptr)
{
  ModelCalculator *p = (ModelCalculator *)ptr;
  return p->Hr(r);
}


/******************************************************************************
Effective finite size factor in r direction, H_r(r)
The upper limit is chosen to be 10 times the cutoff in r integration.
The factor of 10 is somewhat arbitrary.
******************************************************************************/
double ModelCalculator::Hr(double r)
{
  curr_r = r;
  double result, abserr;
  gsl_function F;
  F.function = &s_HrIntegrand;
  F.params = this;
  double lowerLimit = r;
  double upperLimit = 10 * cutoff_r;
  gsl_integration_qag(&F, lowerLimit, upperLimit, g_epsabs, g_epsrel,
                      WORKSPACE_SIZE, KEY, workspace, &result, &abserr);
  return PI * result;
}


/******************************************************************************
The integrand in H_r(r) integration
******************************************************************************/
double ModelCalculator::s_HrIntegrand(double Lr, void *params)
{
  ModelCalculator *p = (ModelCalculator *)params;
  return Pr(Lr, p->LrCutOff) * Lr * Lr * finiteSizeFactor(p->curr_r/Lr);
}


/******************************************************************************
In-plane domain size distribution
Can be changed to another function easily
******************************************************************************/
double Pr(double Lr, double LrCutOff)
{

  double LrStar = LrCutOff;
// LrStar = avgLr; /* uncomment to revert to older version */

/* Below is code for gaussian in-plane distribution */
//  double sigma_r = LrStar/3.0; //could choose a different value instead of "3.0"
//  double  a = (Lr - LrStar) / sqrt(2.0) / sigma_r;
//  return exp(-a*a) / sigma_r;

/* Below is code for exponential in-plane distribution */
  return exp(-Lr / LrStar) / LrStar; /* exponential domain distribution */
}


/******************************************************************************
Finite size factor in r direction, F_r(r/Lr)
******************************************************************************/
double finiteSizeFactor(double x)
{
  if (x > 1) {
    return 0;
  } else {
    return acos(x) - x*sqrt(1-x*x);
  }
}


/******************************************************************************
Effective finite size factor in z direction, H_z(nD)
For either exponential or gaussian out-of-plane domain distributions
the finite size factor integral can be computed analytically.
******************************************************************************/
double ModelCalculator::Hz(int n)
{
	double MzStar = MzCutOff; //new 2/4/2015
// MzStar = avgMz; /* uncomment to revert to older version */
//	double sigma_n = MzStar/3.0; //could choose a different value instead of "3.0"

/* Below is code for gaussian out-of-plane finite size factor */
//	double a = (n-MzStar) / sqrt(2.0) / sigma_n; 
//	double ret = (MzStar-n)*sqrt(PI)*gsl_sf_erfc(a);
//	ret += sqrt(2.0)*sigma_n*exp(-a*a);
//		return sqrt(2.0) * ret / 2.0 / MzStar;

/* Below is code for exponential out-of-plane finite size factor */
	return exp(-n/MzStar); /* exponential domain distribution */
}

/******************************************************************************
update value of onlyZERO
If onlyZERO == 1 then S_0(q) is calculated
otherwise S(q) is calculated.
See MJ thesis for definition of S_0(q)
******************************************************************************/
void ModelCalculator::set_onlyZERO(int _onlyZERO)
{
  onlyZERO = _onlyZERO;
}


/******************************************************************************
given a slice with qz = tqz, update necessary temporary variables
and tables.
******************************************************************************/
void ModelCalculator::setSliceParameter(double _qz)
{
  qz = _qz;
  MzCutOffChanged();
  spsumn.findPoints( log(0+SMALLNUM), log(10*cutoff_r+SMALLNUM) );
}

/* appears to be vestigial 7/28/15 . . . */
void ModelCalculator::KcBDTChanged()
{
  //double lambda = 1e16 * sqrt(Kc/B) / dspacing;
  //double eta = 2.17 * (T+273.15) / sqrt(Kc*B) / dspacing /dspacing;
  //in = sqrt(utable.lambda * utable.D / lambda / dspacing);
  //out = eta * dspacing * dspacing / utable.D / utable.D / utable.eta;
}


/******************************************************************************
update domain size parameters
When they are modified, the cutoffs and the interpolant must also be modified
******************************************************************************/
void ModelCalculator::LrCutOffChanged()
{
	LrCutOff = min(avgLr,min(XiTx,XiTz));
	cutoff_r = 20 * LrCutOff; // was 20
	spHr.findPoints(0, 10*cutoff_r);
}

void ModelCalculator::MzCutOffChanged()
{
	MzCutOff = min(avgMz,min(XiTz/dspacing,2*PI*XiL/qz/wavelength/dspacing));
	cutoff_n = 20 * floor(MzCutOff); // was 20
}

// 2/26/25 several new updating functions to keep track of beam coherence lengths
void ModelCalculator::XiTxChanged()
{
	XiTx = wavelength / 2 / divergeX;
	LrCutOffChanged();
}

void ModelCalculator::XiTzChanged()
{
	XiTz = wavelength / 2 / divergeZ;
	LrCutOffChanged();
	MzCutOffChanged();
}

void ModelCalculator::XiLChanged()
{
	XiL = wavelength / 2. / edisp;
	MzCutOffChanged();
}

/******************************************************************************
Converts beam FWHM in pixels to beam Gaussian sigma in q-space units
******************************************************************************/
void ModelCalculator::set_beamSigma()
{
	double deltaQPerPixel = 2 * PI * pixelSize / wavelength / sdistance;
	// bFWHM / 2.3548 is equal to beam Gaussian sigma
	beamSigma = deltaQPerPixel * bFWHM / 2.3548;
}


/* mosaic is in radian, but the program takes the input in degrees */
void ModelCalculator::set_mosaic(double a)
{
  mosaic = PI * a / 180;
}


/******************************************************************************
This function sets a model parameter to an input value.
A model parameter is specified by std::string name.
See stringToVarEnum function for a list of strings that can be
passed as an input name.

value: an input value
name: an input name specifying which parameter to be set to an input value
******************************************************************************/
void ModelCalculator::setModelParameter(double value, const char *_name)
{
	string name(_name);
	setModelParameter(value, name);
}


void ModelCalculator::setModelParameter(double value, const string& name)
{
  setpara(value, stringToVarEnum(name));
}


/******************************************************************************
set parameters for Sccd calculation

value : the value to be set
idx: the index for the corresponding parameter
******************************************************************************/
void ModelCalculator::setpara(double value, int idx)
{
  switch(idx) {
		case Var_Kc: Kc = value; break;
		case Var_B: B = value; break;
		case Var_D: dspacing = value; break;
		case Var_T: T = value; break;
		case Var_avgLr: avgLr = value; LrCutOffChanged(); break;
		case Var_avgMz: avgMz = value; MzCutOffChanged(); break;
		case Var_mosaic: set_mosaic(value); break;
		case Var_edisp: edisp = value; XiLChanged(); break;
		case Var_wavelength: wavelength = value; XiTxChanged(); XiTzChanged();
							 XiLChanged(); set_beamSigma(); break; //new 2/26/15
		case Var_pixelSize: pixelSize = value; set_beamSigma(); break;
		case Var_bFWHM: bFWHM = value; set_beamSigma(); break;
		case Var_sdistance: sdistance = value; set_beamSigma(); break;
    case Var_Kt: Kt = value; break; //new
    case Var_at: at = value; break; //new
	case Var_Ls: Ls = value; break; //new 2/16/15
    case Var_divergeX: divergeX = value; XiTxChanged(); break; //new 2/26/15
    case Var_divergeZ: divergeZ = value; XiTzChanged(); break; //new 2/26/15
    case Var_teff: teff = value; break; //new 3/30/15
    case Var_rst: rst = value; break; // mfe
		default: ;
  }

}


/******************************************************************************
make the ModelCalculator in sync with parameter set p
******************************************************************************/
void ModelCalculator::paraSet(Para *p)
{
  Kc = p->Kc;
  B = p->B;
  Kt = p->Kt; //new
  at = p->at; //new
  Ls = p->Ls; //new

  divergeX = p->divergeX; //new
  divergeZ = p->divergeZ; //new

  teff = p->teff; //new

  dspacing = p->D;
  T = p->T;
  avgLr = p->avgLr;
  avgMz = p->avgMz;
  edisp = p->edisp;
  rst = p->rst; // mfe
  wavelength = p->setup.wavelength;
  pixelSize = p->setup.pz;
  bFWHM = p->bFWHM;
  sdistance = p->setup.s;
  XiTxChanged(); //new
  XiTzChanged(); //new
  XiLChanged(); //new
  LrCutOffChanged();
  MzCutOffChanged();
  set_beamSigma();
  set_mosaic(p->mosaic);
}


void ModelCalculator::get_spStrFctPoints(vector<double>& x, vector<double>& y)
{
	spStrFct.getPoints(x, y);
}


void ModelCalculator::get_spsumnPoints(vector<double>& x, vector<double>& y)
{
	spsumn.getPoints(x, y);
}


void ModelCalculator::get_spHrPoints(vector<double>& x, vector<double>& y)
{
	spHr.getPoints(x, y);
}


void ModelCalculator::get_spRotatedPoints(vector<double>& x,vector<double>& y)
{
	spRotated.getPoints(x, y);
}


void ModelCalculator::getCCDStrFct(double qxlow, double qxhigh, double qz,
                                   vector<double>& qxv, vector<double>& sfv)
{
  setSliceParameter(qz);
  QxSlice(qxlow, qxhigh);
  
  for (double qx = qxlow; qx < qxhigh; ) {
    qxv.push_back(qx);
    sfv.push_back(getCCDStrFct(qx));
    qx += 0.001;
  }
}


void ModelCalculator::getRotatedStrFct(double qxlow, double qxhigh, double _qz, 
                                       vector<double>& qxvec, vector<double>& rsf)
{
	setSliceParameter(_qz);
	buildInterpForRotatedStrFct(qxlow, qxhigh);

	for (double qx = qxlow; qx < qxhigh; ) {
		qxvec.push_back(qx);
		rsf.push_back( exp( spRotated.val(log(qx+SMALLNUM)) ) );
		qx += 0.001;
	}
}


void ModelCalculator::getMosaicStrFct(double qrlow, double qrhigh, double _qz,
                                      vector<double>& qr, vector<double>& sf)
{
	if (qrlow < 0 || qrhigh < 0)
		throw domain_error("qrlow and qrhigh must both take positive values");
	setSliceParameter(_qz);
  buildInterpForMosaicStrFct(qrlow, qrhigh);

	for (double q = qrlow; q < qrhigh;) {
		qr.push_back(q);
		sf.push_back(getMosaicStrFct(q));
		q += 0.001;
	}
}


void ModelCalculator::getStrFct(double qrlow, double qrhigh, double _qz, 
                                vector<double>& qrvec, vector<double>& sf)
{
	if (qrlow < 0 || qrhigh < 0)
		throw domain_error("qrlow and qrhigh must both take positive values");
	setSliceParameter(_qz);
	buildInterpForStrFct(qrlow, qrhigh);

	for (double qr = qrlow; qr < qrhigh;) {
		qrvec.push_back(qr);
		sf.push_back( exp( spStrFct.val(log(qr + SMALLNUM)) ) );
		qr += 0.0005; //change back to 0.001 6/19/2014
/* above line sets qr step size for stucture factor test calculations */
	}
}


void ModelCalculator::getSumn(double rlow, double rhigh, double _qz, 
                              vector<double>& rvec, vector<double>& zvec)
{
	setSliceParameter(_qz);
	for (double r = rlow; r < rhigh; ) {
		rvec.push_back(r);
		zvec.push_back( exp( spsumn.val(log(r+SMALLNUM)) ) );
		r += 10;
	}
}


void ModelCalculator::getHr(double rlow, double rhigh, 
                            vector<double>& rvec, vector<double>& hvec)
{
	LrCutOffChanged();
	for (double r = rlow; r < rhigh; ) {
		rvec.push_back(r);
		hvec.push_back(spHr.val(r));
		r += 100;
	}
}


void saveDoubleColumns(vector<double>& xvec, vector<double>& yvec, 
                       const char *filename)
{
  cout << "xvec.size(): " << xvec.size() << " yvec.size(): " << yvec.size() 
       << endl;
  if (xvec.size() != yvec.size()) {
    cerr << "\nThe size of input vectors must be identical for " 
         << "saveDoubleColumns function to work\n" << endl;
    return;
  }
  ofstream myfile;
  myfile.open(filename);
  typedef vector<double>::size_type vec_sz;
  for (vec_sz i = 0; i < xvec.size(); i++) {
    myfile << xvec[i] << " " << yvec[i] << endl;
  }
}
