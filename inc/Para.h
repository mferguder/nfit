#ifndef GUARD_PARA_H
#define GUARD_PARA_H

#include <tcl.h>
#include "tvds.h"
#include "nfit.h"
#include "globalVariables.h"

//#define TOTAL_NUM_OF_FIT_PARAMS 15
enum Var stringToVarEnum(const std::string &);


/******************************************************************************
parameter class
******************************************************************************/
struct Para {
public:
  // The followings are the fitting parameters ////////////////////////////////
  double Kc; // bending modulus
  double B; // bulk modulus
  double avgLr; // average domain size in-plane
  double avgMz; // average domain size out-of-plane
  double D; // D-spacing
  double mosaic; // mosaic spread
  double edisp; // energy dispersion and beam divergence
  double bFWHM; // beam full width half maximum in pixel
  TvDs2 setup;   // piece borrowed from DS
  double T; // temperature
  double Kt; // tilt modulus (in ergs / cm^2) NEW
  double at; // smallest real-space in-plane length scale (in Angstroms) NEW
  double Ls; // width of sample along beam direction, typically 5 (in mm) 2/16/15
  double divergeX; // beam divergence in x-direction, typically 0.0001 rad 2/26/15
  double divergeZ; // beam divergence in z-direction, typically 0.0001 rad 2/26/15
  double teff; // absorption length divided by sample thickness, ~260 3/30/15
  double rst;
  /////////////////////////////////////////////////////////////////////////////
  
  double beamSigma; // beamFWHM in q-space (inverse Angstrom)
  double lambda; // fluctuation parameter
  double eta; // fluctuation parameter

  std::vector<double> epsfcn; // used in lmdif
  std::vector<int> idx; // the index of free parameters
  std::vector<double> xp; // initial values of free parameters
  //bool isFree[g_totalNumOfFitParams]; // This line does not work.
  bool isFree[TOTAL_NUM_OF_FIT_PARAMS];
  int nfit; // the number of free parameters for a current fit

  void set_index(const std::vector<bool>&);
  void set_xp();
  void set_epsfcn();
  void setValue(double, int);
  void setValue(double, const char *);
  void setValue(double, const std::string&);
  void setLambdaEta();
  double getValue(int);
  void set_beamSigma();
  void set_beamSigma(double x) {beamSigma = x;}
  double get_beamSigma() {return beamSigma;}
};

#endif
