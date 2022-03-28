#include <tcl.h>
#include <iostream>

#include "Para.h"
#include "nfit.h"

#define PI 3.1415926535897932384626433832795

using std::cout; using std::endl; using std::vector; using std::string;


/****************************************************************************************
This function sets a model parameter to an input value.
A model parameter is specified by the name string.
See ModelCalculator::stringToVarEnum for a list of strings that can be
passed as an input name.

When a new parameter needs to be added to the model, add a new case,
add a new if statement in ModelCalculator::stringToVarEnum, and
add a new enum value in enum Var

value: an input value
name: an input name specifying which parameter to be set to an input value
****************************************************************************************/
void Para::setValue(double value, const char *_name)
{
	string name(_name);
	setValue(value, name);
}


void Para::setValue(double value, const string& name)
{
  setValue(value, stringToVarEnum(name));
}


/****************************************************************************************
an accessor function
Given an input index, a corresponding variable gets set to an input value
When either Kc, B, D, or T gets modified, lambda and eta get recalculated.
This function accepts int instead of enum itself so that a user of this object
does not have to define a separate enum.
****************************************************************************************/
void Para::setValue(double value, int index)
{
  switch (index) {
    case Var_Kc:                       Kc = value;  setLambdaEta(); break;
    case Var_B:                         B = value;  setLambdaEta(); break;
    case Var_avgLr:                 avgLr = value;                  break;
    case Var_avgMz:                 avgMz = value;                  break;
    case Var_D:                         D = value;  setLambdaEta(); break;
    case Var_mosaic:               mosaic = value;                  break;
    case Var_edisp:                 edisp = value;                  break;
    case Var_bFWHM:                 bFWHM = value; set_beamSigma(); break;
    case Var_sdistance:           setup.s = value; set_beamSigma(); break;
    case Var_bc2b:             setup.bc2b = value;                  break;
    case Var_wavelength: setup.wavelength = value; set_beamSigma(); break;
    case Var_pixelSize:          setup.pz = value; set_beamSigma(); break;
    case Var_qxzero:         setup.qrzero = value;                  break;
    case Var_nindex:         setup.nindex = value;                  break;
    case Var_T:                         T = value;  setLambdaEta(); break;
    case Var_Kt:                       Kt = value;                  break; //new
    case Var_at:                       at = value;                  break; //new
    case Var_Ls:                       Ls = value;                  break; //new
    case Var_divergeX:           divergeX = value;                  break; //new
    case Var_divergeZ:           divergeZ = value;                  break; //new
    case Var_teff:                   teff = value; 		            break; //new
    case Var_WrongInput: cerr << "Var_WrongInput at Para::setValue" << endl; break;
    default: cerr << "wrong index input to Para::setValue" << endl; break;
  }
}


double Para::getValue(int index)
{
  double value = 0;
  switch (index) {
    case Var_Kc:         value = Kc; break;
    case Var_B:          value = B; break;
    case Var_avgLr:      value = avgLr; break;
    case Var_avgMz:      value = avgMz; break;
    case Var_D:          value = D; break;
    case Var_mosaic:     value = mosaic; break;
    case Var_edisp:      value = edisp; break;
    case Var_bFWHM:      value = bFWHM; break;
    case Var_sdistance:  value = setup.s; break;
    case Var_bc2b:       value = setup.bc2b; break;
    case Var_wavelength: value = setup.wavelength; break;
    case Var_pixelSize:  value = setup.pz; break;
    case Var_qxzero:     value = setup.qrzero; break;
    case Var_nindex:     value = setup.nindex; break;
    case Var_T:          value = T; break;
    case Var_Kt:         value = Kt; break; //new
    case Var_at:         value = at; break; //new
    case Var_Ls:         value = Ls; break; //new
	case Var_divergeX:   value = divergeX; break; //new
	case Var_divergeZ:   value = divergeZ; break; //new
	case Var_teff:       value = teff; break; //new
    case Var_WrongInput: cerr << "Var_WrongInput at Para::getValue" << endl; break;
    default: cerr << "wrong index input to Para_getValue" << endl; break;
  }
  return value;
}


/****************************************************************************************
This function updates values of lambda and eta.
Should be called whenever a value of Kc, B, D, or T gets modified
****************************************************************************************/
void Para::setLambdaEta()
{
  lambda = 1e16 * sqrt(Kc/B) / D;
  eta = 2.17 * (T+273.15) / sqrt(Kc*B) / D / D;
}


void Para::set_beamSigma()
{
	double deltaQPerPixel = 2 * PI * setup.pz / setup.wavelength / setup.s;
	beamSigma = deltaQPerPixel * bFWHM;
}


/******************************************************************************
This function determines which parameter is free based on the input vector.
It then stores the corresponding index for each free parameter in idx vector.
Below is a list of free parameters and their corresponding indices. These
are defined through enum Var (see nfit.h).
0. Kc
1. B
2. avgLr
3. avgMz
4. D
5. mosaic
6. edisp
7. bFWHM
8. sdistance
90. bc2b
10. wavelength
11. pixelSize
12. qxzero
13. nindex
14. T
15. Kt
16. at
17. Ls
18. divergeX
19. divergez
20. teff
******************************************************************************/
void Para::set_index (const vector<bool>& x)
{
  for (int i = 0; i < g_totalNumOfFitParams; i++) {
    isFree[i] = x[i];
  }
  
  idx.clear();
  for (int i = 0; i < g_totalNumOfFitParams; i++) {
    if (isFree[i] == true) {
      idx.push_back(i);
    }
  }
  
  nfit = idx.size();
  
  set_xp();
  set_epsfcn();
}


void Para::set_xp ()
{
  xp.clear();
  typedef vector<int>::size_type sz;
  for (sz i = 0; i < idx.size(); i++) {
    xp.push_back(getValue(idx[i]));
  }
}


void Para::set_epsfcn()
{
  epsfcn.clear();
  typedef vector<int>::size_type sz;
  for (sz i = 0; i < idx.size(); i++) {
    if (idx[i] == Var_D) {
      epsfcn.push_back(0.0001);
    } else if (idx[i] == Var_bc2b) {
      epsfcn.push_back(0.001);  
    } else {
      epsfcn.push_back(0.01);
    }
  }
}


/******************************************************************************
This function converts an input string to a correspoinding enum value.
See the definition of enum Var in nfit.h

inString: an input string
******************************************************************************/
Var stringToVarEnum(const string& inString)
{
  if      (inString == "Kc")          return Var_Kc;
  else if (inString == "B")           return Var_B;
  else if (inString == "D")           return Var_D;
  else if (inString == "T")           return Var_T;
  else if (inString == "Lr")          return Var_avgLr;
  else if (inString == "Mz")          return Var_avgMz;
  else if (inString == "edisp")       return Var_edisp;
  else if (inString == "mosaic")      return Var_mosaic;
  else if (inString == "wavelength")  return Var_wavelength;
  else if (inString == "pixelSize")   return Var_pixelSize;
  else if (inString == "bFWHM")       return Var_bFWHM;
  else if (inString == "s")           return Var_sdistance;
  else if (inString == "bc2b")        return Var_bc2b;
  else if (inString == "qxzero")      return Var_qxzero;
  else if (inString == "nindex")      return Var_nindex;
  else if (inString == "Kt")          return Var_Kt; //new
  else if (inString == "at")          return Var_at;  //new
  else if (inString == "Ls")          return Var_Ls;  //new
  else if (inString == "divergeX")    return Var_divergeX;  //new
  else if (inString == "divergeZ")    return Var_divergeZ;  //new
  else if (inString == "teff")        return Var_teff;  //new
  else {
	  cerr << "Wrong input to stringToVarEnum" << endl;
    return Var_WrongInput;
  }
}
