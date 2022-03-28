#ifndef GUARD_NFIT_H
#define GUARD_NFIT_H

enum Var {
  Var_Kc, // 0 
  Var_B, // 1
  Var_avgLr,// 2
  Var_avgMz, // 3
  Var_D, // 4
  Var_mosaic, // 5
  Var_edisp, // 6
  Var_bFWHM, // 7
  Var_sdistance, // 8
  Var_bc2b, // 9
  Var_wavelength, // 10
  Var_pixelSize, // 11
  Var_qxzero, // 12
  Var_nindex, // 13
  Var_T, // 14
  Var_Kt, // 15
  Var_at, // 16
  Var_Ls, // 17 (new 2/16/15)
  Var_divergeX, // 18 (new 2/26/15)
  Var_divergeZ, // 19 (new 2/26/15)
  Var_teff, // 20 (new 3/30/15)
  Var_rst, // 21 mfe
  Var_WrongInput // 22 // was 21 mfe
};

extern const char *Varname[];

#endif

