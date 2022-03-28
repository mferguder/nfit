#ifndef GUARD_UTABLE_H
#define GUARD_UTABLE_H

#include <gsl/gsl_spline.h>
#include "interp2d_spline.h"
#include <vector>
#include <math.h>

struct f_params {
  	double rho;
  	double ell;
	int n;
};

double bessel_J0(double);
double bessel_J1(double);
double C_n_gr_one(double, double);
double expIntegral(double);

class Utable {
public:
    Utable();
    ~Utable(){cleanup();}
    void cleanup();
    void interp2d_splines_accs_cleanup();
    double heightDiff(double, double, int, double); 
    double dointeg(double, double, int, double);
    static double my_integrand(double, void *);
    double getHeightDiffFunc(int, double, double, double, double, double, double);
    double nzero_small_rho(double, double, double);
	double upperlim_correct(int, double, double, double);
    double AsymForm(double, double, int, double);
    double do_interp2d(double, double, int, double);
    void writeUtableFile();
	void readUtableFile(const std::string&);
	void readUtableFile(const char *);

//	void write_nstableFile(); //new
//	void read_nstableFile(const std::string&); //new
//	void read_nstableFile(const char *); //new
private:
    static const double pi;
    std::vector<gsl_interp_accel *> ellaccs;
    std::vector<gsl_interp_accel *> rhoaccs;
    std::vector<interp2d_spline *> splines;
    double minrho;   
    double maxrho;
    double minell;
    double maxell;
    double maxn;
    double upperlim_set;

//	std::vector<interp2d_spline *> splines_ns; //new
//	double numlim_values; //new
//	std::vector<double> limvec; //new: need for later interpolation

    friend class ModelCalculator;
};

#endif
