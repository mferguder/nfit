#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "utable.h"
#include <gsl/gsl_integration.h>
#include "interp2d_spline.h"

using std::cout; using std::endl; using std::cin; 
using std::string; using std::getline; using std::istringstream;
using std::vector; using std::ofstream; using std::ifstream;
using std::cerr;

#define WORKSPACE_SIZE 100000
#define KEY 6
const double Utable::pi = 3.1415926535897932384626433832795;

//Constructor
Utable::Utable()
{
  //Turn off the error handler in gsl_integration_qag(s)
  gsl_set_error_handler_off();
}

//Clean up resouces, should be called before the object is destructed
void Utable::cleanup()
{
    interp2d_splines_accs_cleanup();
}

/* Cleans up the three vectors, splines, rhoaccs, and ellaccs */
/* Also cleans up splines_ns */
void Utable::interp2d_splines_accs_cleanup()
{
	typedef vector<interp2d_spline *>::size_type sz;
	for (sz i = 0; i < splines.size(); i++) {
		if(splines[i] != NULL) interp2d_spline_free(splines[i]);
	}
	
	typedef vector<gsl_interp_accel *>::size_type sz2;
	for (sz2 i = 0; i < rhoaccs.size(); i++) {
		if(rhoaccs[i] != NULL) gsl_interp_accel_free(rhoaccs[i]);
	}

	typedef vector<gsl_interp_accel *>::size_type sz3;
	for (sz3 i = 0; i < ellaccs.size(); i++) {
		if(ellaccs[i] != NULL) gsl_interp_accel_free(ellaccs[i]);
	}

/*	typedef vector<interp2d_spline *>::size_type sz4;
	for (sz4 i = 0; i < splines_ns.size(); i++) {
		if(splines_ns[i] != NULL) interp2d_spline_free(splines_ns[i]);
	}*/
}

//Compute desired integral
double Utable::heightDiff(double rho, double ell, int n, double upperlim)
{
    return dointeg(rho, ell, n, upperlim);
}

// Computes the desired integral as a function of rho, ell, n, and chosen upper limit
double Utable::dointeg(double rho, double ell, int n, double upperlim)
{
  double memory = 100000;
  gsl_integration_workspace *work_ptr //work_ptr points to the workspace
    = gsl_integration_workspace_alloc (memory); //Should I allocate more memory for workspace?

  // KEY = 1 -> 15 point Gauss-Kronrod rule
  // KEY = 2 -> 21          "
  // KEY = 3 -> 31          "
  // KEY = 4 -> 41          "
  // KEY = 5 -> 51          "
  // KEY = 6 -> 61          "

  double lower_limit = 0;	/* lower limit */
  double upper_limit = upperlim;/* upper limit */
  double abs_error = 0;		/* to avoid round-off problems */
  double rel_error = 1.0e-4;	/* the result will usually be much better */
  double result;		/* the result from the integration */
  double error;			/* the estimated error from the integration */
	
  f_params pars;
  pars.rho = rho;	//integrand parameter
  pars.ell = ell;	//integrand parameter
  pars.n = n;		//integrand parameter

  gsl_function My_function;

  My_function.function = my_integrand;
  My_function.params = &pars;

  if (n != 0)
  {
  gsl_integration_qag (&My_function, lower_limit, upper_limit,
			abs_error, rel_error, memory, KEY, work_ptr, &result,
			&error);
  } else {
  
    if ( rho > 0.005) {
      gsl_integration_qags (&My_function, lower_limit, upper_limit,
			abs_error, rel_error, memory, work_ptr, &result,
			&error);
    } else {
      result = nzero_small_rho(rho, ell, upperlim);
    }
  }
  
  gsl_integration_workspace_free(work_ptr); /* Frees memory */ 

  return result;
}

// Defines integrand for later integration
double Utable::my_integrand (double x, void *params)
{
// The next line recovers the parameter values in f_params from the passed pointer
	f_params *pars = (f_params *) params;
	double part1, part2, part3, temp, ret;
	part1 = bessel_J0(sqrt(2*x)*pars->rho);
	temp = (x*x)/(1+x*pars->ell);
	part2 = pow(sqrt(1 + temp) - sqrt(temp),2*pars->n);
	part3 = sqrt(temp)*sqrt(1+temp);
	ret = (1-part1*part2)/part3;	
  return ret;
}

/******************************************************************************
zeroth order Bessel function of the first kind, 
taken from Numerical Recipes
******************************************************************************/
double bessel_J0(double x)
{
	double ax, z;
	double xx, y, ans, ans1, ans2;

	if ((ax = fabs(x)) < 8.0) {
		y = x * x;
		ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7 + y 
		* (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
		ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718 + y 
		* (59272.64853 + y * (267.8532712 + y * 1.0))));
		ans = ans1 / ans2;
	} else {
		z = 8.0 / ax;
		y = z * z;
		xx = ax - 0.785398164;
		ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 + y 
		* (-0.2073370639e-5 + y * 0.2093887211e-6)));
		ans2 = -0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5 +
		y * (0.7621095161e-6 - y * 0.934935152e-7)));
		ans = sqrt(0.636619772 / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
	}
	return ans;
}

/******************************************************************************
first order Bessel function of the first kind, 
taken from Numerical Recipes
******************************************************************************/
double bessel_J1(double x)
{
	double ax, z;
	double xx, y, ans, ans1, ans2;

	if ((ax = fabs(x)) < 8.0) {
		y = x * x;
		ans1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1
             + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
		ans2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
             + y * (99447.43394 + y * (376.9991397 + y * 1.0))));
		ans = ans1 / ans2;
	} else {
		z = 8.0 / ax;
		y = z * z;
		xx = ax - 2.356194491;
		ans1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
             + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
		ans2 = 0.04687499995 + y * (-0.2002690873e-3
             + y * (0.8449199096e-5 + y * (-0.88228987e-6
             + y * 0.105787412e-6)));
		ans = sqrt(0.636619772 / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
	}
	return ans;
}

/******************************************************************************
Returns the exponential integral. Use GSL special function.
******************************************************************************/
double expIntegral(double x)
{
	// for x = 0, exponential integral is not defined
	// for x > 700, gsl_sf_expint_E1 results in underflow error
	return (x > 0 && x < 700) ? gsl_sf_expint_E1(x) : 0;
}

/******************************************************************************
This function returns the height difference function evaluated at 
the input r, ell, and n, given the input xi, eta, and D.

rho = r / xi

There are three situations:
1) The input rho, ell, n are within the range of the utab.dat, so the 
   height difference function gets interpolated.

2) The input rho, ell, n are outside the range of the utab.dat and 
   the input rho is very close to zero. 
   In this case, the asymptotic form is not accurate
   and the interpolation is unavailable, so the full
   integration gets calculated on the fly. This should not 
   reduce the performance of NFIT because this case occurs
   very infrequently, namely only near r = 0.

   The cases for which n=0 are special. n=0 is not tabled at all
   and additionally, for small r, the numerical integration scheme
   returns incorrect values. Therefore a special function is used
   when n=0 and r is small.
   
3) The input rho, ell, n are outside the range of the utab.dat but
   the input rho is sufficietly large that the asymptotic form is 
   valid. The asymptotic form is used to do the calculation.
******************************************************************************/
double Utable::getHeightDiffFunc(int n, double r, double xi, 
                                 double ell, double eta, double D, double a)
{ 
	double ret;
	double rho = r / xi;
	double upperlim = 0.5 * pow(pi*xi/a,2);

//cout << "rho: " << rho << " ell: " << ell << endl;
//cout << "rho: " << rho << " ell: " << ell << " n:" << n << " ul:" << upperlim << 
//" eta:" << eta << endl;
	// If rho, ell, and n are within the tabled ranges use interpolation . . .
	// But also have to complete upper limit correction.
	if (rho >= minrho && rho <= maxrho &&
        ell >= minell && ell <= maxell &&
          n <= maxn) {
    
       ret = upperlim_correct(n, ell, rho, upperlim);
//cout << "from upperlim_correct --> " << ret << endl;
	// The input rho, ell, or n is outside the interpolation range
	// Resort to the asymptotic form . . .
	} else {
		// If rho is too small, the asymptotic form is not accurate.
		// --> Do the integration on the fly.
		if ( rho < minrho ) {
//cout << "full blown calculation" << endl;
			ret = dointeg(rho, ell, n, upperlim);		
		} else {
//cout << "Asymptotic form" << endl;
			ret = AsymForm(rho, ell, n, upperlim);
		}
	}
	return ret * eta * D * D / 2 / pi / pi;
}

// Calculates the height difference function 
// for n=0 and rho <= 0.005
double Utable::nzero_small_rho(double rho, double ell, double upperlim)
{
    double ret;
    double tau = upperlim;
    double beta = sqrt(1+(tau*tau)+tau*ell);
	ret = 2*tau + ell + 2*beta;
	ret = log(ret / (ell+2));
    ret = ret * (ell*ell - 2) / 2;
    ret = -ret + (ell * (beta - 1));
    ret = rho*rho * ret / 2;
    //cout << "rho: " << rho << " ell: " << ell << " ret: " << ret << endl;

    return ret;
}

// Calculate correction because of the choice of upper limit used to create the table.
// for n=0 and rho<20 calculate integral (done in do_interp2d)
// for n=1 and rho<10 calculate integral (done in do_interp2d)
double Utable::upperlim_correct(int n, double ell, double rho, double upperlim)
{
    double ret;

	ret = do_interp2d(rho, ell, n, upperlim);
	if ( n > 1 ) {
		double correction = C_n_gr_one(upperlim,ell) - C_n_gr_one(upperlim_set,ell);
    	ret += correction; 
	}

    return ret;
}

// Calculates upper limit correction for n > 1
// Defined in MJ thesis.
double C_n_gr_one(double x, double ell)
{
  double ret, f_c;
  
  f_c = sqrt( 1 + x*ell + x*x );
  ret = log(x) + ell*log( ell + 2*x + 2*f_c );
  ret -= log( 2 + x*ell + 2*f_c );

  return ret;
}

// Asymptotic form of the height-height correlation function.
// Defined in MJ thesis.
double Utable::AsymForm(double rho, double ell, int n, double upperlim)
{       
    //Euler's constant = 0.5772156649    
    double gamma = 0.5772156649;
    double ret = 0;
    double tau = upperlim;
    double beta;
    double temp;

    beta = sqrt(1 + ell*tau + tau*tau);
    temp = tau * pow(ell + 2*tau + 2*beta, ell);
    temp = log(temp / (2 + ell*tau + 2*beta));
    temp = temp - log(pow(ell+2, ell) / 2);
    
    if (n > 0) {
	double z = rho * rho / 4 / n;    
        ret += expIntegral(z) - ell * exp(-z) * (1 - z) / 4 / n;
    }
    ret += log(rho * rho) + (2 * gamma) + temp;
    
    return ret;
}

// Do interpolation of the 2-dimensional splines.
// If n==0 or n==1 compute height-height function on the fly
// otherwise use interpolants.
double Utable::do_interp2d(double rho, double ell, int n, double upperlim)
{       
	double ret;    
	double logrho = log(rho);
	if ( n==0 || n==1 ) {
		ret = dointeg(rho, ell, n, upperlim);
	} else {
		ret = interp2d_spline_eval(splines[n], logrho, ell, rhoaccs[n], ellaccs[n]);
	}

    return ret;
}

/******************************************************************************/
/* Create the utable. By default the utable file is named utab.dat */
void Utable::writeUtableFile()
{ 
    int numlayers = 2000;
    string input;
    int numrhoPoints = 400;
    double ellmax = 4;
    int numellPoints = 40;
    double upperlim = 5;
    
/*    cout << "Enter the number of layers (default is 100): " ;
    getline(cin, input);
    if (!input.empty()) {
        istringstream stream(input);
        stream >> numlayers;
        //utab.setN(num);
    }
    cout << "Enter the number of points in rho (default is 100): " ;
    getline(cin, input);
    if (!input.empty()) {
        istringstream stream(input);
        stream >> numrhoPoints;
    }
    cout << "Enter maximum ell value (Should be < 2; default is 2): " ;
    getline(cin, input);
    if (!input.empty()) {
        istringstream stream(input);
        stream >> ellmax;
    }
    cout << "Enter the number of points in ell (default is 100): " ;
    getline(cin, input);
    if (!input.empty()) {
        istringstream stream(input);
        stream >> numellPoints;
    }
    cout << "Enter the upper limit of the height-height correlation integral (default is 1000): " ;
    getline(cin, input);
    if (!input.empty()) {
        istringstream stream(input);
        stream >> upperlim;
    }*/

    // Specify rho min. and max.
    // For rho < rhomin compute integral at runtime
    double rhomin = pow(10,-4);
    double rhomax = 200;
    double logrhostep = log(rhomax/rhomin)/numrhoPoints;
    double logrhomin = log(rhomin);
    
    //Construct rhovec, which is a holder for log(rho) values used for the utable
    double rho;
    vector<double> rhovec;
    for (int i = 0; i<=numrhoPoints; i++) {
	rho = logrhomin + (i * logrhostep); 
    	rhovec.push_back(rho);
    }

    //Construct ellvec, which is a holder for all the ell values used for the utable
    double ell;
    double ellstep = ellmax / numellPoints;
    vector<double> ellvec;
    for (int i = 0; i <= numellPoints; i++) {
	ell = i * ellstep; 
    	ellvec.push_back(ell);
    }

    ofstream outfile("utab.dat");	
    //Write the first line, which contains parmeter values used to construct utable
    typedef vector<double>::size_type vec_sz;
    vec_sz sizerho = rhovec.size();
    vec_sz sizeell = ellvec.size();
    outfile << rhovec[0] << " " << rhovec[sizerho-1] << " " << logrhostep << " " << sizerho << " "
            << ellvec[0] << " " << ellvec[sizeell-1] << " " << ellstep << " " << sizeell << " " 
            << numlayers << " " << upperlim << endl;
    
    //Write the second line of utab.dat, which is a list of log(rho) values
    for (int i = 0; i < sizerho; i++) {
    	outfile << rhovec[i] << " ";
    }
    outfile << endl;

    //Write the third line of utab.dat, which is a list of ell values
    for (int i = 0; i < sizeell; i++) {
    	outfile << ellvec[i] << " ";
    }
    outfile << endl;

    /*for (int i = 0; i < sizerho; i++) {
    cout << exp(rhovec[i]) << endl;
    cout << utab.heightDiff(exp(rhovec[i]), ellvec[0], 0, 100) << endl;
    }*/
    
    //Create the utable by looping through all log(rho) and ell values for each n value
    //There will be sizeell number of rows per n value (sizerho number of columns)
    for (int n = 0; n <= numlayers; n++) {
    	for (int j = 0; j < sizeell ; j++) {
    		for (int i = 0; i < sizerho; i++) {
       			outfile << heightDiff(exp(rhovec[i]), ellvec[j], n, upperlim) << " ";
    		}
    		outfile << endl;
    	}
    }
    outfile.close();

    //Notify the user that the program has finished
    cout << "The output file should have " << 3 + (numlayers * sizeell) << " lines." << endl;

    cout << "A new table has been created using\n"
         << " N max = " << numlayers
         << "\nThe program will exit now" << endl;
    cout << "ell max = " << ellmax << endl;
}

//*********************************************************************//

/******************************************************************************
readUtableFile reads an ASCII formatted file for the table of
the integration defined inside the height difference function. 
It expects as the second (third) line,
the values of rho (ell) at which the integration was evaluated.

The following lines will be matrices from n=1 to n=(n max) of integration values.
The matrix rows will be the ell values and the columns the rho values. 

After reading the file, this function will build a vector
whose elements are interp2d_spline interpolants for each n as well as
a vector whose elements are interp2d_accel objects (one for rho values and one for ell values).

Note that the interpolation will be done as a function of log(rho).
******************************************************************************/
void Utable::readUtableFile(const string& filename)
{
	readUtableFile(filename.c_str());
}

void Utable::readUtableFile(const char *filename)
{    
    ifstream myfile;
    myfile.open(filename);
    if (!myfile) {
	  cerr << filename << ": file not found. The program exits." << endl;
		exit(1);
	} 

    string line;
    double d;
    // Skip first line. Need to get integral upper limit to 
    // be used in later calculations
	vector<double> utab_params;
    if (getline(myfile, line)) {
        istringstream iss(line);
        while (iss >> d) utab_params.push_back(d); 
    }
    // Get upper limit of integration used to create utable.
    upperlim_set = utab_params[9];
    cout << "utable upper integration limit = " << upperlim_set << endl;

    // Store second line = log(rho) values
    vector<double> rhovec;
    if (getline(myfile, line)) {
        istringstream iss(line);
        while (iss >> d) rhovec.push_back(d); 
    }
    // Store third line = ell values
    vector<double> ellvec;
    if (getline(myfile, line)) {
        istringstream iss(line);
        while (iss >> d) ellvec.push_back(d); 
    }

    // Set variables equal to the size of the rho and ell vectors.
    typedef vector<double>::size_type vec_sz;
    vec_sz sizerho = rhovec.size();
    vec_sz sizeell = ellvec.size();

    // Store minimmum and maximum values of rhovec and ellvec
    minrho = exp(rhovec[0]);
    maxrho = exp(rhovec[sizerho-1]);
    minell = ellvec[0];
    maxell = ellvec[sizeell-1];
//cout << "# of ells = " << sizeell << endl;
//cout << "minrho: " << minrho << " maxrho: " << maxrho << endl;
//cout << "minell: " << minell << " maxell: " << maxell << endl;

    // Builds splines, rhoaccs, and ellaccs vectors
    // Read all lines related to a particular n, build an interpolant and store
    vector<double> y;
    int row = 1;
    gsl_interp_accel *rhoacc, *ellacc;
    interp2d_spline* spline2d;
    while (getline(myfile, line)) {  
        while (row <= sizeell) {
            if (row != 1) getline(myfile, line);
            istringstream iss(line);    
	    while (iss >> d) {
		y.push_back(d);
        //cout << d << endl;
		//utab.push_back(d);
	    }
        row++;
        }
	row = 1;

    spline2d = interp2d_spline_alloc(interp2d_bicubic, sizerho, sizeell);
    interp2d_spline_init(spline2d, &rhovec[0], &ellvec[0], &y[0], sizerho, sizeell);
	
    splines.push_back(spline2d);	
    rhoacc = gsl_interp_accel_alloc();
	rhoaccs.push_back(rhoacc);
	ellacc = gsl_interp_accel_alloc();
    ellaccs.push_back(ellacc);

    y.clear();
    }
    // check the largest value of n
    maxn = splines.size();
    cout << "maximum tabled n: " << maxn << endl;
    cout << "ell max: " << maxell << endl;

    myfile.close();
}

