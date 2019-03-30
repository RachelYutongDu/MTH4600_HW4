// This file declares various functions that are useful over the course
// of the semester.  It must be "included" above the main program.


// These are some standard C function libraries.
#include <stdlib.h>  // "standard library"
#include <math.h>    // various math functions, such as exp()
#include <stdio.h>   // various "input/output" functions
#include <time.h>    // functions for timing computations


// These functions are found in "Definitions.h".

// - Mersenne Twister RNG.
double MTUniform (unsigned int);

// - Histogram functions.
void Histogram (double, double, double, int, int);
void DiscreteHistogram (int, int, int, int);
void NormalHistogram (double, int, int);
void ExponentialHistogram (double, int, int);
void UniformHistogram (double, int, int);

// - Functions related to standard normals.
double Psi (double);
double PsiInv (double);

// - Miscellanceous functions.
void   Pause (void);
void   Exit (void);
double Time (void);
int    Equal (double, double, double);
int    GetInteger (char *);
double GetDouble (char *);

// - Linear algebra functions.
double **Array (int, int);
void Show (double **);
void Write (double **, FILE *);
double **Transpose (double **);
double **Multiply (double **, double **);
double **Invert (double **);
double **Copy (double **);
double **Identity (int);
double Det(double **);
double **Cholesky (double **);
double **ScalarMultiple (double, double **);
double **Add (double **, double **);
int Rows (double **);
int Columns (double **);
void Free (double **);
double **Evalues (double **);
double **Evector (double **, double **);
double **MeanZero (double **);
double **Covariance (double **);
double **Correlation (double **);
double **QRalgorithm (double **);

// - Black-Scholes call option functions.
double BlackScholes (double, double, double, double, double);
double ImpliedVol (double, double, double, double, double);

// - Yield Curve functions.
void   MakeCurve (double, double, double);
void   ComputeCurves (int);
void   AllocateTSRMemory (void);
double *GetTreasuryParCurve (char *);
double Decay (double);
double Hump (double);
double Twist (double);


#define ZERO 0.00000000001





