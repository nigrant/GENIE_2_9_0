//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 06, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <float.h>

#include <TMath.h>

#include "Messenger/Messenger.h"
#include "Utils/MathUtils.h"

//____________________________________________________________________________
double genie::utils::math::KahanSummation(double x[], unsigned int n)
{
// the Kahan summation algorithm - minimizes the error when adding a sequence
// of finite precision floating point numbers (compensated summation)

  double sum = x[0];
  double c   = 0.0;
  for(unsigned int i=1; i<n; i++) {
    double y = x[i]-c;
    double t = sum+y;
    c   = (t-sum) - y;
    sum = t;
  }
  return sum;
}
//____________________________________________________________________________
double genie::utils::math::KahanSummation(const vector<double> & x)
{
// the Kahan summation algorithm - minimizes the error when adding a sequence
// of finite precision floating point numbers (compensated summation)

  double sum = x[0];
  double c   = 0.0;
  for(unsigned int i=1; i<x.size(); i++) {
    double y = x[i]-c;
    double t = sum+y;
    c   = (t-sum) - y;
    sum = t;
  }
  return sum;
}
//____________________________________________________________________________
bool genie::utils::math::AreEqual(double x1, double x2)
{
  double err = 0.001*DBL_EPSILON;
  double dx  = TMath::Abs(x1-x2);
  if(dx<err) {
    LOG("Math", pINFO) << x1 << " := " << x2;
    return true;
  }
  return false;;
}
//____________________________________________________________________________
bool genie::utils::math::AreEqual(float x1, float x2)
{
  float err = FLT_EPSILON;
  float dx  = TMath::Abs(x1-x2);
  if(dx<err) {
    LOG("Math", pINFO) << x1 << " := " << x2;
    return true;
  }
  return false;;
}
//____________________________________________________________________________
bool genie::utils::math::IsWithinLimits(double x, Range1D_t range)
{
  return ( x >= range.min && x <= range.max );
}
//____________________________________________________________________________
bool genie::utils::math::IsWithinLimits(float x, Range1F_t range)
{
  return ( x >= range.min && x <= range.max );
}
//____________________________________________________________________________
bool genie::utils::math::IsWithinLimits(int i, Range1I_t range)
{
  return ( i >= range.min && i <= range.max );
}
//____________________________________________________________________________
double genie::utils::math::NonNegative(double x)
{
// this is used to handle very small numbers in sqrts

  return TMath::Max(0., x);
}
//____________________________________________________________________________
double genie::utils::math::NonNegative(float x)
{
// this is used to handle very small numbers in sqrts

  return TMath::Max( (float)0., x);
}
//____________________________________________________________________________
double genie::utils::math::Skellam(double mean1, double mean2, int k)
{
  // Calculates the probability of k from the Skellam distribution,
  // which is the difference between two Poisson distributions 
  // with means of mean1 and mean2. 

  double BesselI, result;

  if(mean1 >= 1e-6 && mean2 >= 1e-6)
    {
      //first calculate BesselI (this is modified Bessel function of the first kind)
      BesselI = 0.0;
      //sum converges after ~10 terms; use 30 terms to be on the safe side
      for(int n=0; n<30; n++)
        {
          if(n>=0 && n+k>=0)
            BesselI += (1.0 / (TMath::Factorial(n) * TMath::Factorial(n+k))) * pow(sqrt(mean1 * mean2), 2*n+k);
        }
      //exponent in pow term MUST be "0.5*k" 
      //if it is "k/2", it is rounded down to an int when k is an odd integer, e.g. 5/2 becomes 2 and not 2.5
      result = exp(-(mean1 + mean2)) * pow((mean1 / mean2), 0.5 * k) * BesselI;
    }
  else if(mean1 >= 1e-6 && mean2 < 1e-6 && mean2 >= 0)
    {
      //if mean2 is very small, k can only be >= 0; Skellam is effectively a Poisson with mean of mean1
      if(k>=0)
        result = exp(-mean1) * pow(mean1, k) / TMath::Factorial(k);
      else
        result = 0.0;
    }
  else if(mean1 < 1e-6 && mean1 >= 0 && mean2 > 1e-6)
    {
      //if mean1 is very small, k can only be <= 0; take k to be the negative of a Poisson with mean of mean2
      if(k<=0)
        result = exp(-mean2) * pow(mean2, -k) / TMath::Factorial(-k);
      else
        result = 0.0;
    }
  else
    {
      //LOG("Math", pINFO) "At least one of Mean1 and mean2 is zero or negative; the Skellam probability cannot be calculated.";
      std::cerr<<"At least one of Mean1 (= "<< mean1<<") and mean2 (= "<<mean2<<") is zero or negative; the Skellam probability cannot be calculated."<<std::endl;
      exit(1);
    }

  return result;
}
//____________________________________________________________________________
