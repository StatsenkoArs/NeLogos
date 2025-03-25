#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include "NeLogos.h"

using namespace std;

complex<double> F(complex<double> z)
{
  return complex<double>(120) / (z * z * z * z * z * z);
}

double f(double t)
{
  return pow(t, 5);
}

double f_N(double t, double sigma, int N, double T)
{
  complex<double> result(0, 0);
  complex<double> s(sigma, 0);

  result += 0.5 * real(F(s));

  for (int k = 1; k <= N; ++k)
  {
    complex<double> arg = s + complex<double>(0, k * M_PI / T);
    result += real(F(arg)) * cos(k * M_PI * t / T);
    result -= imag(F(arg)) * sin(k * M_PI * t / T);
  }

  return exp(sigma * t) / T * real(result);
}

double Error(double fn, double f)
{
  return abs(f - fn);
}

bool ErrorIsOkay(double fn, double f)
{
  return Error(fn, f) < 0.5 * pow(10, -6);
}

int LookingForAnswer(double t, bool & retFlag)
{
  retFlag = true;
  for (double sigma = 0.3; sigma <= 10; sigma += 0.1)
  {
    for (double N = 2; N <= 50; ++N)
    {
      for (double T = t; T < t * 10; T += t / 10)
      {
        double fVal = f(t);
        double fNVal = f_N(t, sigma, N, T);

        if (ErrorIsOkay(fNVal, fVal))
        {
          cout << "t = " << t << endl;
          cout << "f_N(t, sigma, N, T) = " << fNVal << endl;
          cout << "f(t) = " << fVal << endl;
          cout << "Error = " << Error(fNVal, fVal) << endl;
          cout << "sigma = " << sigma << endl;
          cout << "N = " << N << endl;
          cout << "T = " << T << endl;
          return 0;
        }
      }
    }
  }
  retFlag = false;
  return {};
}

int main()
{
  vector<double> vals = { 0.1, 1, 5, 10, 50, 100 };
  for (int i = 0; i < vals.size(); ++i)
  {
    double t = vals[i];
    bool retFlag;
    int retVal = LookingForAnswer(t, retFlag);
    if (!retFlag) cout << "\nNo answer\n\n";
    else cout << endl;
  }

  return 0;
}