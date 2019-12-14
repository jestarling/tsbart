//  Pseudorandom numbers from a truncated Gaussian distribution.
//
//  Copyright (C) 2018 Andreas Dzemski, andreas.dzemski@gmail.com
//
//  This implements the zigguart algorithm from
//  N. Chopin, "Fast simulation of truncated Gaussian distributions",
//  Stat Comput (2011) 21:275-288
//
//  The code is based on the implementation by Guillaume Dollé, Vincent Mazet
//  available from http://miv.u-strasbg.fr/mazet/rtnorm/. In particular, the
//  rtchopin_twosided and rtchopin_onesided are derived from their code.
//  Andreas Dzemski adapted these functions to the R environment on Sep 13, 2018 by
//  changing the random number generation using the GNU Scientific library to
//  native R functions.
//
//  Licence: GNU General Public License Version 2
//  This program is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2 of the License, or (at your
//  option) any later version. This program is distributed in the hope that
//  it will be useful, but WITHOUT ANY WARRANTY; without even the implied
//  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details. You should have received a
//  copy of the GNU General Public License along with this program; if not,
//  see http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt


#include <Rcpp.h>
#include <cmath>
#include "rtn.hpp"


//------------------------------------------------------------
// Pseudorandom numbers from a truncated Gaussian distribution
// The Gaussian has parameters mu (default 0) and sigma (default 1)
// and is truncated on the interval [a,b].
// Returns the random variable x and its probability p(x).
//' Pseudorandom numbers from a Gaussian distribution that is truncated to an interval.
//'
//' @param lower lower bound of truncation interval
//' @param upper upper bound of truncation interval
//' @param mu mean of the normal random variable (before truncation)
//' @param sigma standard deviation of the normal random variable (before truncation)
//' @export
// [[Rcpp::export]]
double rtnorm(double lower, double upper, double mu, double sigma)
{
   double r;
   // Rcpp::NumericVector a = Rcpp::NumericVector(1);
   // Rcpp::NumericVector b = Rcpp::NumericVector(1);
   //
   // if (lower.size() != 1 || upper.size() != 1)
   // {
   //    Rcpp::stop("bounds must be scalars");
   // }

   double a = (lower - mu)/sigma;
   double b = (upper - mu)/sigma;

   // return(lower[0]);

   if (std::abs(a) > std::abs(b))
   {
      r = -rtnorm(-b, -a);
   }

   else if (a == R_NegInf)
   {
      r = Rcpp::rnorm(1, 0, 1)[0];
   }

   // If a in the right tail (a > xmax), use rejection algorithm with a truncated exponential proposal
   else if (a > Rtnorm::xmax) {
      if (b == R_PosInf)
         r = c_rtexp_onesided(a);
      else
         r = c_rtexp(a, b);
   }

   // If a in the left tail (a < xmin), use rejection algorithm with a Gaussian proposal
   else if (a < Rtnorm::xmin)
   {
      if (b == R_PosInf)
         r = c_rtnaive_onesided(a);
      else
         r = c_rtnaive(a, b);
   }

   // In other cases (xmin < a < xmax), use Chopin's algorithm
   else
   {
      if (b == R_PosInf)
         r = rtchopin_onesided(a);
      else
         r = rtchopin_twosided(a, b);
   }

   return sigma * r + mu;
}

//------------------------------------------------------------
// Compute y_l from y_k
double yl(int k)
{
   double yl0 = 0.053513975472;                  // y_l of the leftmost rectangle
   double ylN = 0.000914116389555;               // y_l of the rightmost rectangle

   if (k == 0)
      return yl0;

   else if(k == Rtnorm::N-1)
      return ylN;

   else if(k <= 1953)
      return Rtnorm::yu[k-1];

   else
      return Rtnorm::yu[k+1];
}

//------------------------------------------------------------
// Rejection algorithm with a truncated exponential proposal
double c_rtexp(double a, double b)
{
   int stop = false;
   double twoasq = 2*pow(a,2);
   double expab = exp(-a*(b-a)) - 1;
   double z, e;

   while(!stop)
   {
      z = log(1 + R::runif(0,1)*expab);
      e = -log(R::runif(0,1));
      stop = (twoasq*e > pow(z,2));
   }
   return a - z/a;
}

//------------------------------------------------------------
// One-sided rejection algorithm with a truncated exponential proposal
double c_rtexp_onesided(double a)
{
   bool stop = false;
   double lambda = 0.5 * (a + sqrt(pow(a, 2) + 4));
   double z;

   while(!stop)
   {
      z = a - log(R::runif(0,1))/lambda;
      stop = (-2* log(R::runif(0,1)) > pow(z - lambda, 2));
   }
   return z;
}

//------------------------------------------------------------
// Naive accept reject (two-sided)
double c_rtnaive(double a, double b)
{
   bool stop = false;
   double r;

   while(!stop)
   {
      r = R::rnorm(0,1);
      stop = (r>=a) && (r<=b);
   }

   return r;
}

//------------------------------------------------------------
// Naive accept reject  (one-sided)
double c_rtnaive_onesided(double a)
{
   bool stop = false;
   double r;

   while(!stop)
   {
      r = R::rnorm(0,1);
      stop = (r>=a);
   }

   return r;
}

//------------------------------------------------------------
// Rejection algorithm with a truncated Rayleigh proposal
double rtrayleigh(double a)
{
   bool stop = false;
   double asq = pow(a, 2)*0.5;
   double v, x;

   while(!stop)
   {
      v = R::runif(0,1);
      x = asq - log(R::runif(0,1));
      stop = (pow(v, 2) * x < a);
   }

   return sqrt (2*x);
}

//------------------------------------------------------------
// Chopin's one-sided ziggurat algorithm
double rtchopin_onesided(double a)
{
   bool stop = false;
   double r, sim, simy, u, d, ylk;
   int k, ka, i;

   if (a < Rtnorm::xmin)  //if a is too small, use simple reject algorithm
   {
      r = c_rtnaive_onesided(a);
      stop = true;
   }

   else if (a > Rtnorm::xmax) // if a is too large, use Devroye's algorithm
   {
      r = c_rtexp_onesided(a);
      stop = true;
   }

   else
   {
      i = Rtnorm::I0+floor(a*Rtnorm::INVH);
      ka = Rtnorm::ncell[i];
      while (!stop) {
         // sample integer between ka and N
         k = floor(R::runif(0,1) * (Rtnorm::N-ka+1)) + ka;
         if (k == Rtnorm::N)  //right tail (last box on the right)
         {
            r = c_rtexp_onesided(Rtnorm::x[Rtnorm::N]);
            stop = true;
         }
         else if (k <= ka+1)
         {
            sim = Rtnorm::x[k] + (Rtnorm::x[k+1] - Rtnorm::x[k]) * R::runif(0, 1);

            if (sim >= a)
            {
               simy = Rtnorm::yu[k]*R::runif(0, 1);
               if ( (simy < yl(k)) || (sim * sim + 2.*log(simy) + Rtnorm::ALPHA) < 0 )
               {
                  r = sim;
                  stop = true;
               }
            }
         }
         else
         {
            u = R::runif(0,1);
            simy = Rtnorm::yu[k] * u;
            d = Rtnorm::x[k+1] - Rtnorm::x[k];
            ylk = yl(k);
            if (simy < ylk)  // That's what happens most of the time
            {
               // sprintf(char_arr, "%.2f", ka);
               // Rcpp::stop(char_arr);
               r = Rtnorm::x[k] + u*d*Rtnorm::yu[k]/ylk;
               stop = true;
            }
            else
            {
               sim = Rtnorm::x[k] + d * R::runif(0,1);

               // Otherwise, check you're below the pdf curve
               if ((sim * sim + 2*log(simy) + Rtnorm::ALPHA) < 0)
               {
                  r = sim;
                  stop = true;
               }
            }
         }
      } // end while loop

   }

   return r;
}

//------------------------------------------------------------
// Chopin's two-sided ziggurat algorithm
double rtchopin_twosided(double a, double b)
{
   const int kmin = 5;                                 // if kb-ka < kmin then use a rejection algorithm
   int stop = false;
   double r, z, e, ylk, simy, lbound, u, d, sim;
   int i, ka, kb, k;

   // Compute ka
   i = Rtnorm::I0 + floor(a*Rtnorm::INVH);
   ka = Rtnorm::ncell[i];

   // Compute kb
   if (b>=Rtnorm::xmax)
      kb = Rtnorm::N;
   else
   {
      i = Rtnorm::I0 + floor(b*Rtnorm::INVH);
      kb = Rtnorm::ncell[i];
   }

   // If |b-a| is small, use rejection algorithm with a truncated exponential proposal
   if (abs(kb-ka) < kmin)
   {
      r = c_rtexp(a, b);
      stop = true;
   }

   while (!stop)
   {
      // Sample integer between ka and kb
      k = floor(R::runif(0,1) * (kb-ka+1)) + ka;

      if(k == Rtnorm::N)
      {
         // Right tail
         lbound = Rtnorm::x[Rtnorm::N];
         z = -log(R::runif(0,1));
         e = -log(R::runif(0,1));
         z = z / lbound;

         if ((pow(z,2) <= 2*e) && (z < b-lbound))
         {
            // Accept this proposition, otherwise reject
            r = lbound + z;
            stop = true;
         }
      }

      else if ((k <= ka+1) || (k >= kb-1 && b < Rtnorm::xmax))
      {

         // Two leftmost and rightmost regions
         sim = Rtnorm::x[k] + (Rtnorm::x[k+1]-Rtnorm::x[k]) * R::runif(0,1);

         if ((sim >= a) && (sim <= b))
         {
            // Accept this proposition, otherwise reject
            simy = Rtnorm::yu[k]*R::runif(0,1);
            if ( (simy<yl(k)) || (sim * sim + 2*log(simy) + Rtnorm::ALPHA) < 0 )
            {
               r = sim;
               stop = true;
            }
         }
      }

      else // All the other boxes
      {
         u = R::runif(0,1);
         simy = Rtnorm::yu[k] * u;
         d = Rtnorm::x[k+1] - Rtnorm::x[k];
         ylk = yl(k);
         if(simy < ylk)  // That's what happens most of the time
         {
            r = Rtnorm::x[k] + u*d*Rtnorm::yu[k]/ylk;
            stop = true;
         }
         else
         {
            sim = Rtnorm::x[k] + d * R::runif(0,1);

            // Otherwise, check you're below the pdf curve
            if((sim * sim + 2*log(simy) + Rtnorm::ALPHA) < 0)
            {
               r = sim;
               stop = true;
            }
         }

      }
   }

   return r;
}
