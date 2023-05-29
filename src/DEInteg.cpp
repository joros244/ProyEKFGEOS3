#include "../include/SAT_Const.h"
#include "../include/sign_.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <string>

using namespace std;

double *DEInteg(function<double *(double, double *)> func, double t,
                double tout, double relerr, double abserr, int n_eqn,
                double *y) {
  // maxnum = 500;
  double twou = 2 * eps;
  double fouru = 4 * eps;

  map<string, int> DE_STATE;
  DE_STATE["DE_INIT"] = 1;     // Restart integration
  DE_STATE["DE_DONE"] = 2;     // Successful step
  DE_STATE["DE_BADACC"] = 3;   // Accuracy requirement could not be achieved
  DE_STATE["DE_NUMSTEPS"] = 4; // Permitted number of steps exceeded
  DE_STATE["DE_STIFF"] = 5;    // Stiff problem suspected
  DE_STATE["DE_INVPARAM"] = 6; // Invalid input parameters

  int State_ = DE_STATE["DE_INIT"];
  bool PermitTOUT = true; // Allow integration past tout by default
  int told = 0;

  // Powers of two (two(n)=2^n)
  double two[14] = {1.0,   2.0,   4.0,   8.0,    16.0,   32.0,   64.0,
                    128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0};

  double gstr[14] = {1.0,     0.5,     0.0833,  0.0417,  0.0264,
                     0.0188,  0.0143,  0.0114,  0.00936, 0.00789,
                     0.00679, 0.00592, 0.00524, 0.00468};

  double *yy = new double[n_eqn];
  double *wt = new double[n_eqn];
  double *p = new double[n_eqn];
  double *yp = new double[n_eqn];
  double **phi = new double *[n_eqn];
  for (int i = 0; i < n_eqn; i++) {
    phi[i] = new double[17];
  }
  // memset(yy, 0.0, sizeof(double) * n_eqn);
  // memset(wt, 0.0, sizeof(double) * n_eqn);
  // memset(p, 0.0, sizeof(double) * n_eqn);
  // memset(yp, 0.0, sizeof(double) * n_eqn);
  // memset(phi, 0.0, sizeof(double) * n_eqn * 17);
  double *g = new double[14];
  double *sig = new double[14];
  double *rho = new double[14];
  double *w = new double[13];
  double *alpha = new double[13];
  double *beta = new double[13];
  double *v = new double[13];
  double *psi_ = new double[13];
  // memset(g, 0.0, sizeof(double) * 14);
  // memset(sig, 0.0, sizeof(double) * 14);
  // memset(rho, 0.0, sizeof(double) * 14);
  // memset(w, 0.0, sizeof(double) * 13);
  // memset(alpha, 0.0, sizeof(double) * 13);
  // memset(beta, 0.0, sizeof(double) * 13);
  // memset(v, 0.0, sizeof(double) * 13);
  // memset(psi_, 0.0, sizeof(double) * 13);
  //  while(true)

  // Return, if output time equals input time

  if (t == tout) { // No integration
    exit(EXIT_FAILURE);
  }

  // Test for improper parameters

  double epsilon = max(relerr, abserr);

  if ((relerr < 0.0) ||
      // Negative relative error bound
      (abserr < 0.0) ||
      // Negative absolute error bound
      (epsilon <= 0.0) ||
      // Both error bounds are non-positive
      (State_ > DE_STATE["DE_INVPARAM"]) ||
      // Invalid status flag
      ((State_ != DE_STATE["DE_INIT"]) && (t != told))) {
    State_ = DE_STATE["DE_INVPARAM"]; // Set error code
    exit(EXIT_FAILURE);               // Exit
  }

  // On each call set interval of integration and counter for
  // number of steps. Adjust input error tolerances to define
  // weight vector for subroutine STEP.

  double del = tout - t;
  double absdel = fabs(del);

  double tend = t + 100.0 * del;
  if (!PermitTOUT) {
    tend = tout;
  }

  int nostep = 0;
  int kle4 = 0;
  bool stiff = false;
  double releps = relerr / epsilon;
  double abseps = abserr / epsilon;

  bool start;
  double x, delsgn, h;
  if ((State_ == DE_STATE["DE_INIT"]) || (delsgn * del <= 0.0)) {
    // On start and restart also set the work variables x and yy(*),
    // store the direction of integration and initialize the step size
    start = true;
    x = t;
    yy = y;
    delsgn = sign_(1.0, del);
    h = sign_(max(fouru * abs(x), abs(tout - x)), tout - x);
  }

  int ki = 0, kold = 0;
  bool OldPermit;
  double *yout = new double[n_eqn];
  double *ypout = new double[n_eqn];
  while (true) { // Start step loop

    // If already past output point, interpolate solution and return
    if (fabs(x - t) >= absdel) {
      g[1] = 1.0;
      rho[1] = 1.0;
      double hi = tout - x;
      ki = kold + 1;

      // Initialize w[*] for computing g[*]
      double temp1;
      for (int i = 0; i < ki; i++) {
        temp1 = i;
        w[i + 1] = 1.0 / temp1;
      }
      // Compute g[*]
      double term = 0.0;
      for (int j = 1; j < ki; j++) {
        double psijm1 = psi_[j];
        double gamma = (hi + term) / psijm1;
        double eta = hi / psijm1;
        for (int i = 0; i < ki + 1 - j; i++) {
          w[i + 1] = gamma * w[i + 1] - eta * w[i + 2];
        }
        g[j + 1] = w[1];
        rho[j + 1] = gamma * rho[j];
        term = psijm1;
      }

      // Interpolate for the solution yout and for
      // the derivative of the solution ypout
      int i;
      for (int j = 0; j < ki; j++) {
        i = ki + 1 - j;
        for (int q = 0; q < n_eqn; q++) {
          yout[q] += g[i + 1] * phi[q][i + 1];
          ypout[q] += rho[i + 1] * phi[q][i + 1];
        }
      }
      for (int q = 0; q < n_eqn; q++) {
        yout[q] = y[q] + hi * yout[q];
      }
      y = yout;
      State_ = DE_STATE["DE_DONE"]; // Set return code
      t = tout;                     // Set independent variable
      told = t;                     // Store independent variable
      OldPermit = PermitTOUT;
      return y; // Normal exit
    }

    // If cannot go past output point and sufficiently close,
    // extrapolate and return
    if (!PermitTOUT && (fabs(tout - x) < fouru * fabs(x))) {
      h = tout - x;
      yp = func(x, yy); // Compute derivative yp(x)
      for (int q = 0; q < n_eqn; q++) {
        y[q] = yy[q] + h * yp[q]; // Extrapolate vector from x to tout
      }
      State_ = DE_STATE["DE_DONE"]; // Set return code
      t = tout;                     // Set independent variable
      told = t;                     // Store independent variable
      OldPermit = PermitTOUT;
      return y; // Normal exit
    }

    // Test for too much work
    //   if (nostep >= maxnum)
    //       State_ = DE_STATE.DE_NUMSTEPS; // Too many steps
    //       if (stiff)
    //           State_ = DE_STATE.DE_STIFF;// Stiffness suspected
    //       end
    //       y         = yy;                // Copy last step
    //       t         = x;
    //       told      = t;
    //       OldPermit = true;
    //       return;                        // Weak failure exit
    //   end

    // Limit step size, set weight vector and take a step
    h = sign_(min(abs(h), abs(tend - x)), h);
    for (int l = 0; l < n_eqn; l++) {
      wt[l] = releps * fabs(yy[l]) + abseps;
    }

    //   Step
    //
    // Begin block 0
    //
    // Check if step size or error tolerance is too small for machine
    // precision.  If first step, initialize phi array and estimate a
    // starting step size. If step size is too small, determine an
    // acceptable one.
    //

    bool crash;
    if (fabs(h) < fouru * fabs(x)) {
      h = sign_(fouru * abs(x), h);
      crash = true;
      exit(EXIT_FAILURE); // Exit
    }

    double p5eps = 0.5 * epsilon;
    crash = false;
    g[1] = 1.0;
    g[2] = 0.5;
    sig[1] = 1.0;

    int ifail = 0;

    // If error tolerance is too small, increase it to an
    // acceptable value.

    double round = 0.0;
    for (int l = 0; l < n_eqn; l++) {
      round = round + (y[l] * y[l]) / (wt[l] * wt[l]);
    }
    round = twou * sqrt(round);
    if (p5eps < round) {
      epsilon = 2.0 * round * (1.0 + fouru);
      crash = true;
      exit(EXIT_FAILURE); // Exit
    }

    bool phase1, nornd;
    double hnew, absh, hold;
    int k;
    if (start) {
      // Initialize. Compute appropriate step size for first step.
      yp = func(x, y);
      double sum = 0.0;
      for (int l = 0; l < n_eqn; l++) {
        phi[l][1] = yp[l];
        phi[l][2] = 0.0;
        sum = sum + (yp[l] * yp[l]) / (wt[l] * wt[l]);
      }
      sum = sqrt(sum);
      absh = fabs(h);
      if (epsilon < 16.0 * sum * h * h) {
        absh = 0.25 * sqrt(epsilon / sum);
      }
      h = sign_(max(absh, fouru * abs(x)), h);
      hold = 0.0;
      hnew = 0.0;
      k = 1;
      kold = 0;
      start = false;
      phase1 = true;
      nornd = true;
      if (p5eps <= 100.0 * round) {
        nornd = false;
        for (int l = 0; l < n_eqn; l++) {
          phi[l][15] = 0.0;
        }
      }
    }

    //
    // End block 0
    //

    //
    // Repeat blocks 1, 2 (and 3) until step is successful
    //
    k = 0;
    int kp1, kp2, km1, km2, ns;
    double erk, erkm1;
    int knew;
    while (true) {

      //
      // Begin block 1
      //
      // Compute coefficients of formulas for this step. Avoid computing
      // those quantities not changed when step size is not changed.
      //

      kp1 = k + 1;
      kp2 = k + 2;
      km1 = k - 1;
      km2 = k - 2;

      // ns is the number of steps taken with size h, including the
      // current one. When k<ns, no coefficients change.

      ns = 0;
      double hold = 0.0;
      if (fabs(h - hold) > pow(10, -12)) {
        ns = 0;
      }
      if (ns <= kold) {
        ns = ns + 1;
      }
      int nsp1 = ns + 1;

      if (k >= ns) {
        // Compute those components of alpha[*], beta[*], psi[*], sig[*]
        // which are changed
        beta[ns] = 1.0;
        int realns = ns;
        alpha[ns] = 1.0 / realns;
        double temp1 = h * realns;
        sig[nsp1] = 1.0;
        if (k >= nsp1) {
          int im1, reali;
          double temp2;
          for (int i = nsp1; i < k; i++) {
            im1 = i - 1;
            temp2 = psi_[im1 + 1];
            psi_[im1 + 1] = temp1;
            beta[i + 1] = beta[im1 + 1] * psi_[im1 + 1] / temp2;
            temp1 = temp2 + h;
            alpha[i + 1] = h / temp1;
            reali = i;
            sig[i + 2] = reali * alpha[i + 1] * sig[i + 1];
          }
        }
        psi_[k + 1] = temp1;

        // Compute coefficients g[*]; initialize v[*] and set w[*].
        if (ns > 1) {
          // If order was raised, update diagonal part of v[*]
          if (k > kold) {
            double temp4 = k * kp1;
            v[k + 1] = 1.0 / temp4;
            int nsm2 = ns - 2;
            int i;
            for (int j = 0; j < nsm2; j++) {
              i = k - j;
              v[i + 1] = v[i + 1] - alpha[j + 2] * v[i + 2];
            }
          }

          // Update V[*] and set W[*]
          int limit1 = kp1 - ns;
          double temp5 = alpha[ns + 1];
          for (int iq = 0; iq < limit1; iq++) {
            v[iq + 1] = v[iq + 1] - temp5 * v[iq + 2];
            w[iq + 1] = v[iq + 1];
          }
          g[nsp1 + 1] = w[1];
        } else {
          for (int iq = 0; iq < k; iq++) {
            double temp3 = iq * (iq + 1);
            v[iq + 1] = 1.0 / temp3;
            w[iq + 1] = v[iq + 1];
          }
        }

        // Compute the g[*] in the work vector w[*]
        int nsp2 = ns + 2;
        if (kp1 >= nsp2) {
          for (int i = nsp2 - 1; i < kp1; i++) {
            int limit2 = kp2 - i;
            double temp6 = alpha[i];
            for (int iq = 0; iq < limit2; iq++) {
              w[iq + 1] = w[iq + 1] - temp6 * w[iq + 2];
            }
            g[i + 1] = w[2];
          }
        }
      } // if K>=NS

      //
      // End block 1
      //

      //
      // Begin block 2
      //
      // Predict a solution p[*],evaluate derivatives using predicted
      // solution, estimate local error at order k and
      // errors at orders k, k-1, k-2 as if constant step size
      // were used.

      // Change phi to phi star
      double temp1;
      if (k >= nsp1) {
        for (int i = nsp1; i < k; i++) {
          temp1 = beta[i];
          for (int l = 0; l < n_eqn; l++) {
            phi[l][i] = temp1 * phi[l][i];
          }
        }
      }

      // Predict solution and differences
      for (int l = 0; l < n_eqn; l++) {
        phi[l][kp2] = phi[l][kp1];
        phi[l][kp1] = 0.0;
        p[l] = 0.0;
      }
      int i, ip1;
      double temp2;
      for (int j = 0; j < k; j++) {
        i = kp1 - j;
        ip1 = i;
        temp2 = g[i];
        for (int l = 0; l < n_eqn; l++) {
          p[l] = p[l] + temp2 * phi[l][i];
          phi[l][i] = phi[l][i] + phi[l][ip1];
        }
      }
      if (nornd) {
        for (int i = 0; i < n_eqn; i++) {
          p[i] = y[i] + h * p[i];
        }
      } else {
        double tau;
        for (int l = 0; l < n_eqn; l++) {
          tau = h * p[l] - phi[l][15];
          p[l] = y[l] + tau;
          phi[l][16] = (p[l] - y[l]) - tau;
        }
      }
      double xold = x;
      x = x + h;
      yp = func(x, p);
      absh = fabs(h);

      // Estimate errors at orders k, k - 1, k - 2
      double erkm2 = 0.0;
      erkm1 = 0.0;
      erk = 0.0;

      double temp3, temp4;
      for (int l = 0; l < n_eqn; l++) {
        temp3 = 1.0 / wt[l];
        temp4 = yp[l] - phi[l][1];
        if (km2 > 0) {
          erkm2 = erkm2 + ((phi[l][km1] + temp4) * temp3) *
                              ((phi[l][km1] + temp4) * temp3);
        }
        if (km2 >= 0) {
          erkm1 = erkm1 +
                  ((phi[l][k] + temp4) * temp3) * ((phi[l][k] + temp4) * temp3);
        }
        erk = erk + (temp4 * temp3) * (temp4 * temp3);
      }

      if (km2 > 0) {
        erkm2 = absh * sig[km1] * gstr[km2] * sqrt(erkm2);
      }
      if (km2 >= 0) {
        erkm1 = absh * sig[k] * gstr[km1] * sqrt(erkm1);
      }

      double temp5 = absh * sqrt(erk);
      double err = temp5 * (g[k] - g[kp1]);
      erk = temp5 * sig[kp1] * gstr[k];
      knew = k;

      // Test if order should be lowered
      if (km2 > 0) {
        if (max(erkm1, erkm2) <= erk) {
          knew = km1;
        }
      }
      if (km2 == 0) {
        if (erkm1 <= 0.5 * erk) {
          knew = km1;
        }
      }

      //
      // End block 2
      //

      //
      // If step is successful continue with block 4, otherwise repeat
      // blocks 1 and 2 after executing block 3

      bool success = (err <= epsilon);

      if (!success) {

        //
        // Begin block 3
        //

        // The step is unsuccessful.Restore x,phi[ *, *], psi[*].If
        // 3rd consecutive failure, set order to 1. If step fails more
        // than 3 times, consider an optimal step size.Double error
        // tolerance and return if estimated step size is too small
        // for  machine precision.
        //

        // Restore x, phi[ *, *] and psi[*]
        phase1 = false;
        x = xold;
        for (int i = 0; i < k; i++) {
          temp1 = 1.0 / beta[i + 1];
          ip1 = i + 1;
          for (int l = 0; l < n_eqn; l++) {
            phi[l][i + 1] = temp1 * (phi[l][i + 1] - phi[l][ip1 + 1]);
          }
        }

        if (k >= 2) {
          for (int i = 1; i < k; i++) {
            psi_[i] = psi_[i + 1] - h;
          }
        }

        // On third failure, set order to one.
        // Thereafter, use optimal step size
        ifail = ifail + 1;
        temp2 = 0.5;
        if (ifail > 3) {
          if (p5eps < 0.25 * erk) {
            temp2 = sqrt(p5eps / erk);
          }
        }
        if (ifail >= 3) {
          knew = 1;
        }
        h = temp2 * h;
        k = knew;
        if (fabs(h) < fouru * fabs(x)) {
          crash = true;
          h = sign_(fouru * fabs(x), h);
          epsilon = epsilon * 2.0;
          exit(EXIT_FAILURE);
        }
        // Exit

        //
        // End block 3, return to start of block 1
        //

      } // end if(success)

      if (success) {
        break;
      }
    }

    //
    // Begin block 4
    //
    // The step is successful.Correct the predicted solution,
    // evaluate the derivatives using the corrected solution and
    // update the differences.Determine best order and step size
    // for next step.
    //

    kold = k;
    hold = h;

    // Correct and evaluate
    double temp1 = h * g[kp1];
    if (nornd) {
      for (int l = 0; l < n_eqn; l++) {
        y[l] = p[l] + temp1 * (yp[l] - phi[l][1]);
      }
    } else {
      double aux;
      for (int l = 0; l < n_eqn; l++) {
        aux = temp1 * (yp[l] - phi[l][1]) - phi[l][16];
        y[l] = p[l] + aux;
        phi[l][15] = (y[l] - p[l]) - aux;
      }
    }
    yp = func(x, y);

    // Update differences for next step
    for (int l = 0; l < n_eqn; l++) {
      phi[l][kp1] = yp[l] - phi[l][1];
      phi[l][kp2] = phi[l][kp1] - phi[l][kp2];
    }
    for (int i = 0; i < k; i++) {
      for (int l = 0; l < n_eqn; l++) {
        phi[l][i + 1] = phi[l][i + 1] + phi[l][kp1];
      }
    }

    // Estimate error at order k + 1 unless
    // - in first phase when always raise order,
    // - already decided to lower order,
    // - step size not constant so estimate unreliable
    double erkp1 = 0.0;
    if ((knew == km1) || (k == 12)) {
      phase1 = false;
    }

    if (phase1) {
      k = kp1;
      erk = erkp1;
    } else {
      if (knew == km1) { // lower order
        k = km1;
        erk = erkm1;
      } else {
        if (kp1 <= ns) {
          for (int l = 0; l < n_eqn; l++) {
            erkp1 = erkp1 + (phi[l][kp2] / wt[l]) * (phi[l][kp2] / wt[l]);
          }
          erkp1 = absh * gstr[kp1] * sqrt(erkp1);
          // Using estimated error at order k+1, determine
          // appropriate order for next step
          if (k > 1) {
            if (erkm1 <= min(erk, erkp1)) {
              // lower order
              k = km1;
              erk = erkm1;
            } else {
              if ((erkp1 < erk) && (k != 12)) {
                // raise order
                k = kp1;
                erk = erkp1;
              }
            }
          } else if (erkp1 < 0.5 * erk) {
            // raise order
            // Here erkp1 < erk < max(erkm1,ermk2) else
            // order would have been lowered in block 2.
            // Thus order is to be raised
            k = kp1;
            erk = erkp1;
          }
        } // end if kp1<=ns
      }   // end if knew!=km1
    }     // end if !phase1

    // With new order determine appropriate step size for next step
    double temp2;
    double r;
    if (phase1 || (p5eps >= erk * two[k + 1])) {
      hnew = 2.0 * h;
    } else {
      if (p5eps < erk) {
        temp2 = k + 1;
        r = p5eps / pow(erk, (1.0 / temp2));
        hnew = absh * max(0.5, min(0.9, r));
        hnew = sign_(max(hnew, fouru * fabs(x)), h);
      } else {
        hnew = h;
      }
    }
    h = hnew;

    //
    // End block 4
    //

    // Test for too small tolerances
    if (crash) {
      State_ = DE_STATE["DE_BADACC"];
      relerr = epsilon * releps; // Modify relative and absolute
      abserr = epsilon * abseps; // accuracy requirements
      y = yy;                    // Copy last step
      t = x;
      told = t;
      OldPermit = true;
      exit(EXIT_FAILURE); // Weak failure exit
    }

    nostep = nostep + 1; // Count total number of steps

    // Count number of consecutive steps taken with the order of
    // the method being less or equal to four and test for stiffness
    kle4 = kle4 + 1;
    if (kold > 4) {
      kle4 = 0;
    }
    if (kle4 >= 50) {
      stiff = true;
    }
  } // End step loop

  //   if ( State_==DE_STATE.DE_INVPARAM )
  //       error ('invalid parameters in DEInteg');
  //       exit;
  //   end
  //   if ( State_==DE_STATE.DE_BADACC )
  //       warning ('on','Accuracy requirement not achieved in DEInteg');
  //   end
  //   if ( State_==DE_STATE.DE_STIFF )
  //       warning ('on','Stiff problem suspected in DEInteg');
  //   end
  //   if ( State_ >= DE_STATE.DE_DONE )
  //       break;
  //   end
  //
  // end

  for (int i = 0; i < n_eqn; i++) {
    delete[] phi[i];
  }
  delete[] yy;
  delete[] wt;
  delete[] p;
  delete[] yp;
  delete[] phi;
  delete[] g;
  delete[] sig;
  delete[] rho;
  delete[] w;
  delete[] alpha;
  delete[] beta;
  delete[] v;
  delete[] psi_;
  delete[] yout;
  delete[] ypout;
}
