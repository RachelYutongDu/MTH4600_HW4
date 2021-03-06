
#include "Declarations.h"

int main() {

   double r, t, T, mu, sigma, dt, sqrtdt, S, S0, Z, Zbar, Z2bar,
          elapsed_time, t_star, stdhatZ, error, epsilon, n, Discount_factor,
          U, B, B_telda, S_telda, K, C, C_telda;
   int i, N, done, test;

   // Call option struck price
   K = 110;
    
   // Time to expiration.
   T = 0.5;

   // Number of stock price periods.
   N = 50;

   // Time increment per period.
   dt = T / N;

   // Compute this oft-used value once and for all.
   sqrtdt = sqrt(dt);

   // Risk-free interest rate.
   r = 0.05;

   // Compute the oft-used discount factor once and for all.
   Discount_factor = exp(-r*T);

   // Stock price volatility.
   sigma = .30;

   // Drift term.
   mu = r - sigma*sigma/2.0;

   // Initial stock price.
   S0 = 100.0;

   // Specify the 95% error tolerance.
   epsilon = 0.005;

   // Start the clock to time the computations.
   Time ();

   // Seed the RNG.
   MTUniform (1);

   // Print column headings for output to execution window.
   printf ("         n       Zbar       +/-         t       t*\n");

   // Initialize certain values.
   Zbar = Z2bar = done = n = test = 0;

   // Begin the simulation loop.
   while (!done) {

      // Initialize the stock price.
      S = S0;

      // Initialize the Brownian path.
      B = 0;

      // Initialize time.
      t = 0;

      // Simulate a stock price path.  Go forward period-by-period computing
      //   the next stock price.
      for (i = 1; i <= N; i++) {

         // Advance the path.
         U = MTUniform (0);
         Z = PsiInv (U); // Standard normal via inverst transform.
         B += sqrtdt * Z;

         // Advance time by one period.
         t += dt;

         // Compute the next stock price.
         S = S0 * exp (mu*t + sigma*B);
          
         //Compute the call option payoff
         C = fmax(S - K,0);

      }
      // S now is the value of the stock at time T for the simulated price path.
       
      // Compute the antithetic statistic
      B_telda = -B;
      S_telda = S0 * exp (mu*t + sigma*B_telda);
      C_telda = fmax(S_telda - K,0);
       
     //Compute the Z statistic and discount to time 0
      Z = (C + C_telda) * 0.5 * Discount_factor;

      // Update the simulation counter and the sample moments of C at time T.
      n ++;
      Zbar = ((n-1) * Zbar + Z) / n;
      Z2bar = ((n-1) * Z2bar + Z*Z) / n;
    

      // Update the error condition test counter.
      test ++;

      // Test the error condition periodically.
      if (test == 100000) {

         // Estimate the standard deviation of S.
         stdhatZ  = sqrt (Z2bar - Zbar*Zbar);

          // Estimate the error of the Zbar estimator.
         error = 1.96 * stdhatZ / sqrt(n);

         // Compute the elapsed time and estimate the time of completion.
         elapsed_time = Time ();
         t_star = elapsed_time * pow (error / epsilon, 2);

         // Report.
         printf ("%10.0f   %8.4f   %8.4f  %8.3f %8.3f\n", n, Zbar, error, elapsed_time, t_star);

         // Reset the "test" counter and see if error tolerance is met.
         test = 0;
         if (error < epsilon) {
            done = 1;
         }

      }

   }

   Exit ();

}

#include "Definitions.h"

