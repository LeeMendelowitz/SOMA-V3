/* Calculate probability of an R-MAP alignment to the reference restriction map
    Author: Lee Mendelowitz
    Date:   3/17/2012
   
   Assumptions: Adapted from Teague et al., High-resolution human genome structure by single-molecule
                analysis. PNAS 2010 (10.1073/pnas.0914638107)
                and Valouev et al., Alignment of optical maps. (10.1089/cmb.2006.13.442)
     - Missing Cuts: Restriction cuts in the R-MAP are a Bernoulli process with p=0.8.
     - False cuts: False cuts in the R-MAP are a Poisson process with arrival rate
       0.005 1/kb (or 1/200 kb)
     - Desorption: The probability of desorption (loss of small fragments in the R-MAP) of a fragment of
       size x is modeled by an exponential: f(x) = exp(-lambda*x). Let f(1.35 kbp) = 0.5, then
       lambda = log(2)/1.35. In other words, the probability of desorption of a fragment of
       1.35 kb is 0.5.
     - Sizing Error: Let y be the true fragment size in kbp, and X be a random variable
       representing the measured fragment size in kbp in the R-MAP. Then X ~ N(y,sigma^2 y) where
       sigma^2 = 0.300 (can be experiment dependent)
*/

#include <assert.h>
#include <math.h>
#include <iostream>
#include <boost/math/special_functions/factorials.hpp>



// ****************************************************************************
namespace probConstants
{
    const double CUT_PROB = 0.8; // Prob of a true restriction site appearing in R-MAP
    const double FALSE_CUT_RATE = 5E-3*1E-3; // Arrival rate of false cuts (per bp).
    const double MEDIAN_DESORB_FRAG = 1350; // Median desorbed fragment (bp)
    const double DESORB_LAMBDA = log(2)/MEDIAN_DESORB_FRAG; // Prob of desorption of frag x bp is
                                                            // exp(-DESORB_LAMBDA*x)
    const double SIGMA = sqrt(0.300*1000); // SIGMA, units (bp)^(0.5)
    const double PI = 4.0*atan(1.0);

    // Derived constants from parameters above
    const double LOG_SQRT_2PI = 0.5*log(2*PI);
    const double LOG_1_MINUS_CUT_PROB = log(1 - CUT_PROB);
    const double LOG_CUT_PROB = log(CUT_PROB);
}

// ****************************************************************************
// Calculate probability that a reference fragment desorbs
double logProbDesorb(int fragLength)
{   
    using namespace probConstants;
    return -DESORB_LAMBDA*fragLength;
}

// ****************************************************************************
// Calculate probability that a reference fragment does not desorb
double logProbNoDesorb(int fragLength)
{
    using namespace probConstants;
    return log(1 - exp(-DESORB_LAMBDA*fragLength));
}


// ****************************************************************************
// Log of Poisson probability mass function f(x=k; lambda) = exp(-lambda)*lambda^k/k!
double logPoisson(int k, double lambda)
{
    assert(k >= 0);
    double log_kfactorial = log(boost::math::factorial<double>(k)); // Using boost for efficiency
    //double log_kfactorial = 0.0;
    //for (int i=k; i>1; i--) log_kfactorial += log(i);
    double result = -lambda + k*log(lambda) - log_kfactorial;
    assert(result < 0); // Log probability must be negative
    return result;
}

// ****************************************************************************
// Log of normal density function: f(x; mu, sigma) = (sqrt(2*PI)*sigma)^(-1) * exp(-0.5*((x-mu)/sigma)^2)
double logNormal(double x, double mu, double sigma)
{
    using namespace probConstants;
    double one_over_sigma = 1.0/sigma;
    double x_minus_mu = x - mu;
    double result = -0.5 * x_minus_mu * x_minus_mu * one_over_sigma * one_over_sigma;
    result += - LOG_SQRT_2PI - log(sigma);
    assert(result < 0); // Even though this is log of a prob. density, sigma 
                        // should be large enough so that the log is negative.
    return result;
}

// ****************************************************************************
// Calculate the Log Probability of an optical map alignment.

// Inputs:
// lRmap: Length of aggregate R-MAP fragment
// nRmap: number of fragments in aggregate R-MAP fragment
// lRef : Length of aggregate reference fragment
// nRef: number of fragments in aggregate reference fragment
// pRefFragLengths: pointer to entry in integer array with nRef fragment sizes
// Note: if pRefFragLengths not provided, we cannot include the probability that
// the nRef reference fragments do not desorb

// Output: Log probability
double calcOpMapLogProb(int lRmap, int nRmap, int lRef, int nRef, const int * pRefFragLengths = 0)
{
    using namespace probConstants;

    // Check if input indicates desporption (i.e. "lost" fragment from reference)
    if (nRmap==0)
    {
        // Can only calculate probability for one lost reference fragment
        // at a time if pRefFragLength not provided
        assert(nRef==1 || (nRef>1 && pRefFragLengths));
        double logProbFragDesorb = 0.0;
        if(pRefFragLengths)
        {
            const int * const pLastFrag = pRefFragLengths + nRef;
            int aggregateFragLength = 0;
            for (const int * pFragLength = pRefFragLengths;
                 pFragLength < pLastFrag;
                 pFragLength++)
            {
                logProbFragDesorb += logProbDesorb(*pFragLength);
                aggregateFragLength += *pFragLength;
            }
            assert(aggregateFragLength == lRef);
        }
        else 
            logProbFragDesorb = logProbDesorb(lRef);

        // Return probability that the rightmost cut in the RMAP was made and that
        // the fragment desorbed
        return logProbFragDesorb + LOG_CUT_PROB;
    }

    // Being here means that we are scoring an alignment of aggregate fragments
    // (and not desorption).
    assert(nRmap > 0 && nRef > 0);
    assert(lRmap > 0 && lRef > 0);
    int falseCuts = nRmap - 1;
    int missingCuts = nRef - 1;

    // Calculate log probability of false cuts in the R-MAP.
    double meanFalseCuts = FALSE_CUT_RATE*lRmap;
    double logProbFalseCuts = logPoisson(falseCuts, meanFalseCuts);

    // Calculate log probability of missing cuts in the R-MAP.
    // Note: this term includes the probability of the rightmost cut that defines boundary of
    //       the aggregate R-MAP fragment. 
    double logProbMissingCuts = falseCuts*LOG_1_MINUS_CUT_PROB + LOG_CUT_PROB;

    // Calculate the log probability of the sizing error of the R-MAP w.r.t the reference.
    double mu = lRef;
    double sigma = SIGMA*sqrt(lRef);
    double logProbDensitySizingError = logNormal(lRmap, mu, sigma);

    // If reference fragment lengths provided, calculate the probability that none of them
    // desorb.
    double logProbNoDesorption = 0.0;
    if(pRefFragLengths)
    {
        const int * const pLastFrag = pRefFragLengths + nRef;
        int aggregateFragLength = 0;
        for (const int * pFragLength = pRefFragLengths;
             pFragLength < pLastFrag;
             pFragLength++)
        {
            logProbNoDesorption += logProbNoDesorb(*pFragLength);
            aggregateFragLength += *pFragLength;
        }
        assert(aggregateFragLength == lRef);
    }
    return logProbFalseCuts + logProbMissingCuts + logProbDensitySizingError + logProbNoDesorption;
}


int main()
{
    using namespace std;

    int n1 = 1;
    int l1 = 100000;
    int n2 = 1;
    int l2 = 100000;
    int fragLengths[1] = {100000};
    double ans = calcOpMapLogProb(l1, n1 , l2, n2, fragLengths);
    cout << " ans: " << ans << endl;
    return 0;
}



