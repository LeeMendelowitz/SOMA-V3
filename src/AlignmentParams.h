#ifndef ALIGNMENTPARAMS_H
#define ALIGNMENTPARAMS_H
// Parameters used for an alignment
// This is an input to the match function
// TODO: Consider making classes to wrap the functions that create the scoreMatrices
// and compute the scores. This would avoid the need to constantly pass AlignmentParams
// around, or to consistently reference global variables.
class AlignmentParams
{
    public:
    AlignmentParams(double C_r_contig_in, double C_r_optical_in, double C_sigma_in,
                    double sigma2_in, int maxChunkMissesQuery_in, int maxChunkMissesReference_in,
                    double smallFrag_in, double smallFragSlope_in,
                    double H_in, double T_in)
    {
        C_r_contig = C_r_contig_in;
        C_r_optical = C_r_optical_in;
        chi2Max = C_sigma_in * C_sigma_in;
        sigma2 = sigma2_in;
        maxChunkMissesQuery = maxChunkMissesQuery_in;
        maxChunkMissesReference = maxChunkMissesReference_in;

        smallFrag = smallFrag_in;
        smallFragSlope = smallFragSlope_in;

        // Shape parameters for parabola for local scoring function
        // Consider creating a class for 
        T2 = T_in * T_in;
        A = T2 > 0 ? H_in / (T2) : 0;

    }

    double C_r_contig; // Cost of missing a restriction site in the contig insilico map
    double C_r_optical; // Cost of missing a restriction site in the optical map
    double chi2Max; // The maximum allowed length difference for aligned fragments, in variance units.
    double sigma2; // The parameter for sigma^2 for computing the variance
    int maxChunkMissesQuery; // Maximum number of unaligned sites inside an alignment block
    int maxChunkMissesReference; // Maximum number of unaligned sites inside an alignment block

    // Params for penalizing gapped contig fragments
    // or missed contig sites in small fragments
    double smallFrag;
    double smallFragSlope;

    // Params for local alignment
    double A;
    double T2;

};

#endif
