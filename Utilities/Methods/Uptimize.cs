using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace PeptidAce.Utilities.Methods
{
    public class UptimizeOptions
    {
        DBOptions dbOptions;
        Samples samples;
        Ace ace;
        public UptimizeOptions(DBOptions options, string fasta, Samples pSamples)
        {
            dbOptions = options;
            samples = pSamples;

            ace = new Ace(dbOptions, samples);
            ace.Preload(fasta, false, false);
            ace.PrepareQueries();
        }

        public void Run()
        {
            List<double> seed = new List<double>();

            seed.Add(dbOptions.dProduct);
            seed.Add(dbOptions.dPrecursor);
            seed.Add(dbOptions.dMatchingProductFraction);
            seed.Add(dbOptions.dMatchingProduct);
            seed.Add(dbOptions.dIntensityFraction);
            seed.Add(dbOptions.dIntensity);
            seed.Add(dbOptions.dProtein);
            seed.Add(dbOptions.dPeptideScore);
            seed.Add(dbOptions.dFragmentScore);
            List<double> bestOptions = Gradior.Minimize(FCT, seed, 0.00001, 0.0001, 1000, 1);//(long)PrecursorIntensityInCTrap * units.Count, 0.05);

            dbOptions.dProduct = seed[0];
            dbOptions.dPrecursor = seed[1];
            dbOptions.dMatchingProductFraction = seed[2];
            dbOptions.dMatchingProduct = seed[3];
            dbOptions.dIntensityFraction = seed[4];
            dbOptions.dIntensity = seed[5];
            dbOptions.dProtein = seed[6];
            dbOptions.dPeptideScore = seed[7];
            dbOptions.dFragmentScore = seed[8];
        }

        public double FCT(List<double> scoreWeights)
        {
            dbOptions.dProduct = scoreWeights[0];
            dbOptions.dPrecursor = scoreWeights[1];
            dbOptions.dMatchingProductFraction = scoreWeights[2];
            dbOptions.dMatchingProduct = scoreWeights[3];
            dbOptions.dIntensityFraction = scoreWeights[4];
            dbOptions.dIntensity = scoreWeights[5];
            dbOptions.dProtein = scoreWeights[6];
            dbOptions.dPeptideScore = scoreWeights[7];
            dbOptions.dFragmentScore = scoreWeights[8];
            foreach (Query query in ace.AllQueries)
            {
                query.psms.Clear();
                query.precursor.psms_AllPossibilities.Clear();
                if (query.precursor.psms != null)
                    query.precursor.psms.Clear();
            }//*/

            Result rez = ace.LaunchSearch(ace.AllQueries);
            int nbCorrectMatches = 0;

            int nbTarget = 0;
            foreach (Query query in rez.queries)
                if (query.psms.Count > 0)
                    if (query.Target)
                    {
                        nbTarget++;
                        if (query.psms[0].Peptide.Sequence.CompareTo(query.sample.nameColumn) == 0)
                            nbCorrectMatches++;
                    }

            return rez.queries.Count - nbCorrectMatches;// belowZero + tmpErrorOver + tmpErrorUnder;
        }
    }
}
