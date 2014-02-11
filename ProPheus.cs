/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Collections.Concurrent;
using System.IO;
using System.Threading;
using System.Threading.Tasks;

namespace PeptidAce
{

    /// <summary>
    /// Regroups method used to map digested peptide sequences to one or more spectrum
    /// </summary>
    public class ProPheus
    {
        public DBOptions options;
        public ProPheus(DBOptions options)
        {
            this.options = options;
        }

        /// <summary>
        /// Creates a PSM, which computes the score. Only top scorers are saved
        /// </summary>
        /// <param name="query"></param>
        /// <param name="modified_peptide"></param>
        private void ComputePSMs(Query query, Peptide modified_peptide)
        {
            PeptideSpectrumMatch psm = new PeptideSpectrumMatch(query, modified_peptide, options);

            if (psm.MatchingProducts >= options.NbMinProducts)
            {
                if (query.psms.Count < options.NbPSMToKeep)
                    query.psms.Add(psm);
                else if (query.psms[options.NbPSMToKeep-1].ProbabilityScore() < psm.ProbabilityScore())
                {
                    for (int i = 0; i < query.psms.Count; i++)
                        if (query.psms[i].ProbabilityScore() <= psm.ProbabilityScore())
                        {
                            query.psms.Insert(i, psm);
                            break;
                        }
                    if (query.psms.Count > options.NbPSMToKeep)
                        query.psms.RemoveAt(options.NbPSMToKeep-1);
                }
            }
        }

        /// <summary>
        /// Maps all peptide sequences to potential precursors and spectrum
        /// </summary>
        /// <param name="queries"></param>
        /// <param name="fittingPeptides"></param>
        /// <param name="previousProteins"></param>
        /// <returns></returns>
        public Precursors Search(Queries queries, IEnumerable<Tuple<Peptide, int>> fittingPeptides)
        {
            queries.dbOptions.ConSole.WriteLine("Mapping " + queries.Count + " queries to the digested proteome ... ");

            long nbQueryConsidered = 0;
            Parallel.ForEach<Tuple<Peptide, int>>(fittingPeptides, (Tuple<Peptide, int> hit) =>            
            {
                int indexPrecursor = hit.Item2;
                double maximumMass = MassTolerance.MzTop(hit.Item1.MonoisotopicMass, options.precursorMassTolerance);
                double minimumMass = MassTolerance.MzFloor(hit.Item1.MonoisotopicMass, options.precursorMassTolerance);
                
                if (indexPrecursor < queries.Count && queries[indexPrecursor].precursor.Mass >= minimumMass)
                {
                    while (indexPrecursor < queries.Count && queries[indexPrecursor].precursor.Mass <= maximumMass)
                    {

                        lock (queries[indexPrecursor].psms)
                        {
                            //Target (or decoy with enzyme digests)
                            ComputePSMs(queries[indexPrecursor], hit.Item1);

                            //Add Decoy if NoEnzyme digest
                            if (options.DecoyFusion)
                                ComputePSMs(queries[indexPrecursor], hit.Item1.Reverse());

                            indexPrecursor++;
                        }
                    }

                    nbQueryConsidered += indexPrecursor - hit.Item2;
                }
                else
                    options.ConSole.WriteLine("WTF####");
            });

            //Push PSMs to Precursor
            foreach (Query query in queries)
                query.precursor.psms_AllPossibilities.AddRange(query.psms);//No MERGE

            //PeptideSpectrumMatches allPsms = new PeptideSpectrumMatches();
            int nbAssignedPrecursor = 0;
            foreach (Precursor precursor in queries.Precursors)
                if (precursor.psms_AllPossibilities.Count > 0)
                    nbAssignedPrecursor++;

            int nbAssignedQuery = 0;
            foreach (Query query in queries)
                if(query.psms.Count > 0)
                    nbAssignedQuery++;
            options.ConSole.WriteLine(nbAssignedQuery + " queries matched [" + nbAssignedPrecursor + " precursors] out of " + nbQueryConsidered + " psm computed");
                        
            return queries.Precursors;
        }
    }
}