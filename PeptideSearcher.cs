/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using PeptidAce.Utilities;

namespace PeptidAce
{
    /// <summary>
    /// Peptide sequence match to at least one spectrum
    /// Holds the list of clusters associated to this sequence
    /// </summary>
    public class PeptideMatch : ITargetDecoy
    {
        public bool Decoy
        { get { return peptide.Decoy; } }

        public bool Target
        { get { return !peptide.Decoy; } }


        public List<Cluster> clusters;
        public Peptide peptide;
        public PeptideMatch()
        {
            clusters = new List<Cluster>();
        }
        public PeptideMatch(Peptide peptide)
        {
            this.peptide = peptide;
            clusters = new List<Cluster>();
        }
        public void AddOnlyOnce(Cluster cluster)
        {
            if (!clusters.Contains(cluster))
                clusters.Add(cluster);
        }

        public static int DescendingOptimizedScoreComparison(PeptideMatch left, PeptideMatch right)
        {
            return -(left.ProbabilityScore().CompareTo(right.ProbabilityScore()));
        }

        public static int DescendingPrecursorScoreComparison(PeptideMatch left, PeptideMatch right)
        {
            return -(left.CumulPrecursorScore().CompareTo(right.CumulPrecursorScore()));
        }

        public double ProbabilityScore()
        {
            double score = 0;
            foreach (Cluster cluster in clusters)
                score += /*(1 - score) * */cluster.ProbabilityScore(peptide);
            return score;
        }
        
        public double BestPrecursorScore()
        {
            double score = 0;
            foreach(Cluster cluster in clusters)
                foreach (clCondition condition in cluster.conditions)
                    if (condition != null)
                    foreach(clReplicate replicate in condition.replicates)
                        foreach (Precursor precursor in replicate.precursors)
                            if (precursor.ProbabilityScore(peptide) > score)
                                score = precursor.ProbabilityScore(peptide);
            return score;
        }
        
        public double BestPrecursorOptimizedScore()
        {
            double score = 0;
            foreach (Cluster cluster in clusters)
                foreach (clCondition condition in cluster.conditions)
                    if (condition != null)
                    foreach (clReplicate replicate in condition.replicates)
                        foreach (Precursor precursor in replicate.precursors)
                            if (precursor.ProbabilityScore() > score)
                                score = precursor.ProbabilityScore(peptide);
            return score;
        }

        public Precursor BestPrecursor(bool checkMods = false)
        {
            double score = 0;
            Precursor best = null;
            foreach (Cluster cluster in clusters)
            {
                Precursor tmp = cluster.OptimizedBestPrecursor(peptide, checkMods);
                if (tmp != null)
                {
                    double tmpScore = tmp.ProbabilityScore(peptide, checkMods);
                    if (tmpScore > score)
                    {
                        score = tmpScore;
                        best = tmp;
                    }
                }
            }
            return best;
        }
        
        public double CumulPrecursorScore()
        {
            double score = 0;
            foreach (Cluster cluster in clusters)
                foreach (clCondition condition in cluster.conditions)
                    if(condition != null)
                    foreach (clReplicate replicate in condition.replicates)
                        foreach (Precursor precursor in replicate.precursors)
                            score += precursor.ProbabilityScore(peptide);
            return score;
        }

        public double CumulPrecursorOptimizedScore()
        {
            double score = 0;
            foreach (Cluster cluster in clusters)
                foreach (clCondition condition in cluster.conditions)
                    if(condition != null)
                    foreach (clReplicate replicate in condition.replicates)
                        foreach (Precursor precursor in replicate.precursors)
                            score += precursor.ProbabilityScore(peptide);
            return score;
        }//*/

        public double GetPrecursorMassError()
        {
            double cumul = 0;
            int nbSeen = 0;
            foreach(Cluster cluster in clusters)
            {
                cumul += cluster.GetPrecursorMassError(peptide);
                nbSeen++;
            }
            return cumul / (double)nbSeen;
        }
    }

    /// <summary>
    /// Methods to parse the set of clusters and precursors to get a list of expressed peptide sequences
    /// </summary>
    public class PeptideSearcher
    {
        public DBOptions options;
        public PeptideSearcher(DBOptions options)
        {
            this.options = options;
        }
        
        public static int AscendingSpectrumNumberComparison(Cluster left, Cluster right)
        {
            return left.Rt().CompareTo(right.Rt());
        }
        
        public static int DescendingOptimizedScoreComparison(PeptideMatch left, PeptideMatch right)
        {
            return -left.ProbabilityScore().CompareTo(right.ProbabilityScore());
        }

        public PeptideMatches SearchAll(List<Cluster> clusters, List<Precursor> precursors, bool DiffByMod = true)
        {
            options.ConSole.WriteLine("Creating the list of peptide found...");
            Dictionary<string, PeptideMatch> peptideMatches = new Dictionary<string, PeptideMatch>();

            foreach (Precursor prec in precursors)
                foreach (PeptideSpectrumMatch psm in prec.psms_AllPossibilities)
                {
                    Peptide pep = psm.Peptide;
                    string seq = (DiffByMod ? pep.Sequence : pep.BaseSequence);
                    if (!peptideMatches.ContainsKey(seq))
                        peptideMatches.Add(seq, new PeptideMatch(pep));
                    else if (pep.Target)
                        peptideMatches[seq].peptide = pep;
                }

            int nbClusterNewSeq = 0;
            foreach (Cluster cl in clusters)
                foreach (clCondition condition in cl.conditions)
                    if(condition != null)
                    foreach (clReplicate replicate in condition.replicates)
                        foreach (Precursor precursor in replicate.precursors)
                            foreach (PeptideSpectrumMatch psm in precursor.psms_AllPossibilities)
                            {
                                string seq = (DiffByMod ? psm.Peptide.Sequence : psm.Peptide.BaseSequence);
                                if (peptideMatches.ContainsKey(seq))
                                    peptideMatches[seq].AddOnlyOnce(cl);
                            }

            PeptideMatches TotalList = new PeptideMatches(peptideMatches.Values);
            options.ConSole.WriteLine(TotalList.Count + " distinct peptides (based on sequence" + (DiffByMod ? " and modifications)" : ")") + "    [" + nbClusterNewSeq + " added from clusters]");
            return TotalList;
        }

        public PeptideMatches SearchClusters(List<Cluster> clusters, List<Precursor> precursors, bool DiffByMod = true)
        {
            options.ConSole.WriteLine("Creating the list of peptide found...");
            Dictionary<string, PeptideMatch> peptideMatches = new Dictionary<string, PeptideMatch>();

            foreach (Precursor prec in precursors)
                foreach (PeptideSpectrumMatch psm in prec.OptimizedBestPsms())
                {
                    Peptide pep = psm.Peptide;
                    string seq = (DiffByMod ? pep.Sequence : pep.BaseSequence);
                    if (!peptideMatches.ContainsKey(seq))
                        peptideMatches.Add(seq, new PeptideMatch(pep));
                    else if (pep.Target)
                        peptideMatches[seq].peptide = pep;
                }

            int nbClusterNewSeq = 0;
            foreach (Cluster cl in clusters)
                foreach (clCondition condition in cl.conditions)
                    if(condition != null)
                    foreach (clReplicate replicate in condition.replicates)
                        if(replicate != null)
                        foreach (Precursor precursor in replicate.precursors)
                            foreach (PeptideSpectrumMatch psm in precursor.OptimizedBestPsms())
                            {
                                string seq = (DiffByMod ? psm.Peptide.Sequence : psm.Peptide.BaseSequence);
                                if (peptideMatches.ContainsKey(seq))
                                    peptideMatches[seq].AddOnlyOnce(cl);
                                else
                                    nbClusterNewSeq++;
                            }

            PeptideMatches TotalList = new PeptideMatches(peptideMatches.Values);
            options.ConSole.WriteLine(TotalList.Count + " distinct peptides (based on sequence" + (DiffByMod ? " and modifications)" : ")") + "    [" + nbClusterNewSeq + " not unmapped clusters]");
            return TotalList;
        }

        public PeptideMatches Search(List<Cluster> clusters, List<Precursor> precursors, bool DiffByMod = true)
        {
            options.ConSole.WriteLine("Creating the list of peptide found...");
            Dictionary<string, PeptideMatch> peptideMatches = new Dictionary<string, PeptideMatch>();            
         
            foreach(Precursor prec in precursors)
                foreach (PeptideSpectrumMatch psm in prec.OptimizedBestPsms())
                {
                    Peptide pep = psm.Peptide;
                    string seq = (DiffByMod ? pep.Sequence : pep.BaseSequence);
                    if (!peptideMatches.ContainsKey(seq))
                        peptideMatches.Add(seq, new PeptideMatch(pep));
                    else if (pep.Target)
                        peptideMatches[seq].peptide = pep;
                }

            int nbClusterNewSeq = 0;
            foreach(Cluster cl in clusters)
                foreach (clCondition condition in cl.conditions)
                    if (condition != null)
                    foreach (clReplicate replicate in condition.replicates)
                        foreach (Precursor precursor in replicate.precursors)
                            foreach (PeptideSpectrumMatch psm in precursor.OptimizedBestPsms())
                            {
                                string seq = (DiffByMod ? psm.Peptide.Sequence : psm.Peptide.BaseSequence);
                                if (peptideMatches.ContainsKey(seq))
                                    peptideMatches[seq].AddOnlyOnce(cl);
                            }

            PeptideMatches TotalList = new PeptideMatches(peptideMatches.Values);
            options.ConSole.WriteLine(TotalList.Count + " distinct peptides (based on sequence" + (DiffByMod ? " and modifications)" : ")") + "    [" + nbClusterNewSeq + " added from clusters]");
            return TotalList;
        }

        public static void Export(string filename, List<PeptideMatch> peptides)
        {
            vsCSVWriter writer = new vsCSVWriter(filename);
            writer.AddLine("Sequence,Variable Modification,Score,Decoy,Precursor Mass Error");
            foreach (PeptideMatch pm in peptides)
                writer.AddLine(pm.peptide.BaseSequence + "," + pm.peptide.Sequence + "," + pm.ProbabilityScore() + "," + pm.peptide.Decoy + "," + pm.GetPrecursorMassError());
            writer.WriteToFile();
        }
    }
}
