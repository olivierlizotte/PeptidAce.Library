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
    /// Replicates (lists precursors for a given replicate (aka all fractions) 
    /// </summary>
    public class clReplicate
    {
        public int id;
        public List<Precursor> precursors;
        public clReplicate()
        {
            precursors = new List<Precursor>();
        }
        public void Add(Precursor precursor)
        {
            precursors.Add(precursor);
        }
        public double ProbabilityScore(Peptide peptide)
        {
            double score = 0;
            foreach (Precursor precursor in precursors)
                score += /*(1 - score) * */precursor.ProbabilityScore(peptide, true);
            
            return score;
        }
        public double Rt()
        {
            double tmpRt = 0;
            int nbRep = 0;
            foreach (Precursor precursor in precursors)
            {
                nbRep++;
                tmpRt += precursor.Track.RT;
            }
            return tmpRt / (double)nbRep;
        }
        public double Mz()
        {
            double tmpMz = 0;
            int nbRep = 0;
            foreach (Precursor precursor in precursors)
            {
                nbRep++;
                tmpMz += precursor.Track.MZ;
            }
            return tmpMz / (double)nbRep;
        }
    }

    /// <summary>
    /// Condition (list of replicates) were a replicate is found
    /// </summary>
    public class clCondition 
    {
        public int id;
        public List<clReplicate> replicates;
        public clCondition()
        {
            replicates = new List<clReplicate>();
        }
        public clCondition(int nbReplicates)
        {
            replicates = new List<clReplicate>();
            for(int i = 0; i < nbReplicates; i++)
                replicates.Add(null);
        }
        public void Add(Precursor precursor)
        {
            if (replicates[precursor.sample.PROJECT.REPLICATE - 1] == null)
                replicates[precursor.sample.PROJECT.REPLICATE - 1] = new clReplicate();
            replicates[precursor.sample.PROJECT.REPLICATE - 1].Add(precursor);
        }
        public double Reproducibility()
        {
            int nbElem = 0;
            foreach (clReplicate replicate in replicates)
                if (replicate != null)
                    nbElem++;
            return nbElem / (double) replicates.Count;
        }
        public double ProbabilityScore(Peptide peptide)
        {
            double score = 0;
            foreach (clReplicate replicate in replicates)
                if (replicate != null)
                    score += /*(1 - score) * */replicate.ProbabilityScore(peptide);
            return score;
        }
        public double Rt()
        {
            double tmpRt = 0;
            int nbRep = 0;
            foreach (clReplicate replicate in replicates)
                if (replicate != null)
                {
                    nbRep++;
                    tmpRt += replicate.Rt();
                }
            return tmpRt / (double)nbRep;
        }
        public double Mz()
        {
            double tmpMz = 0;
            int nbRep = 0;
            foreach (clReplicate replicate in replicates)
                if (replicate != null)
                {
                    nbRep++;
                    tmpMz += replicate.Mz();
                }
            return tmpMz / (double)nbRep;
        }
    }

    /// <summary>
    /// A cluster is a group of ions, corresponding to the same peptide sequence, seen accross multiple samples
    /// </summary>
    public class Cluster
    {
        public List<clCondition> conditions;
        private Dictionary<int, List<int>> samples;
        public Cluster()
        {
            conditions = new List<clCondition>();
        }
        public Cluster(Dictionary<int, List<int>> samples)
        {
            this.conditions = new List<clCondition>();
            foreach (int c in samples.Keys)
                this.conditions.Add(null);
            this.samples = samples;
        }
        
        public double Reproducibility()
        {
            //Reproducibility is average number of replicates for conditions where seen
            double cumul = 0;
            int nbSeen = 0;
            foreach (clCondition condition in conditions)
                if (condition != null)
                {
                    double coverage = condition.Reproducibility();
                    if (coverage > 0)
                    {
                        cumul += coverage;
                        nbSeen++;
                    }
                }

            return (cumul / (double)nbSeen);
        }//*/

        public double GetPrecursorMassError(Peptide peptide)
        {
            //Reproducibility is average number of replicates for conditions where seen
            double cumul = 0;
            int nbSeen = 0;
            foreach (clCondition condition in conditions)
                if(condition != null)
                foreach(clReplicate replicate in condition.replicates)
                    foreach (Precursor precursor in replicate.precursors)
                    {
                        cumul += peptide.MonoisotopicMass - precursor.Mass;
                        nbSeen++;
                    }

            return (cumul / (double)nbSeen);
        }//*/
        public void Add(Precursor precursor)
        {
            if (conditions[precursor.sample.PROJECT.CONDITION - 1] == null)
                conditions[precursor.sample.PROJECT.CONDITION - 1] = new clCondition(samples[precursor.sample.PROJECT.CONDITION].Count);
            conditions[precursor.sample.PROJECT.CONDITION - 1].Add(precursor);
            //score = -1;
        }
        public double ProbabilityScore(Peptide peptide)
        {
            double score = 0;
            int nbCond = 0;
            foreach (clCondition condition in conditions)
                if (condition != null)
                {
                    double tmp = condition.ProbabilityScore(peptide);
                    if (tmp > 0)
                    {
                        score += tmp;
                        nbCond++;
                    }
                }
            if (nbCond > 0)
                return score;// / (double)nbCond;
            else
                return 0.0;
        }
        public double Rt()
        {
            double tmpRt = 0;
            int nbRep = 0;
            foreach (clCondition condition in conditions)
                if (condition != null)
                {
                    nbRep++;
                    tmpRt += condition.Rt();
                }
            return tmpRt / (double)nbRep;
        }
        public double Mz()
        {
            double tmpMz = 0;
            int nbRep = 0;
            foreach (clCondition condition in conditions)
                if (condition != null)
                {
                    nbRep++;
                    tmpMz += condition.Mz();
                }
            return tmpMz / (double)nbRep;
        }

        public List<PeptideSpectrumMatch> OptimizedBestPsms()
        {
            List<PeptideSpectrumMatch> matches = new List<PeptideSpectrumMatch>();
            foreach (clCondition condition in conditions)
                if(condition != null)
                foreach (clReplicate replicate in condition.replicates)
                    foreach (Precursor precursor in replicate.precursors)
                        foreach (PeptideSpectrumMatch psm in precursor.OptimizedBestPsms())
                            if (matches.Count == 0 || matches[0].ProbabilityScore() == psm.ProbabilityScore())
                                matches.Add(psm);
                            else
                                if (matches[0].ProbabilityScore() < psm.ProbabilityScore())
                                {
                                    matches.Clear();
                                    matches.Add(psm);
                                }
            return matches;
        }

        public Peptide ComputeBestPeptide()
        {            
            Dictionary<Peptide, double> peptides = new Dictionary<Peptide, double>();            
            foreach (clCondition condition in conditions)
                if(condition != null)
                foreach (clReplicate replicate in condition.replicates)
                    foreach (Precursor precursor in replicate.precursors)
                        foreach (PeptideSpectrumMatch psm in precursor.OptimizedBestPsms())
                            if(!peptides.ContainsKey(psm.Peptide))
                                peptides.Add(psm.Peptide, ProbabilityScore(psm.Peptide));
            double best = -1;
            Peptide bestPeptide = null;
            foreach(Peptide key in peptides.Keys)
                if (peptides[key] > best || (peptides[key] == best && bestPeptide.Decoy))
                {
                    best = peptides[key];
                    bestPeptide = key;
                }
            return bestPeptide;
        }
        /*
        public PeptideSpectrumMatch OptimizedBestPsm(Peptide peptide, bool checkMods)
        {
            PeptideSpectrumMatch bestPsm = null;
            foreach (clCondition condition in conditions)
                foreach (clReplicate replicate in condition.replicates)
                    foreach (Precursor precursor in replicate.precursors)
                    {
                        PeptideSpectrumMatch match = precursor.OptimizedBestPsm(peptide, checkMods);
                        if(bestPsm == null || bestPsm.OptimizedScore() < match.OptimizedScore())
                            bestPsm = match;
                    }
            return bestPsm;
        }//*/

        public Precursor OptimizedBestPrecursor(Peptide peptide, bool checkMods = false)
        {
            double score = 0;
            Precursor best = null;
            foreach (clCondition condition in conditions)
                if(condition != null)
                foreach (clReplicate replicate in condition.replicates)
                    foreach (Precursor precursor in replicate.precursors)
                    {
                        double tmpScore = precursor.ProbabilityScore(peptide, checkMods);
                        if (tmpScore > score)
                        {
                            score = tmpScore;
                            best = precursor;
                        }
                    }
            return best;
        }
    }

    /// <summary>
    /// Methods to cluster (aggregate) precursors seen in more than one sample
    /// </summary>
    public class MSSearcher
    {
        public DBOptions options;
        public Dictionary<int, List<int>> samples = new Dictionary<int, List<int>>();
        public Precursors precursors = new Precursors();
        
        public MSSearcher(DBOptions options, Samples project)
        {
            this.options = options;
            
            foreach(Sample sample in project)
            {
                if (!samples.ContainsKey(sample.PROJECT.CONDITION))
                    samples.Add(sample.PROJECT.CONDITION, new List<int>());
                if (!samples[sample.PROJECT.CONDITION].Contains(sample.PROJECT.REPLICATE))
                    samples[sample.PROJECT.CONDITION].Add(sample.PROJECT.REPLICATE);
            }
        }

        public void CumulPsm(List<Precursor> matches)
        {
            foreach (Precursor precursor in matches)
            {
                Sample entry = precursor.sample;
                precursors.Add(precursor);//Remove blank psm (unmatched spectrum/query duos)

             /*   if (!samples.ContainsKey(entry.PROJECT.CONDITION))
                    samples.Add(entry.PROJECT.CONDITION, new List<int>());
                if (!samples[entry.PROJECT.CONDITION].Contains(entry.PROJECT.REPLICATE))
                    samples[entry.PROJECT.CONDITION].Add(entry.PROJECT.REPLICATE);//*/
            }
        }

        private int NbCommonSequences(Precursor a, Precursor b)
        {
            int common = 0;
            foreach (PeptideSpectrumMatch psmA in a.psms)
            {
                bool seen = true;
                foreach (PeptideSpectrumMatch psmB in b.psms)
                    if (psmA.Peptide.IsSamePeptide(psmB.Peptide, true))
                    {
                        seen = true;
                        break;
                    }
                if (seen)
                    common++;
            }
            return common;
        }

        private double Score(Precursor a, Precursor b)
        {
            if (a.sample != b.sample && NbCommonSequences(a, b) > 0)
            {
                //Zero based score
                //double tmp = Math.Abs(MassTolerance.CalculateMassError(a.Track.MZ, b.Track.MZ, MassToleranceUnits.ppm) / options.MzTol);
                double tmp = Math.Abs(Numerics.MzDifference(a.Track.MZ, b.Track.MZ, options.precursorMassTolerance.Units)) / options.precursorMassTolerance.Value;
                if (tmp < 1)
                {
                    tmp = 0.2 * (2 * tmp +
                                 Math.Abs(a.Track.RT - b.Track.RT) / options.ComputedRetentionTimeDiff + //TODO check if it is in seconds?
                                (a.Charge == b.Charge ? 0 : 1) +
                                 0.1 * Math.Abs(Math.Log10(a.Track.INTENSITY) - Math.Log10(b.Track.INTENSITY)));
                    if (tmp < 1)
                        return 1 - tmp;
                }
            }
            return 0;
        }
        
        public static int DescendingScoreComparison(Query left, Query right)
        {
            return - left.Score.CompareTo(right.Score);
        }

        public static int DescendingProteinScoreComparison(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.ProteinScore.CompareTo(right.ProteinScore);
        }

        public List<Cluster> Search(Precursors precursors, bool runCluster)
        {
            options.ConSole.WriteLine("Grouping precursors based on common features...");
            precursors.Sort(Precursor.CompareProbabilityScore);

            List<Cluster> clusters = new List<Cluster>();
            bool[] done = new bool[precursors.Count];
            for(int i = 0; i < done.Length; i++)
                done[i] = false;

            //Step 1 : Regroup psms based on mz/rt/intensity/Sequence proximity score (ProteoProfile Code)
            for (int i = 0; i < precursors.Count; i++)
            {
                if(!done[i])
                {
                    Cluster group = new Cluster(samples);
                    group.Add(precursors[i]);
                    if (runCluster)
                    {
                        for (int j = i + 1; j < precursors.Count; j++)
                        {
                            if (!done[j] && precursors[i].sample != precursors[j].sample)
                            {
                                double score = Score(precursors[i], precursors[j]);
                                //TODO Implement ProteoProfile Clustering algorithm, or anything on the litterature, as long as its backed by the scoring function
                                if (score > 0.75)//TODO Should we put a threshold here? Can it be computed dynamically?
                                {
                                    group.Add(precursors[j]);
                                    done[j] = true;
                                }
                            }
                        }
                    }
                    clusters.Add(group);
                }
            }
            options.ConSole.WriteLine("Created " + clusters.Count + " clusters");
            return clusters;
            //TODO I should not use psms in more than one cluster...
        }

        public static void Export(string filename, List<Precursor> precursors)
        {
            vsCSVWriter writer = new vsCSVWriter(filename);
            writer.AddLine("Index.Mz,Rt,Precursor Mz,Charge,Most Intense Charge,Precursor Mass,Peptide Mass,Sequence,Modified Sequence,Precursor Score,Product Score,Intensity Score,Final Score,Precursor Mass Error,Decoy?,Protein Score");

            foreach (Precursor precursor in precursors)
            {
                string line = precursor.INDEX + "," + precursor.Track.RT + "," + precursor.Track.MZ + "," + precursor.Charge + "," + precursor.GetMostIntenseCharge() + "," + precursor.Mass + ",";
                PeptideSpectrumMatch match = precursor.OptimizedBestPsm();
                if (match != null)
                    line += match.Peptide.MonoisotopicMass + "," + match.Peptide.BaseSequence + "," + match.Peptide.Sequence + "," + match.PrecursorScore + "," + match.ProductScore + "," + match.IntensityScore + "," + precursor.ProbabilityScore(match.Peptide) + "," +
                        match.PrecursorMzError + "," + match.Decoy + "," + match.ProteinScore;
                writer.AddLine(line);
            }

            writer.WriteToFile();
        }

        public static void Export(string filename, IEnumerable<Query> queries)
        {
            vsCSVWriter writer = new vsCSVWriter(filename);
            writer.AddLine("Spectrum Precursor Mz,Rt,Charge,BaseSequence,Sequence,Nb Matched Fragments,Fragment Score,Precursor Score,Product Score,Intensity Score,Final Score,Precursor Mass Error,Decoy?,Protein Score");

            foreach (Query query in queries)
            {
                string line = query.spectrum.PrecursorMZ + "," + query.precursor.Track.RT + "," + query.precursor.Charge + ",";
                PeptideSpectrumMatch match = query.precursor.OptimizedBestPsm();
                if(match != null)
                    line += match.Peptide.BaseSequence + "," + match.Peptide.Sequence + "," + match.MatchingProducts + "," + match.FragmentScore + "," +  match.PrecursorScore + "," + match.ProductScore + "," + match.IntensityScore + "," + query.ScoreFct(match.Peptide) + "," +
                        match.PrecursorMzError + "," + match.Decoy + "," + match.ProteinScore;
                writer.AddLine(line);
            }

            writer.WriteToFile();
        }

        public static void Export(string filename, List<PeptideSpectrumMatch> psms)
        {
            vsCSVWriter writer = new vsCSVWriter(filename);
            writer.AddLine("Mz,Rt,Charge,Sequence,Modifications,Precursor Score,Product Score,Intensity Score,Final Score,Precursor Mass Error,Decoy?,Protein Score");
            
            foreach (PeptideSpectrumMatch psm in psms)
                writer.AddLine(psm.Query.precursor.Track.MZ +","+ psm.Query.spectrum.RetentionTimeInMin + ","+psm.Query.precursor.Charge +
                    "," + psm.Peptide.BaseSequence + "," + psm.Peptide.Sequence + "," + psm.PrecursorScore + "," + psm.ProductScore + "," + psm.IntensityScore + "," + psm.ProbabilityScore() + "," +
                    psm.PrecursorMzError + "," + psm.Decoy + "," + psm.ProteinScore);
            
            writer.WriteToFile();
        }
    }
}
