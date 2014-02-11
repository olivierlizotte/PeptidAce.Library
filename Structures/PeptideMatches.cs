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
    public class PeptideMatches : List<PeptideMatch>
    {
        public PeptideMatches() { }
        public PeptideMatches(IEnumerable<PeptideMatch> list) : base(list) { }

        public List<PeptideMatch> ComputeAtFDR(double desired_fdr)
        {
            this.Sort(PeptideMatch.DescendingOptimizedScoreComparison);

            int index = Utilities.Methods.FDRizer.Extract(this, desired_fdr);
            if (index >= 0)
                return this.GetRange(0, index + 1);
            else
                return new List<PeptideMatch>();
        }

        public static int CompareMatchingProducts(PeptideMatch left, PeptideMatch right)
        {
            return -left.BestPrecursor().OptimizedBestPsm(left.peptide).MatchingProducts.CompareTo(right.BestPrecursor().OptimizedBestPsm(right.peptide).MatchingProducts);
        }
        public static int CompareMatchingProductsFraction(PeptideMatch left, PeptideMatch right)
        {
            return -left.BestPrecursor().OptimizedBestPsm(left.peptide).MatchingProductsFraction.CompareTo(right.BestPrecursor().OptimizedBestPsm(right.peptide).MatchingProductsFraction);
        }
        public static int CompareMatchingIntensityFraction(PeptideMatch left, PeptideMatch right)
        {
            return -left.BestPrecursor().OptimizedBestPsm(left.peptide).MatchingIntensityFraction.CompareTo(right.BestPrecursor().OptimizedBestPsm(right.peptide).MatchingIntensityFraction);
        }
        public static int CompareProductScore(PeptideMatch left, PeptideMatch right)
        {
            return -left.BestPrecursor().OptimizedBestPsm(left.peptide).ProductScore.CompareTo(right.BestPrecursor().OptimizedBestPsm(right.peptide).ProductScore);
        }
        public static int CompareCumulPrecursorOptimizedScore(PeptideMatch left, PeptideMatch right)
        {
            return -left.CumulPrecursorOptimizedScore().CompareTo(right.CumulPrecursorOptimizedScore());
        }
        public static int CompareCumulPrecursorScore(PeptideMatch left, PeptideMatch right)
        {
            return -left.CumulPrecursorScore().CompareTo(right.CumulPrecursorScore());
        }
        public static int CompareBestPrecursorOptimizedScore(PeptideMatch left, PeptideMatch right)
        {
            return -left.BestPrecursorOptimizedScore().CompareTo(right.BestPrecursorOptimizedScore());
        }
        public static int CompareBestPrecursorScore(PeptideMatch left, PeptideMatch right)
        {
            return -left.BestPrecursorScore().CompareTo(right.BestPrecursorScore());
        }
        public static int CompareScore(PeptideMatch left, PeptideMatch right)
        {
            return -left.ProbabilityScore().CompareTo(right.ProbabilityScore());
        }
        public static int CompareScoreDescending(PeptideMatch left, PeptideMatch right)
        {
            return left.ProbabilityScore().CompareTo(right.ProbabilityScore());
        }
        public static int CompareNbCluster(PeptideMatch left, PeptideMatch right)
        {
            return -left.clusters.Count.CompareTo(right.clusters.Count);
        }
        public static int ComparePrecursorMassError(PeptideMatch left, PeptideMatch right)
        {
            return left.GetPrecursorMassError().CompareTo(right.GetPrecursorMassError());
        }
    }
}
