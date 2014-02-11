/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Text;
using System.Collections.Generic;
using PeptidAce.Utilities;

namespace PeptidAce
{
    public class PeptideSpectrumMatch : ITargetDecoy
    {
        public double PeptideScore;
        public double ProteinScore;

        public double highestFragmentIntensity;
        
        public Query Query { get; set; }

        public Peptide Peptide { get; set; }

        public double PrecursorMzError { get; set; }

        public int MatchingProducts { get; set; }
        public double MatchingWeightedProducts { get; set; }
        public double MatchingProductsProbability { get; set; }
        public List<ProductMatch> AllProductMatches{ get; set; }

        public int TotalTheoreticalProducts { get; set; }
        public double TotalWeightedProducts { get; set; }

        public double MatchingProductsFraction { get; set; }

        public double MatchingIntensity { get; set; }

        public double MatchingIntensityFraction { get; set; }

        public double MorpheusScore { get;  set; }
        public double PrecursorScore { get;  set; }
        public double ProductScore { get;  set; }
        public double IntensityScore { get; set; }
        public double FragmentScore { get; set; }

        public double[] ComputeAACoverage()
        {
            double[] aaScoreArray = new double[Peptide.Length];
            foreach (ProductMatch match in AllProductMatches)
            {
                if (match.Fragment == null)
                    aaScoreArray[match.fragmentPos] += match.normalizedIntensity;
                else if (match.Fragment.IsReverse)
                    for (int i = match.fragmentPos - 1; i < aaScoreArray.Length; i++)
                        aaScoreArray[i] += match.normalizedIntensity;
                else
                    for (int i = match.fragmentPos - 1; i >= 0; i--)
                        aaScoreArray[i] += match.normalizedIntensity;
            }
            return aaScoreArray;
        }

        private double ComputeFragmentScore()
        {
            double[] aaScoreArray = ComputeAACoverage();
            double score = 0.0;
            for (int i = 0; i < aaScoreArray.Length; i++)
                if (aaScoreArray[i] <= 1)
                    score += aaScoreArray[i];
                else
                    if (aaScoreArray[i] < 2)
                        score += 2 - aaScoreArray[i];
            return score / (double) aaScoreArray.Length;
        }

        //TODO Check if those var are common with Precursor and remove zeroed out var
        public double ProbabilityScore()
        {
            DBOptions options = Query.options;
            //double score = FragmentScore() * PrecursorScore;// +ProductScore;
            double score = options.dIntensity * MatchingIntensity +
                                options.dIntensityFraction * MatchingIntensityFraction +
                                options.dProduct * ProductScore +
                                options.dPrecursor * PrecursorScore +
                                options.dMatchingProductFraction * MatchingProductsFraction +
                                options.dMatchingProduct * MatchingWeightedProducts +
                                options.dProtein * ProteinScore +
                                options.dPeptideScore * PeptideScore +
                                options.dFragmentScore * FragmentScore;//*/
            //double score = FragmentScore;
            return score;
        }

        public bool Decoy
        {
            get { return Peptide.Decoy; }
        }

        public bool Target
        {
            get { return !Decoy; }
        }

        private static readonly ProductTypes PRODUCT_TYPES = ProductTypes.Instance;

        public PeptideSpectrumMatch()
        {
            AllProductMatches = new List<ProductMatch>();
        }

        public PeptideSpectrumMatch(Query query, Peptide peptide, DBOptions options)
        {
            Query = query;
            highestFragmentIntensity = query.spectrum.MostIntensePeak;
            Peptide = peptide;
            
            UpdatePrecursor(options);

            Initialize(options, GetProductMZs(options, query.spectrum.Peaks));
        }

        public void UpdatePrecursor(DBOptions options)
        {
            //PrecursorMassErrorDa = (spectrum. - spectrum.IsotopicShift) - peptide.MonoisotopicMass;
            PrecursorMzError = Numerics.MzDifference(Numerics.MZFromMass(Peptide.MonoisotopicMass, Query.precursor.Charge), Query.precursor.Track.MZ, options.precursorMassTolerance.Units);
            if (PrecursorMzError > options.precursorMassTolerance.Value)
                PrecursorMzError = options.precursorMassTolerance.Value;

            PrecursorScore = (options.precursorMassTolerance.Value - Math.Abs(PrecursorMzError)) / options.precursorMassTolerance.Value;
        }

        public IEnumerable<ProductMatch> GetProductMZs(DBOptions options, List<MsMsPeak> peaks)//, List<double> theoretical_product_mzs = null)
        {
            // speed optimizations
            //double[] experimental_masses = Query.spectrum.Masses;
            //List<MSPeak> peaks = Query.spectrum.Peaks;
            int num_experimental_peaks = peaks.Count;
            TotalTheoreticalProducts = 0;
            TotalWeightedProducts = 0;
            //New version that should include charged ions
 
            foreach (ProductMatch matchTheo in options.fullFragment.ComputeFragments(Peptide.GetMasses(), Query.precursor.Charge))
            {
                TotalTheoreticalProducts++;
                TotalWeightedProducts += matchTheo.weight;
                double massDiff = options.productMassTolerance.Value;
                double bestMz = -1;
                double bestInt = 0;

                foreach (int index in Query.spectrum.GetIndexOfMZInRange(matchTheo.theoMz, options.productMassTolerance))
                {
                    if (peaks[index].Charge <= 0 || peaks[index].Charge == matchTheo.charge)
                    {

                        double diff = Numerics.CalculateMassError(peaks[index].MZ, matchTheo.theoMz, options.productMassTolerance.Units);
                        if (Math.Abs(diff) < options.productMassTolerance.Value)
                        {
                            if (Math.Abs(diff) < Math.Abs(massDiff))//TODO Priority to intensity, or precision?
                            {
                                massDiff = diff;
                                bestMz = peaks[index].MZ;
                            }
                            bestInt += peaks[index].Intensity;
                        }
                    }
                }
                if (bestMz >= 0)
                {/*
                    double secondIsoIntensity = 0.0;
                    //Check for second isotopic peak
                    double sIsoMass = matchTheo.theoMz + Numerics.IsotopicMassShift(1, matchTheo.charge);
                    foreach (int index in Query.spectrum.GetIndexOfMZInRange(matchTheo.theoMz + Numerics.IsotopicMassShift(1, matchTheo.charge), options.productMassTolerance))
                    {
                        if (peaks[index].Charge <= 0 || peaks[index].Charge == matchTheo.charge)
                        {
                            double diff = Numerics.CalculateMassError(peaks[index].MZ, sIsoMass, options.productMassTolerance.Units);
                            //double diff = matchTheo.theoMz - peaks[index].MZ;// experimental_masses[index];//TODO DALTON ONLY : add product mass tolerance unit test
                            if (Math.Abs(diff) < options.productMassTolerance.Value)
                                secondIsoIntensity += peaks[index].Intensity;
                        }
                    }
                    if (secondIsoIntensity >= bestInt * 0.0107)//Second Isotope should be present... and should not account for less than 0.//*/
                    {
                        ProductMatch pMatch = new ProductMatch();
                        pMatch.weight = matchTheo.weight;
                        pMatch.theoMz = matchTheo.theoMz;// Utilities.MZFromMzSingleCharge(theoMass, charge);
                        pMatch.obsMz = bestMz;// experimental_masses[bestIndex];
                        pMatch.mass_diff = massDiff;
                        pMatch.obsIntensity = bestInt;// Intensities[bestIndex];
                        pMatch.charge = matchTheo.charge;// peaks[bestIndex].Charge;
                        pMatch.Fragment = matchTheo.Fragment;
                        pMatch.fragmentPos = matchTheo.fragmentPos;
                        pMatch.normalizedIntensity = pMatch.obsIntensity / (Query.spectrum.InjectionTime * Query.spectrum.PrecursorIntensityPerMilliSecond);
                        yield return pMatch;
                    }
                    //break;
                }
            }
        }

        public void Initialize(DBOptions options, IEnumerable<ProductMatch> productMZs)
        {
            //List<double> theoretical_product_mz = Peptide.CalculateAllProductMz(PRODUCT_TYPES[Query.spectrum.FragmentationMethod], Query.precursor);
            //TotalProducts = theoretical_product_mz.Count;
            MatchingProducts    = 0;
            MatchingIntensity   = 0.0;
            double cumulDiff    = 0;
            ProductScore        = 0;
            List<ProductMatch> cumulMatch = new List<ProductMatch>();
            MatchingWeightedProducts = 0;
            foreach (ProductMatch match in productMZs)
            {
                MatchingProducts++;
                MatchingWeightedProducts  += match.weight;
                MatchingIntensity += match.obsIntensity;
                cumulDiff += Math.Abs(match.mass_diff) * match.weight;
                
                cumulMatch.Add( match );
                if (match.obsIntensity > highestFragmentIntensity)
                    highestFragmentIntensity = match.obsIntensity;
                    //options.ConSole.WriteLine("fragment intensity higher than most intense fragment ... should not happen!");
            }
            AllProductMatches = cumulMatch;
            MatchingProductsFraction    = (double)MatchingWeightedProducts / (double) TotalWeightedProducts;
            MatchingIntensityFraction   = MatchingIntensity / (double)(highestFragmentIntensity * TotalTheoreticalProducts);
            if (MatchingIntensityFraction > 1)
                MatchingIntensityFraction = 1.0;
            ProductScore                = 1.0 - (cumulDiff / (double)(MatchingProducts * options.productMassTolerance.Value));
            IntensityScore              = MatchingIntensityFraction / (double)Query.spectrum.Peaks.Count;
            MorpheusScore               = MatchingProducts + MatchingIntensityFraction;
            FragmentScore               = ComputeFragmentScore();
        }
        
        public static int AscendingSpectrumNumberComparison(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return left.Query.spectrum.ScanNumber.CompareTo(right.Query.spectrum.ScanNumber);
        }

        public static int DescendingOptimizedScoreComparison(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            int comparison = -(left.ProbabilityScore().CompareTo(right.ProbabilityScore()));
            if (comparison != 0)
            {
                return comparison;
            }
            else
            {
                return left.Target.CompareTo(right.Target);
            }
        }

        public static readonly string Header = "Filename\tScan Number\tRetention Time (min)\tPrecursor m/z\tPrecursor Intensity\tPrecursor Charge\tPrecursor Mass (Da)\tExperimental Peaks\tTotal Intensity"
            + "\tPeptide Sequence\tBase Peptide Sequence\tProtein Description\tStart Residue Number\tStop Residue Number\tMissed Cleavages"
            + "\tTheoretical Mass (Da)\tPrecursor Mass Error (Da)\tPrecursor Mass Error (ppm)"
            + "\tMatching Products\tTotal Products\tRatio of Matching Products\tMatching Intensity\tFraction of Intensity Matching\tMorpheus Score";

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            //sb.Append(Query.spectrum.Filename + '\t');
            sb.Append(Query.spectrum.ScanNumber.ToString() + '\t');
            sb.Append(Query.spectrum.RetentionTimeInMin.ToString() + '\t');
            sb.Append(Query.spectrum.PrecursorMZ.ToString() + '\t');
            sb.Append(Query.spectrum.PrecursorIntensity.ToString() + '\t');
            sb.Append(Query.spectrum.PrecursorCharge.ToString() + '\t');
            sb.Append(Query.spectrum.PrecursorMass.ToString() + '\t');
            sb.Append(Query.spectrum.Peaks.Count.ToString() + '\t');
            sb.Append(Query.spectrum.TotalIntensity.ToString() + '\t');
            sb.Append(Peptide.ExtendedSequence() + '\t');
            sb.Append(Peptide.BaseSequence + '\t');
            sb.Append(Peptide.Parent.Description + '\t');
            sb.Append(Peptide.StartResidueNumber.ToString() + '\t');
            sb.Append(Peptide.EndResidueNumber.ToString() + '\t');
            sb.Append(Peptide.MissedCleavages.ToString() + '\t');
            sb.Append(Peptide.MonoisotopicMass.ToString() + '\t');
            sb.Append(PrecursorMzError.ToString() + '\t');
            sb.Append(MatchingProducts.ToString() + '\t');
            sb.Append(TotalTheoreticalProducts.ToString() + '\t');
            sb.Append(MatchingProductsFraction.ToString() + '\t');
            sb.Append(MatchingIntensity.ToString() + '\t');
            sb.Append(MatchingIntensityFraction.ToString() + '\t');
            sb.Append(ProbabilityScore().ToString());

            return sb.ToString();
        }
    }
}