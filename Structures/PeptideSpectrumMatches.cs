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
    public class PSMScoreThreshold
    {
        public double minProductScore = 0.0;
        public double minPrecursorScore = 0.0;
        public double minMatchingProductFractionScore = 0.0;
        public double minMatchingProductScore = 0.0;
        public double minIntensityFractionScore = 0.0;
        public double minIntensityScore = 0.0;
        public double minProteinScore = 0.0;
        public double minPeptideScore = 0.0;
        public double minFragmentScore = 0.0;
        public double minProbabilityScore = 0.0;
        public double ComputeFDR(PeptideSpectrumMatches matches, out long nbTarget, out long nbDecoy)
        {
            nbTarget = 0;
            nbDecoy = 0;
            foreach(PeptideSpectrumMatch psm in matches)
            {
                if (KeepPSM(psm))
                    if (psm.Target)
                        nbTarget++;
                    else
                        nbDecoy++;
            }
            if (nbTarget == 0)
                return 1;
            else
                return (double)nbDecoy / (double)nbTarget;
        }
        public bool KeepPSM(PeptideSpectrumMatch psm)
        {
            return psm.MatchingIntensity >= minIntensityScore &&
                   psm.MatchingIntensityFraction >= minIntensityFractionScore &&
                   psm.ProductScore >= minProductScore &&
                   psm.PrecursorScore >= minPrecursorScore &&
                   psm.MatchingProductsFraction >= minMatchingProductFractionScore &&
                   psm.MatchingProducts >= minMatchingProductScore &&
                   psm.ProteinScore >= minProteinScore &&
                   psm.PeptideScore >= minPeptideScore &&
                   psm.FragmentScore >= minFragmentScore &&
                   psm.ProbabilityScore() >= minProbabilityScore;
        }
    }

    public class PeptideSpectrumMatches : List<PeptideSpectrumMatch>
    {
        public PSMScoreThreshold ComputeScoreThreshold(double fdr)
        {
            PSMScoreThreshold localMinimas = new PSMScoreThreshold();
            double nbTarget = 0;
            double nbDecoy = 0;
            this.Sort(PeptideSpectrumMatches.CompareMatchingProducts);
            foreach (PeptideSpectrumMatch psm in this)
            {
                if (psm.Target) nbTarget += 1; else nbDecoy += 1;
                if (nbDecoy == 0 || nbDecoy / nbTarget <= fdr)
                    localMinimas.minMatchingProductScore = psm.MatchingProducts;
            }
            this.Sort(PeptideSpectrumMatches.CompareMatchingProductsFraction);
            nbTarget = 0;            nbDecoy = 0;
            foreach(PeptideSpectrumMatch psm in this)
            {
                if (psm.Target) nbTarget += 1; else nbDecoy += 1;
                if (nbDecoy == 0 || nbDecoy / nbTarget <= fdr)
                    localMinimas.minMatchingProductFractionScore = psm.MatchingProductsFraction;
            }
            this.Sort(PeptideSpectrumMatches.CompareMatchingIntensity);
            nbTarget = 0; nbDecoy = 0;
            foreach (PeptideSpectrumMatch psm in this)
            {
                if (psm.Target) nbTarget += 1; else nbDecoy += 1;
                if (nbDecoy == 0 || nbDecoy / nbTarget <= fdr)
                    localMinimas.minIntensityScore = psm.MatchingIntensity;
            }
            this.Sort(PeptideSpectrumMatches.CompareFragmentScore);
            nbTarget = 0; nbDecoy = 0;
            foreach (PeptideSpectrumMatch psm in this)
            {
                if (psm.Target) nbTarget += 1; else nbDecoy += 1;
                if (nbDecoy == 0 || nbDecoy / nbTarget <= fdr)
                    localMinimas.minFragmentScore = psm.FragmentScore;
            }
            this.Sort(PeptideSpectrumMatches.CompareProbabilityScore);
            nbTarget = 0; nbDecoy = 0;
            foreach (PeptideSpectrumMatch psm in this)
            {
                if (psm.Target) nbTarget += 1; else nbDecoy += 1;
                if (nbDecoy == 0 || nbDecoy / nbTarget <= fdr)
                    localMinimas.minProbabilityScore = psm.ProbabilityScore();
            }

            double iterSize = 0.95;
            long nbCumulTarget = 0;
            long nbCumulDecoy = 0;
            double currentFDR = localMinimas.ComputeFDR(this, out nbCumulTarget, out nbCumulDecoy);
            do
            {
                double previousVal = localMinimas.minMatchingProductScore;
                localMinimas.minMatchingProductScore *= iterSize;
                long nbTargetminMatchingProductScore = 0;
                double currentFDRminMatchingProductScore = localMinimas.ComputeFDR(this, out nbTargetminMatchingProductScore, out nbCumulDecoy);
                if (currentFDRminMatchingProductScore > fdr || localMinimas.minMatchingProductScore < 0.001)
                    nbTargetminMatchingProductScore = 0;
                localMinimas.minMatchingProductScore = previousVal;

                previousVal = localMinimas.minMatchingProductFractionScore;
                localMinimas.minMatchingProductFractionScore *= iterSize;
                long nbTargetminMatchingProductFractionScore = 0;
                double currentFDRminMatchingProductFractionScore = localMinimas.ComputeFDR(this, out nbTargetminMatchingProductFractionScore, out nbCumulDecoy);
                if (currentFDRminMatchingProductFractionScore > fdr || localMinimas.minMatchingProductFractionScore < 0.001)
                    nbTargetminMatchingProductFractionScore = 0;
                localMinimas.minMatchingProductFractionScore = previousVal;

                previousVal = localMinimas.minIntensityScore;
                localMinimas.minIntensityScore *= iterSize;
                long nbTargetminIntensityScore = 0;
                double currentFDRminIntensityScore = localMinimas.ComputeFDR(this, out nbTargetminIntensityScore, out nbCumulDecoy);
                if (currentFDRminIntensityScore > fdr || localMinimas.minIntensityScore < 0.001)
                    nbTargetminIntensityScore = 0;
                localMinimas.minIntensityScore = previousVal;

                previousVal = localMinimas.minFragmentScore;
                localMinimas.minFragmentScore *= iterSize;
                long nbTargetminFragmentScore = 0;
                double currentFDRminFragmentScore = localMinimas.ComputeFDR(this, out nbTargetminFragmentScore, out nbCumulDecoy);
                if (currentFDRminFragmentScore > fdr || localMinimas.minFragmentScore < 0.001)
                    nbTargetminFragmentScore = 0;
                localMinimas.minFragmentScore = previousVal;

                previousVal = localMinimas.minProbabilityScore;
                localMinimas.minProbabilityScore *= iterSize;
                long nbTargetminProbabilityScore = 0;
                double currentFDRminProbabilityScore = localMinimas.ComputeFDR(this, out nbTargetminProbabilityScore, out nbCumulDecoy);
                if (currentFDRminProbabilityScore > fdr || localMinimas.minProbabilityScore < 0.001)
                    nbTargetminProbabilityScore = 0;
                localMinimas.minProbabilityScore = previousVal;

                if (nbTargetminMatchingProductScore > 0 && nbTargetminMatchingProductScore >= nbTargetminMatchingProductFractionScore
                                                 && nbTargetminMatchingProductScore >= nbTargetminIntensityScore
                                                 && nbTargetminMatchingProductScore >= nbTargetminFragmentScore
                                                 && nbTargetminMatchingProductScore >= nbTargetminProbabilityScore)
                {
                    iterSize = 0.95;
                    localMinimas.minMatchingProductScore *= iterSize;
                }
                else
                    if (nbTargetminMatchingProductFractionScore > 0 && nbTargetminMatchingProductFractionScore >= nbTargetminMatchingProductScore
                                                                   && nbTargetminMatchingProductFractionScore >= nbTargetminIntensityScore
                                                                   && nbTargetminMatchingProductFractionScore >= nbTargetminFragmentScore
                                                                   && nbTargetminMatchingProductFractionScore >= nbTargetminProbabilityScore)
                    {
                        iterSize = 0.95;
                        localMinimas.minMatchingProductFractionScore *= iterSize;
                    }
                    else
                        if (nbTargetminIntensityScore > 0 && nbTargetminIntensityScore >= nbTargetminMatchingProductFractionScore
                                                         && nbTargetminIntensityScore >= nbTargetminMatchingProductScore
                                                         && nbTargetminIntensityScore >= nbTargetminFragmentScore
                                                         && nbTargetminIntensityScore >= nbTargetminProbabilityScore)
                        {
                            iterSize = 0.95;
                            localMinimas.minIntensityScore *= iterSize;
                        }
                        else
                            if (nbTargetminFragmentScore > 0 && nbTargetminFragmentScore >= nbTargetminMatchingProductFractionScore
                                                            && nbTargetminFragmentScore >= nbTargetminIntensityScore
                                                            && nbTargetminFragmentScore >= nbTargetminMatchingProductScore
                                                            && nbTargetminFragmentScore >= nbTargetminProbabilityScore)
                            {
                                iterSize = 0.95;
                                localMinimas.minFragmentScore *= iterSize;
                            }
                            else
                                if (nbTargetminProbabilityScore > 0 && nbTargetminProbabilityScore >= nbTargetminMatchingProductFractionScore
                                                                && nbTargetminProbabilityScore >= nbTargetminIntensityScore
                                                                && nbTargetminProbabilityScore >= nbTargetminMatchingProductScore
                                                                && nbTargetminProbabilityScore >= nbTargetminFragmentScore)
                                {
                                    iterSize = 0.95;
                                    localMinimas.minProbabilityScore *= iterSize;
                                }
                                else
                                    iterSize -= 0.05;

                currentFDR = localMinimas.ComputeFDR(this, out nbCumulTarget, out nbCumulDecoy);
            } while (iterSize > 0.01);// && currentFDR <= fdr);
            return localMinimas;
        }

        public PeptideSpectrumMatches() { }
        public PeptideSpectrumMatches(IEnumerable<PeptideSpectrumMatch> list) : base(list) { }
        
        public void OptimizePSMScoreRatios(DBOptions options, double desired_fdr, Result results)
        {
            long bestNbTargets = 0;
            double bestProtein = 0.1;
            double bestPeptide = 0.2;
            double bestFragments = 0.4;
            double bestIntensities = 0.1;
            double bestPrecursor = 0.2;

            double incr = 0.05;
            for (options.dProtein = 0.1; options.dProtein < 0.3; options.dProtein += incr)
                for (options.dPeptideScore = 0.1; options.dPeptideScore < 0.4; options.dPeptideScore += incr)
                    for (options.dPrecursor = 0.1; options.dPrecursor < 0.5; options.dPrecursor += incr)
                        for (options.dMatchingProductFraction = 0.1; options.dMatchingProductFraction < 1 - options.dPrecursor - options.dPeptideScore - options.dProtein; options.dMatchingProductFraction += incr)
                        {
                            double cumul = options.dPrecursor + options.dPeptideScore + options.dProtein + options.dMatchingProductFraction;
                            for (options.dIntensityFraction = 0.1; options.dIntensityFraction < 1 - cumul; options.dIntensityFraction += incr)
                            {
                                if (cumul + options.dIntensityFraction == 1)
                                {
                                    long nbTargets = 0;
                                    foreach (Precursor precursor in results.matchedPrecursors)
                                        if (precursor.Target)
                                            nbTargets++;
                                    if (nbTargets > bestNbTargets)
                                    {
                                        bestNbTargets = nbTargets;
                                        bestProtein = options.dProtein;
                                        bestPeptide = options.dPeptideScore;
                                        bestPrecursor = options.dPrecursor;
                                        bestIntensities = options.dIntensityFraction;
                                        bestFragments = options.dMatchingProductFraction;
                                    }
                                }
                            }

                        }

            options.dProtein = bestProtein;
            options.dPeptideScore = bestPeptide;
            options.dPrecursor = bestPrecursor;
            options.dIntensityFraction = bestIntensities;
            options.dMatchingProductFraction = bestFragments;

            options.ConSole.WriteLine("New score ratios   ----------------------------------------------------------------------- ");
            options.ConSole.WriteLine("    PeptideSpectrumMatch.dPrecursor:                     " + options.dPrecursor);
            options.ConSole.WriteLine("    PeptideSpectrumMatch.dMatchingProductFraction:       " + options.dMatchingProductFraction);
            options.ConSole.WriteLine("    PeptideSpectrumMatch.dIntensityFraction:             " + options.dIntensityFraction);
            options.ConSole.WriteLine("    PeptideSpectrumMatch.dPeptideScore:                  " + options.dPeptideScore);
            options.ConSole.WriteLine("    PeptideSpectrumMatch.dProtein:                       " + options.dProtein);
            options.ConSole.WriteLine("------------------------------------------------------------------------------------------ ");
        }
        
        public double ComputePrecursorArea(bool smooth)
        {
            double lastTimeStamp = 0;
            double lastIntensityPerUnitOfTime = 0;
            List<double> timeGap = new List<double>();
            List<double> precursorIntensities = new List<double>();
            this.Sort(PeptideSpectrumMatches.AscendingRetentionTime);

            foreach (PeptideSpectrumMatch psm in this)
            {
                if (psm.Query.spectrum.Ms1InjectionTime > 0)
                {
                    if (psm.Query.spectrum.PrecursorIntensity > 0 && lastTimeStamp > 0)
                    {
                        timeGap.Add((psm.Query.spectrum.RetentionTimeInMin - lastTimeStamp) * 60.0 * 1000.0);
                        precursorIntensities.Add(psm.Query.spectrum.PrecursorIntensityPerMilliSecond);
                    }
                    lastIntensityPerUnitOfTime = psm.Query.spectrum.PrecursorIntensityPerMilliSecond;
                    lastTimeStamp = psm.Query.spectrum.RetentionTimeInMin;
                }
            }
            //Smooth the curve

            if (smooth)
            {
                List<double> timeGapSmooth = new List<double>();
                List<double> precursorIntensitiesSmoothed = new List<double>();
                for (int i = 0; i < timeGap.Count; i++)
                {
                    double cumulIntensity = 0;
                    int nbItems = 0;
                    for (int k = (i - 2 < 0 ? 0 : i - 2); k < timeGap.Count - 2; k++)
                    {
                        nbItems++;
                        cumulIntensity += precursorIntensities[k];
                    }
                    if (nbItems > 0)
                    {
                        timeGapSmooth.Add(timeGap[0]);
                        precursorIntensitiesSmoothed.Add(cumulIntensity / (double)nbItems);
                    }
                }

                double fragSpectrumArea = 0;
                for (int i = 0; i < timeGapSmooth.Count; i++)
                    fragSpectrumArea += timeGapSmooth[i] * precursorIntensitiesSmoothed[i];
                return fragSpectrumArea;
            }
            else
            {
                double fragSpectrumArea = 0;
                for (int i = 0; i < timeGap.Count; i++)
                    fragSpectrumArea += timeGap[i] * precursorIntensities[i];
                return fragSpectrumArea;
            }
        }

        public double ComputePrecursorAreaBKP(bool smooth)
        {
            double lastTimeStamp = 0;
            double lastIntensityPerUnitOfTime = 0;
            List<double> timeGap = new List<double>();
            List<double> precursorIntensities = new List<double>();
            foreach (PeptideSpectrumMatch psm in this)
            {
                if (psm.Query.spectrum.Ms1InjectionTime > 0)
                {
                    if (psm.Query.spectrum.PrecursorIntensity > 0 && lastTimeStamp > 0)
                    {
                        timeGap.Add((psm.Query.spectrum.RetentionTimeInMin - lastTimeStamp) * 60.0 * 1000.0);
                        if (lastIntensityPerUnitOfTime > 0)
                            precursorIntensities.Add((psm.Query.spectrum.PrecursorIntensity / psm.Query.spectrum.Ms1InjectionTime + lastIntensityPerUnitOfTime) * 0.5);
                        else
                            precursorIntensities.Add(psm.Query.spectrum.PrecursorIntensity / psm.Query.spectrum.Ms1InjectionTime);
                    }
                    lastIntensityPerUnitOfTime = psm.Query.spectrum.PrecursorIntensity / psm.Query.spectrum.Ms1InjectionTime;
                    lastTimeStamp = psm.Query.spectrum.RetentionTimeInMin;
                }
            }
            //Smooth the curve

            if (smooth)
            {
                List<double> timeGapSmooth = new List<double>();
                List<double> precursorIntensitiesSmoothed = new List<double>();
                for (int i = 0; i < timeGap.Count; i++)
                {
                    double cumulIntensity = 0;
                    int nbItems = 0;
                    for (int k = (i - 2 < 0 ? 0 : i - 2); k < timeGap.Count - 2; k++)
                    {
                        nbItems++;
                        cumulIntensity += precursorIntensities[k];
                    }
                    if (nbItems > 0)
                    {
                        timeGapSmooth.Add(timeGap[0]);
                        precursorIntensitiesSmoothed.Add(cumulIntensity / (double)nbItems);
                    }
                }

                double fragSpectrumArea = 0;
                for (int i = 0; i < timeGapSmooth.Count; i++)
                    fragSpectrumArea += timeGapSmooth[i] * precursorIntensitiesSmoothed[i];
                return fragSpectrumArea;
            }
            else
            {
                double fragSpectrumArea = 0;
                for (int i = 0; i < timeGap.Count; i++)
                    fragSpectrumArea += timeGap[i] * precursorIntensities[i];
                return fragSpectrumArea;
            }
        }

        public Dictionary<PeptideSpectrumMatch, double> ComputeMsMsNormalizationFactors()
        {
            Dictionary<PeptideSpectrumMatch, double> fragRatio = new Dictionary<PeptideSpectrumMatch, double>();
            double lastIntensity = 0;
            foreach (PeptideSpectrumMatch psm in this)
            {
                if (psm.Query.spectrum.PrecursorIntensity > 0)
                {
                    double intensityFactor = 0;
                    if (psm.Query.spectrum.InjectionTime >= 119.999997317791)//TODO Find instrument/method specific default or max injection time
                        lastIntensity = psm.Query.spectrum.PrecursorIntensity;
                    else
                    {
                        double predictedIntensity = (lastIntensity);// + psm.Query.spectrum.PrecursorIntensity) * 0.5;
                        intensityFactor = (psm.Query.spectrum.PrecursorIntensity - predictedIntensity) / predictedIntensity;// psm.Query.spectrum.PrecursorIntensity;
                    }
                    fragRatio.Add(psm, intensityFactor);
                }
            }
            return fragRatio;
        }    

        public List<ProductMatch> GetCombinedSpectrum(DBOptions dbOptions, Utilities.ElutionCurve eCurvePerMs,  Dictionary<double, int> DicOfCommonPM, ref List<double> unstableMasses)
        {
            Dictionary<PeptideSpectrumMatch, double> DicOfPsmFactor = new Dictionary<PeptideSpectrumMatch, double>();// this.ComputeMsMsNormalizationFactors();
            //Dictionary<ProductMatch, double> DicOfProductMsMsFactor = new Dictionary<ProductMatch, double>();
            Dictionary<string, Dictionary<PeptideSpectrumMatch, ProductMatch>> DicOfProducts = new Dictionary<string, Dictionary<PeptideSpectrumMatch, ProductMatch>>();

            List<double> avgProbList = new List<double>();            
            foreach (PeptideSpectrumMatch psm in this)
            {
                bool usedPsm = false;
                foreach (ProductMatch match in psm.AllProductMatches)
                {
                    if (DicOfCommonPM == null || DicOfCommonPM.ContainsKey(match.theoMz))
                    {
                        string key = match.Fragment.Name + "|" + match.fragmentPos + "|" + match.charge;
                        if (!DicOfProducts.ContainsKey(key))
                            DicOfProducts.Add(key, new Dictionary<PeptideSpectrumMatch, ProductMatch>());
                        DicOfProducts[key].Add(psm, match);
                        //DicOfProductMsMsFactor.Add(match, DicOfPsmFactor[psm]);
                        usedPsm = true;
                    }
                }
                if (usedPsm)
                {
                    DicOfPsmFactor.Add(psm, psm.ProbabilityScore());//.MatchingIntensity);//psm.Query.precursor.Track.INTENSITY);
                    avgProbList.Add(psm.ProbabilityScore());
                }
            }
            double avgProbability = MathNet.Numerics.Statistics.ArrayStatistics.Mean(avgProbList.ToArray());
            avgProbability -= avgProbability * 0.5;
                //MathNet.Numerics.Statistics.ArrayStatistics.Variance(avgProbList.ToArray());

            double avgNormedInt = 0;
            List<ProductMatch> products = new List<ProductMatch>();
            //if (DicOfProductMsMsFactor.Count > 0)
            {
                foreach (Dictionary<PeptideSpectrumMatch, ProductMatch> matchList in DicOfProducts.Values)
                {
                    List<ProductMatch> newPMList = new List<ProductMatch>();
                    if (matchList.Count > 0)
                    {
                        double sumPsmFactor = 0;
                        foreach (PeptideSpectrumMatch psm in matchList.Keys)
                        {
                            if (psm.ProbabilityScore() > avgProbability)//Keep only above average spectrum
                            {
                                ProductMatch pm = matchList[psm];
                                //if (newPM == null)
                                //{
                                //    newPM = new ProductMatch(pm);
                                //    newPM.obsIntensity = 0;
                                //    newPM.normalizedIntensity = 0;
                                //}
                                ProductMatch newPM = new ProductMatch(pm);
                                newPM.obsIntensity = pm.obsIntensity * DicOfPsmFactor[psm];// +pm.obsIntensity * DicOfProductMsMsFactor[pm];
                                //newPM.obsIntensity += pm.obsIntensity + pm.obsIntensity * DicOfProductMsMsFactor[pm];
                                newPM.normalizedIntensity = pm.normalizedIntensity * DicOfPsmFactor[psm];
                                newPMList.Add(newPM);
                                sumPsmFactor += DicOfPsmFactor[psm];
                            }
                        }
                        if (sumPsmFactor > 0 && newPMList.Count > 4)
                        {
                            double[] nIArray = new double[newPMList.Count];
                            for (int i = 0; i < nIArray.Length; i++)
                                nIArray[i] = newPMList[i].normalizedIntensity;
                            double variance = Numerics.Variance(nIArray);
                            double stDev = MathNet.Numerics.Statistics.ArrayStatistics.StandardDeviation(nIArray);
                            double mean = MathNet.Numerics.Statistics.ArrayStatistics.Mean(nIArray);
                            if (mean - stDev > mean * 0.5 || variance < mean * 0.5)
                            {
                                ProductMatch newPM = new ProductMatch(newPMList[0]);
                                newPM.weight = 0.0;
                                newPM.obsIntensity = 0.0;
                                newPM.normalizedIntensity = 0.0;
                                for (int i = 0; i < newPMList.Count; i++)
                                {
                                    newPM.obsIntensity += newPMList[i].obsIntensity;
                                    newPM.normalizedIntensity += newPMList[i].normalizedIntensity;
                                    if (newPMList[i].normalizedIntensity > 0.0)
                                        newPM.weight++;
                                }
                                newPM.obsIntensity /= sumPsmFactor;
                                newPM.normalizedIntensity /= sumPsmFactor;

                                //newPM.weight = matchList.Count * newPM.normalizedIntensity;
                                newPM.weight *= newPM.normalizedIntensity;
                                avgNormedInt += newPM.normalizedIntensity;
                                products.Add(newPM);
                            }
                            else
                                unstableMasses.Add(newPMList[0].theoMz);
                        }
                    }
                }
                if(products.Count > 0)
                    avgNormedInt /= (double)products.Count;

                //Add missed important fragments
                if (DicOfCommonPM != null)
                {
                    foreach (double mz in DicOfCommonPM.Keys)
                    {
                        bool found = false;
                        foreach (ProductMatch match in products)
                            if (match.theoMz == mz)
                                found = true;
                        if (!found)
                        {
                            double sumPsmFactor = 0;
                            ProductMatch newMatch = new ProductMatch();
                            newMatch.theoMz = mz;
                            newMatch.weight = 0;
                            newMatch.obsIntensity = 0;
                            newMatch.normalizedIntensity = 0;
                            foreach (PeptideSpectrumMatch psm in DicOfPsmFactor.Keys)
                            {
                                double timePoint = psm.Query.spectrum.RetentionTimeInMin * 1000.0 * 60.0;
                                double intensityInCTrap = eCurvePerMs.GetLocalArea(timePoint, timePoint + psm.Query.spectrum.InjectionTime);
                                double peakIntensity = 0.0;
                                foreach (MsMsPeak peak in psm.Query.spectrum.Peaks)
                                {
                                    if (peak.Intensity > 0 && Math.Abs(Utilities.Numerics.CalculateMassError(peak.MZ, mz, dbOptions.productMassTolerance.Units)) <= dbOptions.productMassTolerance.Value)
                                        peakIntensity += peak.Intensity;
                                }
                                if (peakIntensity > 0)
                                {
                                    newMatch.obsIntensity += peakIntensity * DicOfPsmFactor[psm];
                                    //MixedPrecursor.GetLocalArea
                                    newMatch.normalizedIntensity += (peakIntensity / intensityInCTrap) * DicOfPsmFactor[psm];
                                    sumPsmFactor += DicOfPsmFactor[psm];
                                    newMatch.weight++;
                                }
                            }
                            if (newMatch.weight > 0)
                            {
                                newMatch.obsIntensity /= sumPsmFactor;
                                newMatch.normalizedIntensity /= sumPsmFactor;
                            }
                            newMatch.weight *= newMatch.normalizedIntensity;
                            products.Add(newMatch);
                        }
                    }
                }
                
                //Keep only most intense fragments (5% of average normalized intensity)
                foreach(ProductMatch pm in products)
                    if (pm.normalizedIntensity < avgNormedInt * 0.1)//0.05
                    {
                        pm.normalizedIntensity = 0;
                        pm.obsIntensity = 0;
                    }
            }

            return products;
        }

        public static int CompareMatchingIntensityFraction(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.MatchingIntensityFraction.CompareTo(right.MatchingIntensityFraction);
        }
        public static int CompareTrapDistance(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.Query.TrapDistance.CompareTo(right.Query.TrapDistance);
        }
        public static int AscendingRetentionTime(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return left.Query.spectrum.RetentionTimeInMin.CompareTo(right.Query.spectrum.RetentionTimeInMin);
        }
        public static int CompareMatchingIntensity(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.MatchingIntensity.CompareTo(right.MatchingIntensity);
        }
        public static int CompareProductScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.ProductScore.CompareTo(right.ProductScore);
        }
        public static int CompareProbabilityScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.ProbabilityScore().CompareTo(right.ProbabilityScore());
        }
        public static int ComparePrecursorScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.PrecursorScore.CompareTo(right.PrecursorScore);
        }
        public static int CompareMatchingProductsFraction(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.MatchingProductsFraction.CompareTo(right.MatchingProductsFraction);
        }
        public static int CompareFragmentScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.FragmentScore.CompareTo(right.FragmentScore);
        }
        public static int CompareMatchingProducts(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.MatchingWeightedProducts.CompareTo(right.MatchingWeightedProducts);
        }
        public static int CompareProteinScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.ProteinScore.CompareTo(right.ProteinScore);
        }
        public static int ComparePeptideScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.PeptideScore.CompareTo(right.PeptideScore);
        }
    }
}
