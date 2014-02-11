using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace PeptidAce.Utilities
{
    /// <summary>
    /// This structure holds curves fot Peptide intensity count and intensity per milliseconds
    /// </summary>
    public class MaxFlowElutionCurve
    {
        public int nbProducts;
        public ElutionCurve eCurveCount;
        public ElutionCurve eCurvePerMs;
        public MaxFlowElutionCurve(int nbProductsUsed)
        {
            this.nbProducts = nbProductsUsed;
            eCurveCount = new ElutionCurve();
            eCurvePerMs = new ElutionCurve();
        }

        public void Compute()
        {
            eCurveCount.Compute(CurveType.LINEAR, true);
            eCurvePerMs.Compute(CurveType.LINEAR, true);
        }
    }

    /// <summary>
    /// Merges a bunch of ElutionCurves that represent the same peptide in order to gain precision
    /// </summary>
    public class ElutionCurveMerger
    {
        //Curves to merge
        private List<MaxFlowElutionCurve> Curves = new List<MaxFlowElutionCurve>();

        //Factor (or wight) that each curve should have
        private List<double> Factor = new List<double>();

        /// <summary>
        /// Add a curve to be merged
        /// </summary>
        /// <param name="newCurve"></param>
        /// <param name="weight"></param>
        public void AddCurve(MaxFlowElutionCurve newCurve, double weight)
        {
            Curves.Add(newCurve);
            Factor.Add(weight);
        }

        /// <summary>
        /// Merges curves added throught the AddCurve function. The result is a more precise curve
        /// </summary>
        /// <returns></returns>
        public MaxFlowElutionCurve Merge()
        {
            if (Curves.Count > 1)
            {
                double sum = 0.0;
                foreach (double val in Factor)
                    sum += val;

                Dictionary<double, int> times = new Dictionary<double, int>();
                foreach (MaxFlowElutionCurve curve in Curves)
                    foreach (double timePoint in curve.eCurvePerMs.time)
                        if (!times.ContainsKey(timePoint))
                            times.Add(timePoint, 1);
                        else
                            times[timePoint]++;
                List<double> sortedTime = new List<double>(times.Keys);
                sortedTime.Sort();

                List<double> interpolTime = ElutionCurve.GetTimePoints(128, false, sortedTime);

                MaxFlowElutionCurve newCurve = new MaxFlowElutionCurve(-1);
                foreach (double pt in interpolTime)
                {
                    double intCount = 0;
                    double intPerMs = 0;
                    for (int i = 0; i < Curves.Count; i++)
                    {
                        intCount += Factor[i] * Curves[i].eCurveCount.InterpolateIntensity(pt);
                        intPerMs += Factor[i] * Curves[i].eCurvePerMs.InterpolateIntensity(pt);
                    }
                    newCurve.eCurveCount.AddPoint(pt, intCount / sum);
                    newCurve.eCurvePerMs.AddPoint(pt, intPerMs / sum);
                }
                newCurve.Compute();
                return newCurve;

                /*
                Dictionary<double, int> times = new Dictionary<double, int>();
                foreach (MaxFlowElutionCurve curve in Curves)
                    foreach (double timePoint in curve.eCurvePerMs.time)
                        if (!times.ContainsKey(timePoint))
                            times.Add(timePoint, 1);
                        else
                            times[timePoint]++;
                List<double> sortedTime = new List<double>(times.Keys);
                sortedTime.Sort();

                MaxFlowElutionCurve newCurve = new MaxFlowElutionCurve(-1);
                foreach (double key in sortedTime)
                    if (times[key] > 1)
                    {
                        double cumulIntensityCount = 0.0;
                        double cumulIntensityPerMs = 0.0;
                        for (int i = 0; i < Curves.Count; i++)
                        {
                            cumulIntensityCount += Curves[i].eCurveCount.InterpolateIntensity(key) * Factor[i] / sum;
                            cumulIntensityPerMs += Curves[i].eCurvePerMs.InterpolateIntensity(key) * Factor[i] / sum;
                        }
                        newCurve.eCurveCount.AddPoint(key, cumulIntensityCount);
                        newCurve.eCurvePerMs.AddPoint(key, cumulIntensityPerMs);
                    }
                return newCurve;//*/
            }
            else if (Curves.Count == 1)
                return Curves[0];
            return new MaxFlowElutionCurve(-1);
        }
    }

}
