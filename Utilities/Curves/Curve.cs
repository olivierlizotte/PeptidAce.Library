using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace PeptidAce.Utilities
{
    public enum CurveType { LINEAR, AKIMA, SPLINE }

    public class ElutionCurve
    {
        public double Area = 0.0;
        public double[] Coefficients = null;        

        public List<double> time = null;
        public List<double> intensityCount = null;

        public List<double> interpolatedTime = null;
        public List<double> interpolatedIntensityCount = null;

        private Dictionary<double, double> dicNoDupe = null;

        private CurveType cType = CurveType.AKIMA;


        public static List<double> GetTimePoints(int nbTimePoints, bool keepOnlyAboveZeroIntensities, List<double> time, List<double> intensityCount = null)
        {
            int minIndex = 0;
            int maxIndex = time.Count - 1;
            if (keepOnlyAboveZeroIntensities && intensityCount != null)
            {
                double highestIntensity = 0;
                int indexHighestIntensity = 0;
                for (int i = 0; i < time.Count; i++)
                    if (intensityCount[i] > highestIntensity)
                    {
                        highestIntensity = intensityCount[i];
                        indexHighestIntensity = i;
                    }

                minIndex = indexHighestIntensity;
                while (minIndex > 0 && !(intensityCount[minIndex] == 0 && intensityCount[minIndex - 1] == 0))
                    minIndex--;

                maxIndex = indexHighestIntensity;
                while (maxIndex < time.Count - 1 && !(intensityCount[maxIndex] == 0 && intensityCount[maxIndex + 1] == 0))
                    maxIndex++;
            }
            //if (maxIndex - minIndex > nbTimePoints)
            {
                double minTime = time[minIndex];
                double maxTime = time[maxIndex];
                List<double> points = new List<double>(nbTimePoints);
                for (int i = 0; i < nbTimePoints; i++)
                    points.Add(minTime + ((i + 0.5) / (double)nbTimePoints) * (maxTime - minTime));
                return points;
            }
            //else
            //    return time;
        }

        public int GetNbPoints()
        {
            return intensityCount.Count;
        }

        private MathNet.Numerics.Interpolation.IInterpolation interpole = null;
        
        public static ElutionCurve Create(Dictionary<double, double> dicOfTimeInMsVsIntensityPerMs)
        {
            ElutionCurve theCurve = new ElutionCurve();
            // -- Test curve fitting function -- //
            theCurve.time = new List<double>();
            theCurve.intensityCount = new List<double>();
            List<double> sortedKeys = new List<double>(dicOfTimeInMsVsIntensityPerMs.Keys);
            sortedKeys.Sort();
            foreach (double time in sortedKeys)
            {
                theCurve.time.Add(time);
                theCurve.intensityCount.Add(dicOfTimeInMsVsIntensityPerMs[time]);
            }
            theCurve.Compute(CurveType.LINEAR, true);
            return theCurve;
        }

        public double InterpolateIntensity(double timePoint)
        {
            if (interpole != null)
            {
                double intensity = interpole.Interpolate(timePoint);
                if (intensity < 0 || double.IsNaN(intensity))
                    return 0;
                else
                    return intensity;
            }
            return 0;
        }

        /// <summary>
        /// Computes the curve and area (filled using AddPoint routine)
        /// </summary>
        /// <param name="useInterpole"></param>
        public void Compute(CurveType type, bool smooth)
        {
            cType = type;
            if (dicNoDupe != null && dicNoDupe.Count > 8)// time != null && time.Count > 8)
            {
                time = new List<double>(dicNoDupe.Keys);
                intensityCount = new List<double>(dicNoDupe.Values);
            }

            if(time != null && (time.Count > 8 || (time.Count >=4 && type == CurveType.LINEAR)))
            {                    
                double[] arrayTime = time.ToArray();
                double[] arrayIntensity = intensityCount.ToArray();
                Array.Sort(arrayTime, arrayIntensity);
                time = new List<double>(arrayTime);
                intensityCount = new List<double>(arrayIntensity);

                double area1 = Utilities.CurveFitter.FitToPolynomial(time.ToArray(), intensityCount.ToArray(), out Coefficients);
                double area2 = Utilities.CurveFitter.AreaUnderTheCurve(time, intensityCount);
                Area = area2;

                try
                {
                    switch(cType)
                    {
                        case CurveType.AKIMA:
                            interpole = MathNet.Numerics.Interpolation.LinearSpline.Interpolate(time, intensityCount);
                            //interpole = new MathNet.Numerics.Interpolation.NevillePolynomialInterpolation(time, intensityCount);
                            //interpole = new MathNet.Numerics.Interpolation.Algorithms.AkimaSplineInterpolation(time, intensityCount);
                            break;
                        case CurveType.LINEAR:
                            interpole = MathNet.Numerics.Interpolation.LinearSpline.Interpolate(time, intensityCount);
                            break;
                        case CurveType.SPLINE:
                            interpole = new MathNet.Numerics.Interpolation.BulirschStoerRationalInterpolation(time, intensityCount);
                            break;
                    }
                    //Area = interpole.Integrate(time[time.Count - 1]);
                    //int nbPoints = time.Count / 2;
                    //if (nbPoints)
                    //    nbPoints = time.Count + 2;
                    if (smooth)
                    {
                        int nbIter = 2;
                        interpolatedTime = time;
                        interpolatedIntensityCount = intensityCount;
                        do
                        {
                            List<double> smallerInterpolatedTime = new List<double>(interpolatedTime.Count - 1);
                            List<double> smallerInterpolatedIntensityCount = new List<double>(interpolatedIntensityCount.Count - 1);
                            for (int i = 0; i < interpolatedTime.Count - 1; i++)
                            {
                                smallerInterpolatedTime.Add(0.5 * (interpolatedTime[i] + interpolatedTime[i + 1]));
                                smallerInterpolatedIntensityCount.Add(0.5 * (interpolatedIntensityCount[i] + interpolatedIntensityCount[i + 1]));
                            }
                            interpolatedTime = smallerInterpolatedTime;
                            interpolatedIntensityCount = smallerInterpolatedIntensityCount;
                            nbIter--;
                        } while (nbIter > 0);

                        switch (cType)
                        {
                            case CurveType.AKIMA:
                                interpole = MathNet.Numerics.Interpolation.LinearSpline.Interpolate(time, intensityCount);
                                //interpole = new MathNet.Numerics.Interpolation.NevillePolynomialInterpolation(time, intensityCount);
                                //interpole = new MathNet.Numerics.Interpolation.Algorithms.AkimaSplineInterpolation(interpolatedTime, interpolatedIntensityCount);
                                break;
                            case CurveType.LINEAR:
                                interpole = MathNet.Numerics.Interpolation.LinearSpline.Interpolate(time, intensityCount);
                                //interpole = MathNet.Numerics.Interpolation.Interpolate.LinearBetweenPoints(interpolatedTime, interpolatedIntensityCount);
                                break;
                            case CurveType.SPLINE:
                                interpole = new MathNet.Numerics.Interpolation.BulirschStoerRationalInterpolation(time, intensityCount);
                                //interpole = new MathNet.Numerics.Interpolation.Algorithms.CubicSplineInterpolation(interpolatedTime, interpolatedIntensityCount);
                                break;
                        }
                        double areaSmooth = interpole.Integrate(interpolatedTime[interpolatedTime.Count - 1]);
                        if (areaSmooth > Area * 0.5 && areaSmooth < Area * 2)
                            Area = areaSmooth;
                        else
                            Console.WriteLine("Contradiction");
                    }
                    //Area = interpole.Integrate(interpolatedTime[interpolatedTime.Length - 1]);
                }
                catch (Exception ex)
                {
                    Console.WriteLine(ex.StackTrace);
                }

            }
            else
            {
                Area = 0;
                Coefficients = new double[0];
            }
        }

        public double GetLocalArea(double timeStart, double timeStop)
        {
            if (interpole != null)
            {
                double cumul = 0.0;
                double stepSize = (timeStop - timeStart) / 100.0;
                for (double k = timeStart; k <= timeStop; k += stepSize)
                    cumul += interpole.Interpolate(k) * stepSize;
                return cumul;
            }
            else
                return 0;
        }

        /// <summary>
        /// Add a point to a curve that is still under construction
        /// </summary>
        /// <param name="newTimePoint"></param>
        /// <param name="newIntensityPerMilliSeconds"></param>
        public void AddPoint(double newTimePoint, double newIntensityPerMilliSeconds)
        {
            if (dicNoDupe == null)
            {
                dicNoDupe = new Dictionary<double, double>();
                dicNoDupe.Add(newTimePoint, newIntensityPerMilliSeconds);
            }
            else
            {
                if (dicNoDupe.ContainsKey(newTimePoint))
                    dicNoDupe[newTimePoint] = 0.5 * (dicNoDupe[newTimePoint] + newIntensityPerMilliSeconds);
                else
                    dicNoDupe.Add(newTimePoint, newIntensityPerMilliSeconds);
             }
        }
    }
}
