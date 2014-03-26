using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using PeptidAce.Utilities.Interfaces;

namespace PeptidAce.Utilities.Methods
{
    public static class GradientDescent
    {
        public static double ComputeOver(Dictionary<double, double> virtualMixed, Dictionary<double, double> mixed)        
        {
            double cumulOver = 0;
            foreach (double key in virtualMixed.Keys)
                if (virtualMixed[key] > mixed[key])
                    cumulOver += virtualMixed[key] - mixed[key];
            return cumulOver;
        }

        public static double ComputeUnder(Dictionary<double, double> virtualMixed, Dictionary<double, double> mixed)
        {
            double cumulUnder = 0;
            foreach (double key in virtualMixed.Keys)
                if (mixed[key] > virtualMixed[key])
                    cumulUnder += mixed[key] - virtualMixed[key];
            return cumulUnder;
        }
        
        public static void SolveMinimaStyle(List<Dictionary<double, double>> units, Dictionary<double, double> mixed,
                                 out List<double> solution, out double underflow, IConSol ConSole)
        {
            List<double> localFlows = new List<double>();
            foreach (Dictionary<double, double> unit in units)
                localFlows.Add(FindLocalMaxima(unit, mixed));

            Dictionary<double, double> virtualMixed = BuildVirtualDic(localFlows, units, mixed.Count);
            double overError = ComputeOver(virtualMixed, mixed);
            double underError = ComputeUnder(virtualMixed, mixed);

            int bestUnit = 0;
            while (overError >= 1 && bestUnit >= 0)
            {
                //double bestFlow = 0;
                double bestMinima = 0;
                bestUnit = -1;
                for (int i = 0; i < units.Count; i++)
                {
                    if (localFlows[i] > 0)
                    {
                        double minima = FindLocalMinima(units[i], mixed, virtualMixed);
                        //double currentFlow = localFlows[i];
                        //localFlows[i] = minima;
                        //Dictionary<double, double> tmpDic = BuildVirtualDic(localFlows, units);
                        //double tmpUnderError = ComputeUnder(tmpDic, mixed);
                        //double tmpOverError  = ComputeOver(tmpDic, mixed);

                        //double tmpFlowRate = Math.Abs(overError - tmpOverError);
                        //if (tmpUnderError > underError)
                        //    tmpFlowRate /= tmpUnderError - underError;
                        //if (tmpFlowRate > bestFlow)
                        if(minima > bestMinima)
                        {
                            //bestFlow = tmpFlowRate;
                            bestMinima = minima;
                            bestUnit = i;
                        }
                        //localFlows[i] = currentFlow;
                    }
                }
                if (bestUnit >= 0)
                {
                    if (bestMinima > 1)
                        localFlows[bestUnit] -= 1.0;// *0.01
                    else 
                        localFlows[bestUnit] -= bestMinima;// *0.01;
                    if (localFlows[bestUnit] < 0)
                        localFlows[bestUnit] = 0.0;

                    virtualMixed = BuildVirtualDic(localFlows, units, mixed.Count);
                    overError = ComputeOver(virtualMixed, mixed);
                    underError = ComputeUnder(virtualMixed, mixed);
                }
            }//End of while overflow > 1

            solution = new List<double>();
            foreach (double localFlow in localFlows)
                if (overError <= 1.0)
                    solution.Add(localFlow);
                else
                    solution.Add(0);

            underflow = underError;
        }//*/
        
        public class GDHelper
        {
            private List<Dictionary<double, double>> units;
            private Dictionary<double, double> mixed;
            public double sumIntensities;
            double PrecursorIntensityInCTrap;
            public GDHelper(List<Dictionary<double, double>> pUnits, Dictionary<double, double> pMixed, double precursorIntensityInCTrap)
            {
                units = pUnits;            
                mixed = pMixed;
                sumIntensities = 0.0;
                foreach (double val in pMixed.Values)
                    sumIntensities += val;
                PrecursorIntensityInCTrap = precursorIntensityInCTrap;
            }

            public double FCT(List<double> ratios)
            {                
                Dictionary<double, double> virtualMixed = BuildVirtualDic(ratios, units, mixed.Count);
                double tmpErrorOver = ComputeOver(virtualMixed, mixed);// / sumIntensities;                
                double tmpErrorUnder = ComputeUnder(virtualMixed, mixed) / sumIntensities;
                double belowZero = 0;
                for (int i = 0; i < ratios.Count; i++)
                {
                    if (ratios[i] < 0)
                        foreach(double val in units[i].Values)
                            belowZero += val * ratios[i];
                }
                belowZero = -belowZero;// / sumIntensities;//*/
                //Scoring goal: Result must have OverFlow < 1; ratios below zero are not permitted; Minimize Underflow
                //score < 1 => only underflow
                //score < 2 => some overflow
                //score < 3 => ratios under zero
                        //belowZero += val;

                
                //return 100.0 * belowZero + 
                //   return 100.0 * tmpErrorOver + tmpErrorUnder;

                return belowZero + tmpErrorOver + tmpErrorUnder;
                /*
                if(belowZero > 0)
                    return 2 + belowZero + tmpErrorOver + tmpErrorUnder;
                if(tmpErrorOver > 0)
                    return 1 + tmpErrorOver + tmpErrorUnder;
                return tmpErrorUnder;//*/

                //return tmpErrorOver + 
                //return PrecursorIntensityInCTrap * (tmpErrorOver - belowZero) + tmpErrorUnder;//0.1675
                //return 100*tmpErrorOver + tmpErrorUnder - belowZero;//0.1675
                //return 100 * tmpErrorOver + tmpErrorUnder - belowZero * 1000;
                //if(tmpErrorUnder > 0)
                //    return tmpErrorOver / tmpErrorUnder;
                //else
                //    return tmpErrorOver;
            }
        }

        public delegate double Function(List<double> x);
        public static void SolveFromGradientDescent(List<Dictionary<double, double>> units, Dictionary<double, double> mixed,
                                double PrecursorIntensityInCTrap,
                                out List<double> solution, out double underflow, IConSol ConSole)
        {
            List<double> seeds = new List<double>();
            foreach (Dictionary<double, double> unit in units)
                seeds.Add(FindLocalMaxima(unit, mixed));

            GDHelper gdh = new GDHelper(units, mixed, PrecursorIntensityInCTrap);
            List<double> localMinima;
            if (seeds.Count == 1)
            {
                localMinima = new List<double>();
                localMinima.Add(FindLocalMaxima(units[0], mixed));
            }
            else
                localMinima = Gradior.Minimize(gdh.FCT, seeds, 0.01, 1, 1000, 0.05);//(long)PrecursorIntensityInCTrap * units.Count, 0.05);
            solution = new List<double>();
            foreach (double localFlow in localMinima)
                solution.Add(localFlow);

            underflow = ComputeUnder(BuildVirtualDic(localMinima, units, mixed.Count), mixed);
        }

        public static void SolveMaxFlowStyle(List<Dictionary<double, double>> units, Dictionary<double, double> mixed,
                                 out List<double> solution, out double underflow, IConSol ConSole, double stepSize)
        {
            List<double> localFlows = new List<double>();
            foreach (Dictionary<double, double> unit in units)
                localFlows.Add(FindLocalMaxima(unit, mixed));

            Dictionary<double, double> virtualMixed = BuildVirtualDic(localFlows, units, mixed.Count);
            double overError = ComputeOver(virtualMixed, mixed);
            double underError = ComputeUnder(virtualMixed, mixed);
            double[] bestIndexes = new double[units.Count];

            double iterSize = 1;
            double bestOverallError = double.MaxValue;
            List<double> bestLocalFlows = new List<double>();
            Random rnd = new Random();
            while (overError >= 1 && iterSize < 10000)//anything less than 1 is an acceptable solution
            {
                for (int index = 0; index < bestIndexes.Length; index++)
                    bestIndexes[index] = -1;

                for (int i = 0; i < units.Count; i++)
                {
                    if (localFlows[i] > 0)
                    {
                        localFlows[i] -= stepSize * iterSize;

                        virtualMixed = BuildVirtualDic(localFlows, units, mixed.Count);
                        double tmpErrorOver = ComputeOver(virtualMixed, mixed);
                        double tmpErrorUnder = ComputeUnder(virtualMixed, mixed);

                        double tmpFlowRate = Math.Abs(overError - tmpErrorOver);
                        double underDiff = 0;
                        if (tmpErrorUnder > underError)
                            underDiff = tmpErrorUnder - underError;
                        if (underDiff >= 1)
                            tmpFlowRate /= underDiff;
                        bestIndexes[i] = tmpFlowRate;

                        localFlows[i] += stepSize * iterSize;
                    }
                }

                //Pick pseudo randomly best index
                double worstFlowRate = 0.0;
                for (int index = 0; index < bestIndexes.Length; index++)
                    if (bestIndexes[index] > worstFlowRate)
                    {
                        worstFlowRate = bestIndexes[index];
                    }

                if (worstFlowRate > 0)
                {
                    int nbMatching = 0;
                    for (int index = 0; index < bestIndexes.Length; index++)
                        if (bestIndexes[index] >= worstFlowRate)
                            nbMatching++;

                    int iterChoice = rnd.Next(0, nbMatching - 1);
                    int iterNb = 0;
                    for (int index = 0; index < bestIndexes.Length; index++)
                        if (bestIndexes[index] >= worstFlowRate)
                        {
                            if (iterChoice == iterNb)
                            {
                                localFlows[index] -= stepSize * iterSize;
                                if (localFlows[index] < 0)
                                    localFlows[index] = 0.0;
                            }
                            iterNb++;
                        }
                    iterSize = 1;
                }
                else
                    iterSize++;

                virtualMixed = BuildVirtualDic(localFlows, units, mixed.Count);
                overError = ComputeOver(virtualMixed, mixed);
                underError = ComputeUnder(virtualMixed, mixed);
                if (overError + underError < bestOverallError)
                {
                    bestLocalFlows = new List<double>(localFlows);
                    bestOverallError = overError + underError;
                }
            }//End of while overflow > 1

            solution = new List<double>();
            foreach (double localFlow in localFlows)
                solution.Add(localFlow);

            underflow = underError;
        }

        private static Dictionary<double, double> BuildVirtualDic(List<double> ratios, List<Dictionary<double, double>> units, int size)
        {
            Dictionary<double, double> virtualMixed = new Dictionary<double, double>(size);
            for (int i = 0; i < ratios.Count; i++)
            {
                if (ratios[i] > 0)
                {
                    foreach (double key in units[i].Keys)
                        if (!virtualMixed.ContainsKey(key))
                            virtualMixed.Add(key, units[i][key] * ratios[i]);
                        else
                            virtualMixed[key] += units[i][key] * ratios[i];
                }
            }
            return virtualMixed;
        }

        private static double FindLocalMaxima(Dictionary<double, double> unit, Dictionary<double, double> mixed)
        {
            double minUnitRatio = double.MaxValue;
            foreach (double key in unit.Keys)
            {
                if (unit[key] > 0)
                {
                    double local = mixed[key] / unit[key];
                    if (local < minUnitRatio)
                        minUnitRatio = local;
                }
            }
            if (minUnitRatio == double.MaxValue)
                return 0;
            else
                return minUnitRatio;
        }
        
        private static double FindLocalMinima(Dictionary<double, double> unit, Dictionary<double, double> mixed, Dictionary<double, double> virtualMixed)
        {
            double minUnitRatio = double.MaxValue;
            foreach (double key in unit.Keys)
            {
                if (unit[key] > 0 && virtualMixed[key] > mixed[key])
                {
                    double local = (virtualMixed[key] - mixed[key]) / unit[key];
                    if (local < minUnitRatio)
                        minUnitRatio = local;
                }
            }
            if (minUnitRatio == double.MaxValue)
                return 0;
            else
                return minUnitRatio;
        }
    }
}
