using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra.Generic;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Double.Factorization;

namespace PeptidAce.Utilities
{
    public static class CurveFitter
    {
        /*
        public static void FitToSin(double[] xdata, double[] ydata, out Vector<double> p, out double a, out double b, out double c)
        {
            double minTime = double.MaxValue;
            double maxTime = double.MinValue;
            foreach (double val in xdata)
            {
                if (val > maxTime)
                    maxTime = val;
                if (val < minTime)
                    minTime = val;
            }
            var omega = 1.0 / (maxTime - minTime);

            // build matrices
            var X = DenseMatrix.OfColumnVectors(new[] {
                                new DenseVector(xdata.Length, 1),
                                new DenseVector(xdata.Select(t => Math.Sin(omega*t)).ToArray()),
                                new DenseVector(xdata.Select(t => Math.Cos(omega*t)).ToArray())});
            var y = new DenseVector(ydata);

            // solve
            p = X.QR().Solve(y);
            a = p[0];
            b = SpecialFunctions.Hypotenuse(p[1], p[2]);
            c = Math.Atan2(p[2], p[1]);
        }//*//*
        public static void FitToSin2(double[] xdata, double[] ydata, out double a, out double b, out double c)
        {// p = [ -0.287, 4.02, -1.46 ], hence f: x -> -0.287 + 4.02*sin(x) - 1.46*cos(x)
            var p = Fit.LinearCombination(xdata, ydata, z => 1.0, Math.Sin, Math.Cos);

            a = p[0];
            b = p[1];
            c = p[2];
        }//*/
        public static void FitToSin2(double[] xdata, double[] ydata, out double a, out double b, out double c)
        {// p = [ -0.287, 4.02, -1.46 ], hence f: x -> -0.287 + 4.02*sin(x) - 1.46*cos(x)
            double minTime = double.MaxValue;
            double maxTime = double.MinValue;
            foreach (double val in xdata)
            {
                if (val > maxTime)
                    maxTime = val;
                if (val < minTime)
                    minTime = val;
            }
            var omega = 1.0 / (maxTime - minTime);
            double[] xNormed = new double[xdata.Length];
            for (int i = 0; i < xdata.Length; i++)
                xNormed[i] = xdata[i] * omega;

            var p = Fit.LinearCombination(xNormed, ydata, z => 1.0, Math.Sin, Math.Exp);

            a = p[0];
            b = p[1];
            c = p[2];
        }

        public static double AreaUnderTheCurve(double xTimeStart, double timeStop, double[] coefficients)
        {
            double cumul = 0.0;
            if (timeStop > xTimeStart)
            {
                double iterSize = (timeStop - xTimeStart) / 100.0;
                for (double timePoint = xTimeStart; timePoint <= timeStop; timePoint += iterSize)
                {
                    double localIntensity = Evaluate.Polynomial(timePoint, coefficients);
                    if (localIntensity > 0)
                        cumul += localIntensity * iterSize;
                    //else
                    //    break;
                }
            }
            return cumul;
        }

        public static double AreaUnderTheCurve(IList<double> xTime, IList<double> yIntensity)
        {            
            //MathNet.Numerics.Interpolation.IInterpolation interpole = MathNet.Numerics.Interpolation.Interpolate.LinearBetweenPoints(xTime, yIntensity);

            try
            {
            MathNet.Numerics.Interpolation.IInterpolation interpole = new MathNet.Numerics.Interpolation.Algorithms.AkimaSplineInterpolation(xTime, yIntensity);
            //MathNet.Numerics.Interpolation.IInterpolation interpole = new MathNet.Numerics.Interpolation.Algorithms.CubicSplineInterpolation(xTime, yIntensity);
            double maxTime = 0;
            for (int i = 0; i < xTime.Count; i++)
                if (xTime[i] > maxTime && yIntensity[i] > 0.0)
                    maxTime = xTime[i];
                if (interpole.SupportsIntegration && maxTime > 0)
                    return interpole.Integrate(maxTime);
            }catch(Exception ex)
            {
                Console.WriteLine(ex.Message);
            }
            return 0;
            //MathNet.Numerics.Interpolation.IInterpolation interpole = MathNet.Numerics.Interpolation.Interpolate.(xTime, yIntensity);
            /*
            double minTime = double.MaxValue;
            double maxTime = double.MinValue;
            foreach (double val in xTime)
            {
                if (val > maxTime)
                    maxTime = val;
                if (val < minTime)
                    minTime = val;
            }
            return interpole.Integrate(maxTime);
            /*
            double iterSize = (maxTime - minTime) / 100.0;
            for(double i = minTime; i <= maxTime; i += iterSize)
            {
                double localIntensity = interpole.Integrate(maxTime);Evaluate.Polynomial(timePoint, coefficients);
                    if (localIntensity > 0)
                        cumul += localIntensity * iterSize;
                    //else
                    //    break;
                }
            }
            return cumul;//*/
        }
        
        public static double FitToPolynomial(double[] xdata, double[] ydata, out double[] coeff)
        {
            double minTime = double.MaxValue;
            double maxTime = double.MinValue;
            foreach (double val in xdata)
            {
                if (val > maxTime)
                    maxTime = val;
                if (val < minTime)
                    minTime = val;
            }
            var omega = 1.0 / (maxTime - minTime);
            double[] xNormed = new double[xdata.Length];
            for (int i = 0; i < xdata.Length; i++)
                xNormed[i] = xdata[i] * omega;
            
            coeff = Fit.Polynomial(xdata, ydata, 2);
            return AreaUnderTheCurve(minTime, maxTime, coeff);
        }
    }
}
