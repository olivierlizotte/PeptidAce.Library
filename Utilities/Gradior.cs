/**
 * Gradient implementation by Tyler Bishel : 
 * http://bleedingedgemachine.blogspot.ca/2012/12/gradient-descent.html

 * Integrated to Iso-PeptidAce by Olivier Caron-Lizotte
 * Aim: To test result on a typical Gradient Descent approach
 * This gradient descent stops when the distance between last visited coordinates is below precision threshold
 * /
 **/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PeptidAce.Utilities
{
    public class Gradior
    {
        public static List<double> Minimize(Function fct, List<double> seed, double precision, double stepSize, long maxIter, double minScore)
        {
            var gd = new Gradior(precision, stepSize, maxIter);
            return gd.Optimize(fct, seed, minScore);
        }

        public double Precision { get; set; }
        public double Eps { get; set; }
        public long MaxIterations { get; set; }
        public Gradior(double precision, double stepSize, long maxIter)
        {
            Precision = precision;// 0.001;
            Eps = stepSize;// 1;// 0.1;
            MaxIterations = maxIter;
        }
        public double Distance(List<double> a, List<double> b)
        {
            return Math.Sqrt(a.Zip(b, (i, j) => Math.Pow(i - j, 2.0)).Sum());
        }
        public delegate double Function(List<double> x);
        public delegate List<double> Derivative(List<double> x);
        
        public Derivative Gradient(Function f)
        {
            return x => Enumerable.Range(0, x.Count)
                            .Select(i => (f(x.Select((y, j) => j == i ? y + Precision : y).ToList()) -
                                          f(x.Select((y, j) => j == i ? y - Precision : y).ToList()))
                                              / (2 * Precision))                            
                            .ToList();
        }
        public List<double> Optimize(Function f, List<double> seed, double minScore)
        {
            var fPrime = Gradient(f);
            var cnt = 0;
            var x_old = seed.Select(x => x + Eps).ToList();
            var x_new = seed;
            double lastScore;// = f(x_new);
            do
            {
                x_old = x_new;
                x_new = fPrime(x_old).Zip(x_old, (i, j) => j - Eps * i).ToList();
                lastScore = f(x_new);
                cnt++;
            } while (lastScore > minScore && Distance(x_new, x_old) > Precision && cnt < MaxIterations);
            return x_new;
        }
    }
}
