/*
 * Wenger CD, Coon JJ. A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra, Journal of Proteome Research, 2013; 12(3): 1377-86
 * http://www.chem.wisc.edu/~coon/software.php#morpheus
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 * Altered by Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 */
using PeptidAce.Utilities;

namespace PeptidAce
{
    public class MassTolerance
    {
        public double Value { get; set; }

        public MassToleranceUnits Units { get; set; }

        public MassTolerance()
        {
        }
        public MassTolerance(double value, MassToleranceUnits units)
        {
            Value = value;
            Units = units;
        }
        
        public static double MzTop(double left, MassTolerance right)
        {
            if (right.Units == MassToleranceUnits.Da)
                return left + right.Value;
            else
                return left + left * right.Value / 1e6;
        }

        public static double MzFloor(double left, MassTolerance right)
        {
            if (right.Units == MassToleranceUnits.Da)
                return left - right.Value;
            else
                return left - left * right.Value / 1e6;
        }

        public static double operator +(double left, MassTolerance right)
        {
            if(right.Units == MassToleranceUnits.Da)
            {
                return left + right.Value;
            }
            else if(right.Units == MassToleranceUnits.ppm)
            {
                return left + left * right.Value / 1e6;
            }
            else
            {
                return double.NaN;
            }
        }

        public static double operator -(double left, MassTolerance right)
        {
            if(right.Units == MassToleranceUnits.Da)
            {
                return left - right.Value;
            }
            else if(right.Units == MassToleranceUnits.ppm)
            {
                return left - left * right.Value / 1e6;
            }
            else
            {
                return double.NaN;
            }
        }
    }
}