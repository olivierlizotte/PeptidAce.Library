/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace PeptidAce.Utilities
{
    public static class Constants
    {
        public const double HYDROGEN_MASS = 1.007825;
        public const double OXYGEN_MASS = 15.994915;
        public const double NITROGEN_MASS = 14.003074;
        public const double CARBON_MASS = 12.0;
        public const double PROTON_MASS = 1.00727647;
        public const double AMONIA_MASS = 17.026549;
        public const double C12_C13_MASS_DIFFERENCE = 1.0033548378;
        public const double WATER_MONOISOTOPIC_MASS = 18.0105646942;
        public const double WATER_AVERAGE_MASS = 18.01528;
        public const double PEPTIDE_N_TERMINAL_MONOISOTOPIC_MASS = 1.0078250321;
        public const double PEPTIDE_N_TERMINAL_AVERAGE_MASS = 1.00794;
        public const double PEPTIDE_C_TERMINAL_MONOISOTOPIC_MASS = 17.0027396621;
        public const double PEPTIDE_C_TERMINAL_AVERAGE_MASS = 17.00734;
    }
    /// <summary>
    /// Mass Tolerance Units enum, as expressed in the Morpheus code
    /// </summary>
    public enum MassToleranceUnits
    {
        Da,
        ppm
    }

    /// <summary>
    /// Static Routines used in Proteomic studies. 
    /// TODO In .Net 4.5, you can specify InLining, which should be benefical for these methods
    /// </summary>
    public static class Numerics
    {
        public static double StandardDeviation(IEnumerable<double> values)
        {
            return MathNet.Numerics.Statistics.Statistics.StandardDeviation(values);
        }

        public static double Variance(IEnumerable<double> values)
        {
            return MathNet.Numerics.Statistics.Statistics.Variance(values);
        }

        public static double MassFromMZ(double mz, int charge)
        {
            return mz * Math.Abs(charge) - charge * Constants.PROTON_MASS;
        }

        public static int ChargeFromMassAndMZ(double mass, double mz)
        {
            return (int)Math.Round(mass / mz);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double MZFromMass(double mass, int charge)
        {
            return (mass + charge * Constants.PROTON_MASS) / (double) Math.Abs(charge);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double MZFromMassShift(double mass, int charge)
        {
            return mass / (double)Math.Abs(charge);
        }

        public static double MZFromMzSingleCharge(double mz, int charge)
        {
            return (mz + (charge - 1) * Constants.PROTON_MASS) / (double)Math.Abs(charge);
        }

        public static double IsotopicMassShift(int nbIsotope, int charge)
        {
            return (Constants.C12_C13_MASS_DIFFERENCE * nbIsotope) / (double)charge;
        }

        public static double CalculateMassError(double experimental, double theoretical, MassToleranceUnits massErrorUnits)
        {
            if (massErrorUnits == MassToleranceUnits.Da)
            {
                return experimental - theoretical;
            }
            else if (massErrorUnits == MassToleranceUnits.ppm)
            {
                return ((experimental - theoretical) / theoretical) * 1e6;
            }
            else
            {
                return double.NaN;
            }
        }

        public static double AvgVolSideChain(char aa)
        {
            switch (aa)
            {
                case 'A': return 92;
                case 'C': return 106;
                case 'D': return 125;
                case 'E': return 155;
                case 'F': return 203;
                case 'G': return 66;
                case 'H': return 167;
                case 'I': return 169;
                case 'K': return 171;
                case 'L': return 168;
                case 'M': return 171;
                case 'N': return 135;
                case 'P': return 129;
                case 'Q': return 161;
                case 'R': return 225;
                case 'S': return 99;
                case 'T': return 122;
                case 'V': return 142;
                case 'W': return 240;
                case 'Y': return 203;
                default: return 0.0;
            }
        }

        public static double CalculateVSideChain(string sequence)
        {
            double phil = 0;
            foreach (char c in sequence)
                phil += AvgVolSideChain(c);
            return phil;
        }

        public static double Hydrophilicity(char aa)
        {
            switch(aa)
            {
                case 'A': return -0.5;
                case 'C': return -1.0;
                case 'D': return 3.0;
                case 'E': return 3.0;
                case 'F': return -2.5;
                case 'G': return 0.0;
                case 'H': return -0.5;
                case 'I': return -1.8;
                case 'K': return 3.0;
                case 'L': return -1.8;
                case 'M': return -1.3;
                case 'N': return 0.2;
                case 'P': return 0.0;
                case 'Q': return 0.2;
                case 'R': return 3.0;
                case 'S': return 0.3;
                case 'T': return -0.4;
                case 'V': return -1.5;
                case 'W': return -3.4;
                case 'Y': return -2.3;
                default: return 0.0;
            }
        }

        public static double CalculateHydrophilicity(string sequence)
        {
            double phil = 0;
            foreach (char c in sequence)
                phil += Hydrophilicity(c);
            return phil;
        }

        public static double CalculateHydrophobicity(string sequence)
        {
            int nbHydro = 0;
            foreach(char c in sequence)
            {
                switch (c)
                {
                    case 'I':
                    case 'L':
                    case 'V':
                    case 'M':
                    case 'C':
                    case 'A':
                    case 'T':
                    case 'F':
                    case 'Y':
                    case 'W':
                    case 'H':
                    case 'K':
                        nbHydro++;
                        break;
                }
            }
            return nbHydro / (double)sequence.Length;
        }

        public static double CalculateAromaticity(string sequence)
        {

            int nbAroma = 0;
            foreach (char c in sequence)
            {
                switch (c)
                {
                    case 'F':
                    case 'Y':
                    case 'W':
                    case 'H':
                        nbAroma++;
                        break;
                }
            }
            return nbAroma / (double)sequence.Length;
        }

        public static double CalculatePolarity(string sequence)
        {

            int nbPolar = 0;
            foreach (char c in sequence)
            {
                switch (c)
                {
                    case 'C':
                    case 'T':
                    case 'S':
                    case 'N':
                    case 'Q':
                    case 'D':
                    case 'E':
                    case 'H':
                    case 'K':
                    case 'R':
                    case 'W':
                    case 'Y':
                        nbPolar++;
                        break;
                }
            }
            return nbPolar / (double)sequence.Length;
        }

        public static double CalculateSulfuricity(string sequence)
        {
            int nbSulfure = 0;
            foreach (char c in sequence)
            {
                switch (c)
                {
                    case 'M':
                    case 'C':
                        nbSulfure++;
                        break;
                }
            }
            return nbSulfure / (double)sequence.Length;
        }

        public static double CalculateAliphaticity(string sequence)
        {
            int nbAli = 0;
            foreach (char c in sequence)
            {
                switch (c)
                {
                    case 'I':
                    case 'L':
                    case 'V':
                        nbAli++;
                        break;
                }
            }
            return nbAli / (double)sequence.Length;
        }

        public static double CalculateHydroxylicity(string sequence)
        {
            int nbHydroxyl= 0;
            foreach (char c in sequence)
            {
                switch (c)
                {
                    case 'T':
                    case 'S':
                        nbHydroxyl++;
                        break;
                }
            }
            return nbHydroxyl / (double)sequence.Length;
        }

        public static int CalculateAcidity(string sequence)
        {
            int nbAcid = 0;
            foreach (char c in sequence)
            {
                switch (c)
                {
                    case 'N':
                    case 'Q':
                        nbAcid++;
                        break;
                }
            }
            return nbAcid;
        }

        public static int CalculateBasic(string sequence)
        {
            int nbAcid = 0;
            foreach (char c in sequence)
            {
                switch (c)
                {
                    case 'H':
                    case 'K':
                    case 'R':
                        nbAcid++;
                        break;
                }
            }
            return nbAcid;
        }

        public static int CalculateNetCharge(string sequence)
        {
            int sum = 0;
            foreach (char c in sequence)
            {
                switch (c)
                {
                    case 'D':
                    case 'E':
                        sum--;
                        break;
                    case 'R':
                    case 'K':
                    case 'H':
                        sum++;
                        break;
                }
            }
            return sum;
        }

        public const char Asp = 'D';
        public const char Glu = 'E';
        public const char Cys = 'C';
        public const char Tyr = 'Y';
        public const char His = 'H';
        public const char Lys = 'K';
        public const char Arg = 'R';

        public static double CalculateIsoElectricPoint(string sequence)
        {
            //Source:
            //http://isoelectric.ovh.org/files/practise-isoelectric-point.html#mozTocId763352
            int ProtLength = sequence.Length;

            int AspNumber = 0;
            int GluNumber = 0;
            int CysNumber = 0;
            int TyrNumber = 0;
            int HisNumber = 0;
            int LysNumber = 0;
            int ArgNumber = 0;


            foreach (char aa in sequence)//for (i = 0; i <= protein.length() - 1; ++i)              //  we are looking for charged amino acids
            {
                switch (aa)
                {
                    case 'D':
                        ++AspNumber; break;
                    case 'E':
                        ++GluNumber; break;
                    case 'C':
                        ++CysNumber; break;
                    case 'Y':
                        ++TyrNumber; break;
                    case 'H':
                        ++HisNumber; break;
                    case 'K':
                        ++LysNumber; break;
                    case 'R':
                        ++ArgNumber; break;
                }
            }


            double NQ = 0.0; //net charge in given pH

            double QN1 = 0;  //C-terminal charge
            double QN2 = 0;  //D charge
            double QN3 = 0;  //E charge
            double QN4 = 0;  //C charge
            double QN5 = 0;  //Y charge
            double QP1 = 0;  //H charge
            double QP2 = 0;  //NH2 charge
            double QP3 = 0;  //K charge
            double QP4 = 0;  //R charge

            double pH = 0.0;

            //Normal, long method (about 650 iterations)
            //if NQ <= 0, we succeedded
            /*
            for (pH = 0.0; NQ > 0 && pH < 14.0; pH += 0.01)
            {
                // we are using pK values form Wikipedia as they give quite good approximation
                // if you want you can change it

                QN1 = -1 / (1 + Math.Pow(10, (3.65 - pH)));
                QN2 = -AspNumber / (1 + Math.Pow(10, (3.9 - pH)));
                QN3 = -GluNumber / (1 + Math.Pow(10, (4.07 - pH)));
                QN4 = -CysNumber / (1 + Math.Pow(10, (8.18 - pH)));
                QN5 = -TyrNumber / (1 + Math.Pow(10, (10.46 - pH)));
                QP1 = HisNumber / (1 + Math.Pow(10, (pH - 6.04)));
                QP2 = 1 / (1 + Math.Pow(10, (pH - 8.2)));
                QP3 = LysNumber / (1 + Math.Pow(10, (pH - 10.54)));
                QP4 = ArgNumber / (1 + Math.Pow(10, (pH - 12.48)));

                NQ = QN1 + QN2 + QN3 + QN4 + QN5 + QP1 + QP2 + QP3 + QP4;
            }//*/

            //Bisection method (about 12 iterations)
            double E = 0.01;             //epsilon means precision [pI = pH ± E]
            double pHprev = 0.0;         //of finding the solution
            double pHnext = 14.0;        //0-14 is possible pH range
            double temp = 0.0;

            for (pH = 6.5; !((pH - pHprev < E) && (pHnext - pH < E)) && pH < 14.0; )
            {
                // we are using pK values form Wikipedia as they give quite good approximation
                // if you want you can change it

                QN1 = -1 / (1 + Math.Pow(10, (3.65 - pH)));
                QN2 = -AspNumber / (1 + Math.Pow(10, (3.9 - pH)));
                QN3 = -GluNumber / (1 + Math.Pow(10, (4.07 - pH)));
                QN4 = -CysNumber / (1 + Math.Pow(10, (8.18 - pH)));
                QN5 = -TyrNumber / (1 + Math.Pow(10, (10.46 - pH)));
                QP1 = HisNumber / (1 + Math.Pow(10, (pH - 6.04)));
                QP2 = 1 / (1 + Math.Pow(10, (pH - 8.2)));
                QP3 = LysNumber / (1 + Math.Pow(10, (pH - 10.54)));
                QP4 = ArgNumber / (1 + Math.Pow(10, (pH - 12.48)));

                NQ = QN1 + QN2 + QN3 + QN4 + QN5 + QP1 + QP2 + QP3 + QP4;

                if (NQ < 0)              //we are out of range, thus the new pH value must be smaller    
                {
                    temp = pH;
                    pH = pH - ((pH - pHprev) / 2);
                    pHnext = temp;
                }
                else                  //we used to small pH value, so we have to increase it
                {
                    temp = pH;
                    pH = pH + ((pHnext - pH) / 2);
                    pHprev = temp;
                }

                if ((pH - pHprev < E) && (pHnext - pH < E)) //terminal condition, finding isoelectric point with given precision
                    break;
            }
            return pH;
        }

        public static double CalculateTolerance(double experimental, MassTolerance tolerance)
        {
            if (tolerance.Units == MassToleranceUnits.Da)
                return tolerance.Value;
            else
                return ((MassTolerance.MzTop(experimental, tolerance) - experimental) / experimental) * 1e6;
        }
        
        public static double MzDifference(double mz1, double mz2, MassToleranceUnits tol)
        {
            return Math.Abs(CalculateMassError(mz1, mz2, tol));
        }//*/
        /*
        public static double MzTop(double mz1, MassTolerance tol)
        {
            return MassTolerance.MzTop(mz1, tol);
        }//*/

        public static IEnumerable<double> GetValuesInRange(List<double> sortedMzs, double mz, MassTolerance tolerance)
        {
            double minimum_precursor_mz = MassTolerance.MzFloor(mz, tolerance);
            double maximum_precursor_mz = MassTolerance.MzTop(mz, tolerance);
            int low_index = BinarySearchMZ(sortedMzs, minimum_precursor_mz);

            if (low_index >= 0 && low_index < sortedMzs.Count && sortedMzs[low_index] >= minimum_precursor_mz)
                for (int i = low_index; i < sortedMzs.Count && sortedMzs[i] <= maximum_precursor_mz; i++)
                    yield return sortedMzs[i];
        }

        private static int BinarySearchMZ(List<double> sortedMzs, double lowestPrecursorMz)
        {
            int low_index = 0;
            int high_index = sortedMzs.Count - 1;
            while (low_index <= high_index)
            {
                int mid_index = low_index + ((high_index - low_index) / 2);
                int comparison = sortedMzs[mid_index].CompareTo(lowestPrecursorMz);
                if (comparison == 0)
                {
                    while (mid_index > low_index && sortedMzs[mid_index - 1].CompareTo(lowestPrecursorMz) == 0)
                        mid_index--;

                    return mid_index;
                }
                if (comparison < 0)
                    low_index = mid_index + 1;
                else
                    high_index = mid_index - 1;
            }
            return low_index;
        }
    }
}