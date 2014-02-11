/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System.Collections.Generic;
using PeptidAce.Utilities;

namespace PeptidAce
{
    public class ProductSpectrum
    {
        public static int COMPTEUR = 0;
        //public string Filename { get;  set; }
        public int ScanNumber { get;  set; }
        public double RetentionTimeInMin { get;  set; }
        public string FragmentationMethod { get;  set; }
        public double PrecursorMZ { get; set; }
        public double PrecursorIntensity { get; set; }
        public double PrecursorIntensityPerMilliSecond { get; set; }        
        public int PrecursorCharge { get; set; }
        public double PrecursorMass { get; set; }
        public List<MsMsPeak> Peaks { get; set; }// Masses;
        //public double[] Intensities { get; private set; }
        public double TotalIntensity { get; set; }
        public double MostIntensePeak { get; set; }
        public double IsolationWindow { get; set; }
        public double InjectionTime { get; set; }
        public double Ms1InjectionTime { get; set; }
        public int INDEX;

        public ProductSpectrum()
        {
            Peaks = new List<MsMsPeak>();
        }
        public ProductSpectrum(int scanNumber, double retentionTimeInMin, string fragmentationMethod, double precursorMZ, double precursorIntensity, int precursorCharge, double precursorMass, List<MsMsPeak> peaks, double isolationWindow, double injectionTime, double ms1InjectionTime)
        {
            this.INDEX = COMPTEUR++;
            //this.Filename = filename;
            this.ScanNumber = scanNumber;
            this.RetentionTimeInMin = retentionTimeInMin;
            this.FragmentationMethod = fragmentationMethod;
            this.PrecursorMZ = precursorMZ;
            this.PrecursorIntensity = precursorIntensity;
            this.PrecursorCharge = precursorCharge;
            this.PrecursorMass = precursorMass;
            this.IsolationWindow = isolationWindow;
            this.InjectionTime = injectionTime;
            this.Ms1InjectionTime = ms1InjectionTime;

            this.TotalIntensity = 0.0;
            this.MostIntensePeak = 0.0;
            if(peaks != null)
            {
                this.Peaks = peaks;
                for (int p = 0; p < peaks.Count; p++)
                {
                    if (peaks[p].Intensity > this.MostIntensePeak)
                        this.MostIntensePeak = peaks[p].Intensity;
                    this.TotalIntensity += peaks[p].Intensity;
                }
            }
            this.PrecursorIntensityPerMilliSecond = precursorIntensity / ms1InjectionTime;
        }
        public static string TITLE
        {
            get { return "ScanNumber,RetentionTimeInMin,FragmentationMethod,PrecursorMZ,PrecursorIntensity,PrecursorCharge,PrecursorMass,TotalIntensity,IsolationWindow,Nb Peaks,InjectionTime"; }
        }
        public override string ToString()
        {
            System.Text.StringBuilder sb = new System.Text.StringBuilder();
            sb.AppendLine(ScanNumber + "," + RetentionTimeInMin + "," + FragmentationMethod + "," + PrecursorMZ + "," + PrecursorIntensity + "," + PrecursorCharge + "," + PrecursorMass + "," + TotalIntensity + "," + IsolationWindow + "," + Peaks.Count + "," + InjectionTime + "," + Ms1InjectionTime);
            foreach (MsMsPeak peak in Peaks)
                sb.AppendLine(peak.ToString());
            return sb.ToString();
        }
        /*
        public ProductSpectrum(ProductSpectrum spectrum)
        {
            this.Filename = spectrum.Filename;
            this.ScanNumber = spectrum.ScanNumber;
            this.RetentionTimeInMin = spectrum.RetentionTimeInMin;
            this.FragmentationMethod = spectrum.FragmentationMethod;
            this.PrecursorMZ = spectrum.PrecursorMZ;
            this.PrecursorIntensity = spectrum.PrecursorIntensity;
            this.PrecursorCharge = spectrum.PrecursorCharge;
            this.PrecursorMass = spectrum.PrecursorMass;
            this.TotalIntensity = spectrum.TotalIntensity;
            this.Masses = spectrum.Masses;
            this.Intensities = spectrum.Intensities;
        }*/

        public static int AscendingPrecursorMassComparison(ProductSpectrum left, ProductSpectrum right)
        {
            return left.PrecursorMass.CompareTo(right.PrecursorMass);
        }

        public static int AscendingRetentionTimeComparison(ProductSpectrum left, ProductSpectrum right)
        {
            return left.RetentionTimeInMin.CompareTo(right.RetentionTimeInMin);
        }

        public void AlignSpectrum(double decalMz, double decalRt)
        {
            this.RetentionTimeInMin += decalRt;
            this.PrecursorMZ += decalMz;
            this.PrecursorMass = Numerics.MassFromMZ(this.PrecursorMZ, this.PrecursorCharge);
        }
        /*
        public IEnumerable<int> GetIndexOfMassesInRange(double theoMass, MassTolerance productMassTolerance)
        {
            double minimum_precursor_mass = theoMass - productMassTolerance;
            double maximum_precursor_mass = theoMass + productMassTolerance;
            int mid_index = BinarySearch(theoMass);
            if (mid_index < Peaks.Count)//.Length)
            {
                for (int i = mid_index; i >= 0 && Peaks[i].Mass >= minimum_precursor_mass; i--)
                {
                    if (Peaks[i].Mass <= maximum_precursor_mass)
                        yield return i;
                }

                for (int i = mid_index + 1; i < Peaks.Count && Peaks[i].Mass <= maximum_precursor_mass; i++)
                {
                    if (Peaks[i].Mass >= minimum_precursor_mass)
                        yield return i;
                }
            }
        }//*/
        //TODO 
        public IEnumerable<int> GetIndexOfMZInRange(double theoMz, MassTolerance productMassTolerance)
        {
            double minimum_precursor_mass = theoMz - productMassTolerance;
            double maximum_precursor_mass = theoMz + productMassTolerance;
            int mid_index = BinarySearchMz(theoMz);
            if (mid_index < Peaks.Count)//.Length)
            {
                for (int i = mid_index; i >= 0 && Peaks[i].MZ >= minimum_precursor_mass; i--)
                {
                    if (Peaks[i].MZ <= maximum_precursor_mass)
                        yield return i;
                }

                for (int i = mid_index + 1; i < Peaks.Count && Peaks[i].MZ <= maximum_precursor_mass; i++)
                {
                    if (Peaks[i].MZ >= minimum_precursor_mass)
                        yield return i;
                }
            }
        }

        private int BinarySearchMz(double theoMz)
        {
            int low_index = 0;
            int high_index = Peaks.Count - 1;// Masses.Length - 1;
            while (low_index <= high_index)
            {
                int mid_index = low_index + ((high_index - low_index) / 2);
                int comparison = Peaks[mid_index].MZ.CompareTo(theoMz);
                if (comparison == 0)
                {
                    return mid_index;
                }
                if (comparison < 0)
                {
                    low_index = mid_index + 1;
                }
                else
                {
                    high_index = mid_index - 1;
                }
            }
            return low_index;
        }
        /*
        private int BinarySearch(double theoMass)
        {
            int low_index = 0;
            int high_index = Peaks.Count - 1;// Masses.Length - 1;
            while (low_index <= high_index)
            {
                int mid_index = low_index + ((high_index - low_index) / 2);
                int comparison = Peaks[mid_index].Mass.CompareTo(theoMass);
                if (comparison == 0)
                {
                    return mid_index;
                }
                if (comparison < 0)
                {
                    low_index = mid_index + 1;
                }
                else
                {
                    high_index = mid_index - 1;
                }
            }
            return low_index;
        }//*/
    }
}