/*
 * Wenger CD, Coon JJ. A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra, Journal of Proteome Research, 2013; 12(3): 1377-86
 * http://www.chem.wisc.edu/~coon/software.php#morpheus
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 * Altered by Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 */
using System;
using System.Collections.Generic;
using System.Xml;
using PeptidAce.Utilities;

namespace PeptidAce
{
    /// <summary>
    /// Loads and stores the sample information sored in a Raw file
    /// </summary>
    public partial class Spectra : List<ProductSpectrum>
    {
        public double MaxMZ;
        public double MinMZ;
        public double MaxProductMass;
        public double MinProductMass;
        public double MaxRt;
        public double MinRt;

        private const bool HARMONIC_CHARGE_DETECTION = false;
        public int NbScans;

        public Tracks tracks;
        public List<MS1Spectrum> MS1s;

        public Spectra(int nbScans = -1) : base() 
        {
            this.NbScans = nbScans;
            MaxMZ = double.MinValue;
            MinMZ = double.MaxValue;
            MaxProductMass = double.MinValue;
            MinProductMass = double.MaxValue;
            MaxRt = double.MinValue;
            MinRt = double.MaxValue;
            tracks = new Tracks();
            MS1s = new List<MS1Spectrum>();
        }

        public void ExportTracks(string filename)
        {
            tracks.Export(filename);
        }

        public void ExportMSMS(string filename)
        {
            vsCSVWriter writer = new vsCSVWriter(filename);
            writer.AddLine(ProductSpectrum.TITLE);

            foreach (ProductSpectrum spectrum in this)
                writer.AddLine(spectrum.ToString());
            writer.WriteToFile();
        }

        public static Spectra Import(string filenameMSMS, string filenameTracks, DBOptions dbOptions)
        {
            Spectra spectra = new Spectra();
            vsCSV csv = new vsCSV(filenameMSMS);
            if (csv.LINES_LIST.Count == 0 || csv.LINES_LIST[0].CompareTo(ProductSpectrum.TITLE) != 0)
                return null;
            for (int i = 1; i < csv.LINES_LIST.Count; i++)
            {
                string[] splits = csv.LINES_LIST[i].Split(vsCSV._Generic_Separator);
                double mz = double.Parse(splits[3]);
                int charge = int.Parse(splits[5]);
                int nbPeaks = int.Parse(splits[9]);
                List<MsMsPeak> peaks = new List<MsMsPeak>(nbPeaks);
                i++;
                for(int j = 0; j < nbPeaks; i++,j++)
                {
                    try
                    {
                        string[] splitPeaks = csv.LINES_LIST[i].Split('\t');
                        if (splitPeaks.Length > 2)
                            peaks.Add(new MsMsPeak(double.Parse(splitPeaks[0]), double.Parse(splitPeaks[1]), int.Parse(splitPeaks[2])));
                        else
                            peaks.Add(new MsMsPeak(double.Parse(splitPeaks[0]), double.Parse(splitPeaks[1]), 0));
                    }
                    catch (Exception)
                    {
                        dbOptions.ConSole.WriteLine("Error parsing line : " + csv.LINES_LIST[i]);
                    }
                }
                spectra.AddMSMS(new ProductSpectrum(int.Parse(splits[0]), double.Parse(splits[1]), splits[2], mz, double.Parse(splits[4]), charge, Utilities.Numerics.MassFromMZ(mz, charge), peaks, double.Parse(splits[8]), double.Parse(splits[10]), double.Parse(splits[11])));
            }
            if(!string.IsNullOrEmpty(filenameTracks))
                spectra.tracks = Tracks.Import(filenameTracks, dbOptions);
            return spectra;
        }

        public void AddMSMS(ProductSpectrum spectrum)
        {
            if (spectrum.PrecursorMZ > MaxMZ)
                MaxMZ = spectrum.PrecursorMZ;
            if (spectrum.PrecursorMZ < MinMZ)
                MinMZ = spectrum.PrecursorMZ;
            if (spectrum.RetentionTimeInMin > MaxRt)
                MaxRt = spectrum.RetentionTimeInMin;
            if (spectrum.RetentionTimeInMin < MinRt)
                MinRt = spectrum.RetentionTimeInMin;

            foreach (MsMsPeak peak in spectrum.Peaks)
            {
                if (peak.MZ > MaxProductMass)
                    MaxProductMass = peak.MZ;
                if (peak.MZ < MinProductMass)
                    MinProductMass = peak.MZ;
            }
            Add(spectrum);
        }

        public static Spectra Load(pwiz.CLI.msdata.MSDataFile msFile, DBOptions options, string filePath, bool loadMS = true, bool filterMS2 = true)
        {
            //Find file name in msFile;
            string mzMlFilepath = filePath;
            int num_spectra = msFile.run.spectrumList.size();
            Spectra spectra = new Spectra(num_spectra);
            //List<Trail> trails = new List<Trail>();       
            MS1Spectrum previousMS1 = null;
            try
            {

                //TODO DONT forget to remove the limiter
                //int maxNbMSMS = 10;
                double LastMs1InjectionTime = 0;
                for (int i = 0; i < num_spectra/* && i < 200*/; i++)//TODO Fix that later!
                {
                    //Spectrum
                    pwiz.CLI.msdata.Spectrum spec = msFile.run.spectrumList.spectrum(i, true);

                    if (spec.precursors.Count > 0 || spec.cvParam(pwiz.CLI.cv.CVID.MS_ms_level).value > 1)//is an MSMS
                    {
                        double retention_time = spec.scanList.scans[0].cvParam(pwiz.CLI.cv.CVID.MS_scan_start_time).timeInSeconds() / 60.0;

                        //List precursors and their intensities
                        double precursor_mz = 0;//Is there a value for the time a scan took to complete?
                        int charge = 2;
                        double precursor_intensity = 0;
                        string fragmentation_method = "unknown";
                        double isolationWindow = 1.0;
                        double injectionTime =  spec.scanList.scans[0].cvParam(pwiz.CLI.cv.CVID.MS_ion_injection_time).value;
                        if (injectionTime <= 0.0)
                            injectionTime = 120;
                        foreach (pwiz.CLI.msdata.Precursor precursor in spec.precursors)
                        {
                            fragmentation_method = precursor.activation.cvParam(pwiz.CLI.cv.CVID.MS_fragmentation_information).name;
                            if(string.IsNullOrEmpty(fragmentation_method) || fragmentation_method.Contains("nknown"))
                                if(precursor.activation.cvParams.Count > 0)
                                    fragmentation_method = precursor.activation.cvParams[0].name;

                            isolationWindow = precursor.isolationWindow.cvParam(pwiz.CLI.cv.CVID.MS_isolation_width).value;
                            //if (precursor.isolationWindow.cvParam(pwiz.CLI.cv.CVID.MS_isolation_width).value > 0)//s.Count > 2 && (double)precursor.isolationWindow.cvParams[1].value == (double)precursor.isolationWindow.cvParams[2].value)
                            //    isolationWindow = precursor.isolationWindow.cvParams[1].value;
                            //else if (precursor.isolationWindow.cvParams.Count > 2)
                            //    options.ConSole.WriteLine("Weird Isolation Window");
                            if (isolationWindow <= 0)
                                isolationWindow = 2;

                            foreach (pwiz.CLI.msdata.SelectedIon ion in precursor.selectedIons)
                            {
                                //Cycle through MS to get real precursor intensities
                                precursor_mz = ion.cvParam(pwiz.CLI.cv.CVID.MS_selected_ion_m_z).value;//.MS_isolation_width).value
                                charge = (int)ion.cvParam(pwiz.CLI.cv.CVID.MS_charge_state).value;
                                precursor_intensity = ion.cvParam(pwiz.CLI.cv.CVID.MS_peak_intensity).value;

                                if (precursor_intensity <= 0)
                                {
                                    precursor_intensity = precursor.cvParam(pwiz.CLI.cv.CVID.MS_intensity_of_precursor_ion).value;
                                    //precursor_mz = ion.cvParams[0].value;
                                    //charge = (int)ion.cvParams[1].value;
                                    //precursor_intensity = ion.cvParams[2].value;
                                }
                            }
                        }

                        int scan_index = i;
                        int scan_number = scan_index + 1;

                        pwiz.CLI.msdata.BinaryDataArray mz = spec.getMZArray();
                        pwiz.CLI.msdata.BinaryDataArray intensity = spec.getIntensityArray();

                        int num_peaks = mz.data.Count;
                        if (num_peaks != intensity.data.Count)
                        {
                            options.ConSole.WriteLine("PreoteWizard reports peaks arrays (mz/intensity) of different sizes : (" + num_peaks + "/" + intensity.data.Count + ")");
                            if (intensity.data.Count < num_peaks)
                                num_peaks = intensity.data.Count;
                        }
                        List<MsMsPeak> peaks = new List<MsMsPeak>(num_peaks);
                        for (int k = 0; k < num_peaks; k++)
                        {
                            if (intensity.data[k] > 0)
                            {
                                MsMsPeak peak = new MsMsPeak(mz.data[k], intensity.data[k], 0);
                                peaks.Add(peak);
                            }
                        }
                        mz.Dispose(); mz = null;
                        intensity.Dispose(); intensity = null;

                        peaks.Sort(MsMsPeak.AscendingMzComparison);

                        if (filterMS2)
                        {
                            //peaks = AssignChargeStates(peaks, options.maximumAssumedPrecursorChargeState, options.precursorMassTolerance);
                            //peaks = Deisotopebkp(peaks, options.maximumAssumedPrecursorChargeState, options.precursorMassTolerance);
                            peaks = AssignChargeStatesAndDeisotope(peaks, options.MaximumPrecursorChargeState, new MassTolerance(options.productMassTolerance.Value * 0.5, options.productMassTolerance.Units));
                            peaks = FilterPeaks(peaks, options.MaximumNumberOfFragmentsPerSpectrum);

                            //TODO Add Contaminant removal 
                            //peaks = ContaminantMasses.RemoveContaminantsFromMzSortedList(peaks, options.productMassTolerance);

                            //Can sometime be sorted by intensity after this call
                            //peaks = FilterPeaksV2(peaks);
                            peaks.Sort(MsMsPeak.AscendingMzComparison);
                        }

                        /*//TODO Validate that in most cases, next steps can calculate missing charge
                        if (charge == 0)
                        {
                            for (int c = options.minimumAssumedPrecursorChargeState; c <= options.maximumAssumedPrecursorChargeState; c++)
                            {
                                if (options.assignChargeStates)
                                {
                                    peaks = AssignChargeStates(peaks, c, options.productMassTolerance);
                                    if (options.deisotope)
                                    {
                                        peaks = Deisotope(peaks, c, options.productMassTolerance);
                                    }
                                }

                                double precursor_mass = Utilities.MassFromMZ(precursor_mz, c);

                                ProductSpectrum spectrum = new ProductSpectrum(mzMlFilepath, scan_number, retention_time, fragmentation_method, precursor_mz, precursor_intensity, c, precursor_mass, peaks);
                                spectra.Add(spectrum);
                            }
                        }
                        else//*/
                        {/*
                        if (options.assignChargeStates)
                        {
                            peaks = AssignChargeStatesbkp(peaks, charge, options.productMassTolerance);
                            if (options.deisotope)
                            {
                                peaks = Deisotopebkp(peaks, charge, options.productMassTolerance);
                            }
                        }//*/
                            //peaks = AssignChargeStatesAndDeisotope(peaks, options.maximumAssumedPrecursorChargeState, options.productMassTolerance);

                            double precursor_mass = Numerics.MassFromMZ(precursor_mz, charge);

                            ProductSpectrum spectrum = new ProductSpectrum(scan_number, retention_time, fragmentation_method, precursor_mz, precursor_intensity, charge, precursor_mass, peaks, isolationWindow, injectionTime, LastMs1InjectionTime);
                            spectra.AddMSMS(spectrum);
                            //zones.Add(new Zone(precursor_mz - isolationWindow, precursor_mz + isolationWindow, retention_time));
                        }

                        //if (spectra.Count >= maxNbMSMS)
                        //    i = 10000000;
                    }
                    else //Is an MS
                    {
                        LastMs1InjectionTime = spec.scanList.scans[0].cvParam(pwiz.CLI.cv.CVID.MS_ion_injection_time).value;
                        if (loadMS)
                        {
                            double retention_time = spec.scanList.scans[0].cvParam(pwiz.CLI.cv.CVID.MS_scan_start_time).timeInSeconds() / 60.0;

                            pwiz.CLI.msdata.BinaryDataArray mz = spec.getMZArray();
                            pwiz.CLI.msdata.BinaryDataArray intensity = spec.getIntensityArray();

                            if (previousMS1 != null)
                            {
                                previousMS1.ScanDuration = retention_time - previousMS1.RetentionTimeInMin;
                                spectra.MS1s.Add(previousMS1);
                            }
                            previousMS1 = new MS1Spectrum(i, retention_time, intensity.data, mz.data, 1);
                            //Trail.Follow(mz.data, intensity.data, retention_time, ref trails, options);
                            //Trail.RemoveFinished(ref trails, spectra, 1);
                        }
                    }
                    spec.Dispose(); spec = null;
                    Console.Write("\r{0}%   ", ((100 * i) / num_spectra));
                }
                if (previousMS1 != null)
                    spectra.MS1s.Add(previousMS1);

                if (spectra.MS1s.Count > 0)
                    spectra.tracks = new Tracks();//TODO ReImplement Precursor finding routine (from MassSense?)
                else
                    spectra.tracks = new Tracks();
                spectra.tracks.Sort(Tracks.AscendingPrecursorMassComparison);
                Console.Write("\r{0}%   ", 100);

                //ContaminantMasses.DisplayContaminants();
            }
            catch (Exception ex)
            {
                options.ConSole.WriteLine(ex.StackTrace);
                options.ConSole.WriteLine(ex.Message);
            }
            return spectra;
        }

        private static List<MsMsPeak> AssignChargeStates(List<MsMsPeak> peaks, int maxCharge, MassTolerance isotopicMzTolerance)
        {
            List<MsMsPeak> new_peaks = new List<MsMsPeak>(peaks);

            for (int i = 0; i < peaks.Count - 1; i++)
            {
                double massMax = (peaks[i].MZ + Constants.C12_C13_MASS_DIFFERENCE) + isotopicMzTolerance;
                int j = i + 1;
                List<int> charges = new List<int>();
                while (j < peaks.Count)
                {
                    if (peaks[j].MZ > massMax)
                        break;

                    for (int c = maxCharge; c >= 1; c--)
                    {
                        if (Math.Abs(Numerics.CalculateMassError(peaks[j].MZ, peaks[i].MZ + Constants.C12_C13_MASS_DIFFERENCE / (double)c, isotopicMzTolerance.Units)) <= isotopicMzTolerance.Value)
                        {
                            new_peaks.Add(new MsMsPeak(peaks[i].MZ, peaks[i].Intensity, c));
                            charges.Add(c);
                        }
                    }

                    j++;
                }
                if (charges.Count == 0)
                {
                    new_peaks.Add(new MsMsPeak(peaks[i].MZ, peaks[i].Intensity, 0));
                }
            }

            return new_peaks;
        }

        private static List<MsMsPeak> AssignChargeStatesbkp(List<MsMsPeak> peaks, int maxCharge, MassTolerance isotopicMzTolerance)
        {
            List<MsMsPeak> new_peaks = new List<MsMsPeak>();

            for(int i = 0; i < peaks.Count - 1; i++)
            {
                int j = i + 1;
                List<int> charges = new List<int>();
                while(j < peaks.Count)
                {
                    if(peaks[j].MZ > (peaks[i].MZ + Constants.C12_C13_MASS_DIFFERENCE) + isotopicMzTolerance)
                    {
                        break;
                    }

                    for(int c = maxCharge; c >= 1; c--)
                    {
                        if(Math.Abs(Numerics.CalculateMassError(peaks[j].MZ, peaks[i].MZ + Constants.C12_C13_MASS_DIFFERENCE / (double)c, isotopicMzTolerance.Units)) <= isotopicMzTolerance.Value)
                        {
                            new_peaks.Add(new MsMsPeak(peaks[i].MZ, peaks[i].Intensity, c));
                            charges.Add(c);
                        }
                    }

                    j++;
                }
                if(charges.Count == 0)
                {
                    new_peaks.Add(new MsMsPeak(peaks[i].MZ, peaks[i].Intensity, 0));
                }
            }

            return new_peaks;
        }

        private static List<MsMsPeak> AssignChargeStatesAndDeisotope(List<MsMsPeak> peaks, int maxCharge, MassTolerance isotopicMzTolerance)
        {
            List<MsMsPeak> new_peaks = new List<MsMsPeak>(peaks);
            //peaks.Sort(MSPeak.AscendingMzComparison);

            int[] bestIsotopes = new int[4];
            int[] currentIsotopes = new int[4];
            for (int lowMassIndex = 0; lowMassIndex < new_peaks.Count - 1; lowMassIndex++)
            {
                double bestChargeScore = 0;
                int bestCharge = 0;
                bestIsotopes[0] = 0; bestIsotopes[1] = 0; bestIsotopes[2] = 0; bestIsotopes[3] = 0;
                for(int charge = maxCharge; charge > 0; charge--)
                {
                    currentIsotopes[0] = 0; currentIsotopes[1] = 0; currentIsotopes[2] = 0; currentIsotopes[3] = 0;
                    double score = 0;
                    int potentialIsotopeIndex = lowMassIndex + 1;
                    for(int isotope = 1; isotope <= 4; isotope++)
                    {
                        double bestMassError = isotopicMzTolerance.Value;
                        double aim = Numerics.IsotopicMassShift(isotope, charge) + new_peaks[lowMassIndex].MZ;

                        while (potentialIsotopeIndex < new_peaks.Count && new_peaks[potentialIsotopeIndex].MZ < aim + bestMassError)
                        {
                            if (new_peaks[lowMassIndex].Intensity > new_peaks[potentialIsotopeIndex].Intensity)
                            {
                                double massError = Math.Abs(Numerics.CalculateMassError(new_peaks[potentialIsotopeIndex].MZ, aim, isotopicMzTolerance.Units));
                                if (massError < bestMassError)
                                {
                                    bestMassError = massError;
                                    currentIsotopes[isotope-1] = potentialIsotopeIndex;
                                }
                            }
                            potentialIsotopeIndex++;
                        }
                        score += isotopicMzTolerance.Value - bestMassError;
                        if (score == 0)
                            break;;
                    }
                    if (score > bestChargeScore)
                    {
                        bestIsotopes[0] = currentIsotopes[0];
                        bestIsotopes[1] = currentIsotopes[1];
                        bestIsotopes[2] = currentIsotopes[2];
                        bestIsotopes[3] = currentIsotopes[3];
                        bestChargeScore = score;
                        bestCharge = charge;
                    }
                }

                new_peaks[lowMassIndex].Charge = bestCharge;
                for(int i = 3; i >= 0; i--)
                    if (bestIsotopes[i] > 0)
                    {
                        new_peaks[lowMassIndex].Intensity += new_peaks[bestIsotopes[i]].Intensity;
                        new_peaks.RemoveAt(bestIsotopes[i]);
                    }                
            }
            return new_peaks;
        }

        private static List<MsMsPeak> Deisotopebkp(List<MsMsPeak> peaks, int maxCharge, MassTolerance isotopicMzTolerance)
        {
            List<MsMsPeak> new_peaks = new List<MsMsPeak>(peaks);
            peaks.Sort(MsMsPeak.AscendingMzComparison);

            for(int lowMassIndex = 0; lowMassIndex < new_peaks.Count - 1; lowMassIndex++)
            {
                if(new_peaks[lowMassIndex].Charge > 0)
                {
                    int toRemove = -1;
                    double bestMassError = isotopicMzTolerance.Value;
                    double aim = Numerics.IsotopicMassShift(1, new_peaks[lowMassIndex].Charge) + new_peaks[lowMassIndex].MZ;

                    int potentialIsotopeIndex = lowMassIndex + 1;
                    while (potentialIsotopeIndex < new_peaks.Count && new_peaks[potentialIsotopeIndex].MZ < aim + bestMassError)
                    {
                        if(new_peaks[lowMassIndex].Intensity > new_peaks[potentialIsotopeIndex].Intensity)
                        {
                            double massError = Math.Abs(Numerics.CalculateMassError(new_peaks[potentialIsotopeIndex].MZ, aim, isotopicMzTolerance.Units));
                            if(massError < bestMassError)
                            {
                                bestMassError = massError;
                                toRemove = potentialIsotopeIndex;
                            }
                        }
                        potentialIsotopeIndex++;
                    }
                    if(toRemove > 0)
                    {
                        new_peaks[lowMassIndex].Intensity += new_peaks[toRemove].Intensity;
                        new_peaks.RemoveAt(toRemove);
                    }
                }
            }
            return new_peaks;
        }

        private static List<MsMsPeak> FilterPeaks(List<MsMsPeak> peaks, int maximumNumberOfPeaks)
        {
            List<MsMsPeak> filtered_peaks = new List<MsMsPeak>(peaks);

            if (maximumNumberOfPeaks > 0 && filtered_peaks.Count > maximumNumberOfPeaks)
            {
                filtered_peaks.Sort(MsMsPeak.DescendingIntensityComparison);
                filtered_peaks.RemoveRange(maximumNumberOfPeaks, filtered_peaks.Count - maximumNumberOfPeaks);
            }

            return filtered_peaks;
        }
    }
}