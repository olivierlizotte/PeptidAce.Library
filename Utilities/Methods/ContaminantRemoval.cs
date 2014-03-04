using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace PeptidAce.Utilities.Methods
{
    public class ContaminantRemoval
    {
        public static void FromRaw(string rawContaminants, string[] raws)
        {
            double retentiontimTol = 0.4;
            MassTolerance mt = new MassTolerance(2, MassToleranceUnits.ppm);

            //Extract MS Peaks
            Dictionary<double, List<double>> contaminants = new Dictionary<double, List<double>>();

            try
            {
                pwiz.CLI.msdata.MSDataFile msFileContaminant = new pwiz.CLI.msdata.MSDataFile(rawContaminants);

                int num_spectra = msFileContaminant.run.spectrumList.size();

                for (int i = 0; i < num_spectra; i++)
                {
                    //Spectrum
                    pwiz.CLI.msdata.Spectrum spec = msFileContaminant.run.spectrumList.spectrum(i, true);

                    if (!(spec.precursors.Count > 0 || spec.cvParam(pwiz.CLI.cv.CVID.MS_ms_level).value > 1))//is an MSMS
                    {
                        double retention_time = spec.scanList.scans[0].cvParam(pwiz.CLI.cv.CVID.MS_scan_start_time).timeInSeconds() / 60.0;

                        pwiz.CLI.msdata.BinaryDataArray mz = spec.getMZArray();
                        pwiz.CLI.msdata.BinaryDataArray intensities = spec.getIntensityArray();
                        List<double> mzOk = new List<double>();
                        for (int k = 0; k < mz.data.Count; k++)
                            if(intensities.data[k] > 500)
                                mzOk.Add(mz.data[k]);

                        if (mzOk.Count > 0)
                        {
                            mzOk.Sort();
                            contaminants.Add(retention_time, mzOk);
                        }
                    }
                    spec.Dispose(); spec = null;
                }

                msFileContaminant.Dispose();
                msFileContaminant = null;
                                
                foreach (string file in raws)
                {
                    pwiz.CLI.msdata.MSDataFile msFile = new pwiz.CLI.msdata.MSDataFile(file);

                    pwiz.CLI.msdata.SpectrumListSimple editedSpec = new pwiz.CLI.msdata.SpectrumListSimple();
                    int nbSpec = msFile.run.spectrumList.size();

                    for (int i = 0; i < nbSpec; i++)
                    {
                        //Spectrum
                        pwiz.CLI.msdata.Spectrum spec = msFile.run.spectrumList.spectrum(i, true);

                        if (!(spec.precursors.Count > 0 || spec.cvParam(pwiz.CLI.cv.CVID.MS_ms_level).value > 1))//is an MSMS
                        {
                            double retention_time = spec.scanList.scans[0].cvParam(pwiz.CLI.cv.CVID.MS_scan_start_time).timeInSeconds() / 60.0;
                            
                            List<double> timePoints = new List<double>();
                            foreach(double time in contaminants.Keys)
                            {
                                if(Math.Abs(time - retention_time) < retentiontimTol)
                                    timePoints.Add(time);
                            }

                            pwiz.CLI.msdata.BinaryDataArray mzs = spec.getMZArray();
                            pwiz.CLI.msdata.BinaryDataArray intensities = spec.getIntensityArray();

                            List<double> mzFiltered = new List<double>(mzs.data.Count);
                            List<double> intFiltered = new List<double>(intensities.data.Count);
                            for (int k = 0; k < mzs.data.Count; k++)
                            {
                                bool found = false;
                                double theMz = mzs.data[k];
                                foreach (double t in timePoints)
                                    foreach (double mz in Numerics.GetValuesInRange(contaminants[t], theMz, mt))
                                        if (Utilities.Numerics.MzDifference(mz, theMz, mt.Units) < mt.Value)
                                            found = true;
                                if (!found)
                                {
                                    mzFiltered.Add(theMz);
                                    intFiltered.Add(intensities.data[k]);
                                }
                            }
                            spec.setMZIntensityArrays(mzFiltered, intFiltered);
                        }

                        editedSpec.setDataProcessing(msFile.run.spectrumList.dataProcessing());
                        editedSpec.spectra.Add(spec);
                    }

                    pwiz.CLI.msdata.MSDataFile.WriteConfig config = new pwiz.CLI.msdata.MSDataFile.WriteConfig();
                    config.format = pwiz.CLI.msdata.MSDataFile.Format.Format_mzML;
                    //            config.compression = pwiz.CLI.msdata.MSDataFile.Compression.Compression_None;
                    //            config.precision = pwiz.CLI.msdata.MSDataFile.Precision.Precision_32;

                    msFile.run.spectrumList = editedSpec;
                    msFile.write(vsCSV.GetFolder(file) + vsCSV.GetFileName_NoExtension(file) + "_DeContaminated.mzml", config);
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.StackTrace);
            }
        }
    }
}
