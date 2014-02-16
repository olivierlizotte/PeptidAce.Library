/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using PeptidAce.Structures;

namespace PeptidAce
{
    public class MS1Spectrum
    {
        public Spectrum Peaks { get; private set; }// Masses;

        public int ScanNumber { get; set; }

        public double RetentionTimeInMin { get; set; }
        public double ScanDuration       { get; set; }

        public double MinMz { get; set; }

        public double MaxMz { get; set; }
		/*
		public MS1Spectrum(int scanNumber, double retentionTimeInMin, pwiz.CLI.msdata.BinaryData intensities, pwiz.CLI.msdata.BinaryData masses, double scanDuration)
        {
            this.Peaks = new Spectrum(masses, intensities);
            this.ScanNumber = scanNumber;
            this.RetentionTimeInMin = retentionTimeInMin;

            this.MinMz = double.MaxValue;
            this.MaxMz = double.MinValue;
            foreach (double mass in masses)
            {
                if (mass < MinMz)
                    MinMz = mass;
                if (mass > MaxMz)
                    MaxMz = mass;
            }
            this.ScanDuration = scanDuration;
        }//*/
    }
}
