/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using PeptidAce.Utilities;

namespace PeptidAce
{

    public class Query : ITargetDecoy
    {
        public DBOptions                    options;
        public Sample                       sample;
        //public Dictionary<string, PeptideSpectrumMatch>       psms;
        public PeptideSpectrumMatches       psms;
		public List<PSM> psmsFast = new List<PSM>();
        public ProductSpectrum              spectrum;
        public Precursor                    precursor;
        public int                          spectrumIndex;
        public double                       TrapDistance;

        public double Score
        {            get { return ScoreFct(); }        }

        public bool Decoy
        {            get { return precursor.Decoy; }        }

        public bool Target
        {
            get { return !Decoy; }
        }

        public Query()
        {
            this.psms = new PeptideSpectrumMatches();
        }
        public Query(DBOptions dbOptions, Sample entry, ProductSpectrum spectrum, Precursor precursor, int spectraIndex = -1)
        {
            this.options    = dbOptions;
            this.sample     = entry;
            this.psms       = new PeptideSpectrumMatches();
            this.spectrum   = spectrum;
            this.precursor  = precursor;
            this.spectrumIndex = spectraIndex;
            if (spectrum != null && precursor != null)
                this.TrapDistance = Math.Abs(spectrum.PrecursorMZ - precursor.Track.MZ);
        }

        public double ScoreFct(Peptide peptide = null)
        {
            return precursor.ProbabilityScore(peptide);
        }

        public static int AscendingRetentionTimeComparison(Query left, Query right)
        {
            return left.spectrum.RetentionTimeInMin.CompareTo(right.spectrum.RetentionTimeInMin);
        }
        /*
        public double Intensity(Peptide peptide)
        {
            return precursor.GetIntensity(peptide);
        }//*/
    }
}
