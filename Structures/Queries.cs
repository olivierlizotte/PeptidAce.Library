/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Xml;
using PeptidAce.Utilities;

namespace PeptidAce
{
    public partial class Queries : List<Query>
    {
        public double MaxMZ;
        public double MinMZ;
        public double MaxRt;
        public double MinRt;
		public double MaxMass = 0;

        private const bool HARMONIC_CHARGE_DETECTION = false;

        //public Sample sample;
        public DBOptions dbOptions;
        public int NbSpectrum = 0;

        public Precursors Precursors;

        public Queries()
        {
            Precursors = new Precursors();
        }
        /// <summary>
        /// Constructor used for Unit tests
        /// </summary>
        /// <param name="masses"></param>
        public Queries(DBOptions options, double[] masses)
        {
            Precursors = new Precursors();
            foreach (double mass in masses)
            {
                Precursor precursor = new Precursor();
                precursor.Mass = mass;
                Add(new Query(options, null, null, precursor));
                Precursors.Add(precursor);
            }
        }

        public List<Query> ComputeAtFDR(double desired_fdr)
        {
            this.Sort(Queries.CompareScore);

            int index = Utilities.Methods.FDRizer.Extract(this, desired_fdr);
            if (index >= 0)
                return this.GetRange(0, index + 1);
            else
                return new List<Query>();
        }

        public static int CompareScore(Query left, Query right)
        {
            return -left.Score.CompareTo(right.Score);
        }

        public Queries(DBOptions dbOptions)
            : base() 
        {
            MaxMZ = double.MinValue;
            MinMZ = double.MaxValue;
            MaxRt = double.MinValue;
            MinRt = double.MaxValue;
            //this.sample = sample;
            this.dbOptions = dbOptions;
            Precursors = new Precursors();
        }

        public void GenerateQueries(Sample entry, Spectra spectra, Tracks tracks)//, double mz, double rt, double intensity)
        {
            Dictionary<Track, Precursor> Tracks = new Dictionary<Track, Precursor>();
            Dictionary<Track, Precursor> Isotopes = new Dictionary<Track, Precursor>();

            //Create one query per Spectrum-Precursor duo, including Isotopes in the process to ease search
            //For further analysis, maintain a list of precursors (excluding isotopes)
            int nbMissedTrack = 0;
            //vsSDF sdf = entry.GetSDF();// Samples.LoadSDF(entry);
            //tracks.PrepareRtSort();
            //sdf.TRACKS_LIST.PrepareRtSort();
            spectra.Sort(ProductSpectrum.AscendingPrecursorMassComparison);

            foreach (ProductSpectrum spectrum in spectra)
            {
                NbSpectrum++;
                double intensityCumul = 0.0;
                bool foundCharge = false;
                Track closestTrack = null;

                List<Query> newQueries = new List<Query>();

                //TODO No threshold on sdf files, and preferably a C# routine that does what MassSense do
                foreach (Track track in tracks.GetTracksInMzRange(spectrum.PrecursorMZ, spectrum.IsolationWindow * dbOptions.EffectiveIsolationWindowRatio))//TODO Optimize this value
                {
                    Precursor prec = null;

                    if ( track.RT_Min <= spectrum.RetentionTimeInMin &&
                         track.RT_Max >= spectrum.RetentionTimeInMin )
                    {
                        if (closestTrack == null || Math.Abs(track.MZ - spectrum.PrecursorMZ) < Math.Abs(closestTrack.MZ - spectrum.PrecursorMZ))
                            closestTrack = track;

                        if (Isotopes.ContainsKey(track))
                            break;

                        if (Tracks.ContainsKey(track))
                            prec = Tracks[track];
                        else
                        {
                            List<Precursor> isotopes = GetIsotopes(track, dbOptions, tracks, entry);
                            if (isotopes.Count > 0)
                            {
                                prec = new Precursor(track, isotopes[0].Charge, entry, 0.0, isotopes);
                                Tracks.Add(track, prec);
                                prec.OtherCharges = GetOtherCharges(prec, dbOptions, tracks, entry);

                                foreach (Precursor isotope in prec.Isotopes)
                                    if (!Isotopes.ContainsKey(isotope.Track))
                                        Isotopes.Add(isotope.Track, isotope);
                            }
                        }
                        if (prec != null)
                        {
                            intensityCumul += track.INTENSITY;
                            newQueries.Add(new Query(dbOptions, entry, spectrum, prec, NbSpectrum));

                            if (prec.Charge == spectrum.PrecursorCharge)
                                foundCharge = true;                            
                        }
                    }
                }

                if (!foundCharge)
                {
                    /*if (closestTrack != null && Tracks.ContainsKey(closestTrack) && Math.Abs(Numerics.CalculateMassError(closestTrack.MZ, spectrum.PrecursorMZ, dbOptions.precursorMassTolerance.Units)) < dbOptions.precursorMassTolerance.Value)
                    {
                        if(closestTrack.RT_Min > (float)(spectrum.RetentionTimeInMin - dbOptions.ComputedRetentionTimeDiff))
                            closestTrack.RT_Min = (float)(spectrum.RetentionTimeInMin - dbOptions.ComputedRetentionTimeDiff);
                        if (closestTrack.RT_Max < (float)(spectrum.RetentionTimeInMin + dbOptions.ComputedRetentionTimeDiff))
                            closestTrack.RT_Max = (float)(spectrum.RetentionTimeInMin + dbOptions.ComputedRetentionTimeDiff);
                        if (closestTrack.INTENSITY < spectrum.PrecursorIntensity)
                            closestTrack.INTENSITY = spectrum.PrecursorIntensity;

                        Precursor prec = Tracks[closestTrack];
                        if (prec.Charge == spectrum.PrecursorCharge)
                        {
                            Add(new Query(dbOptions, entry, spectrum, prec, NbSpectrum));
                        }
                        else
                        {
                            Precursor newPrec = new Precursor(closestTrack, spectrum.PrecursorCharge, entry);
                            Add(new Query(dbOptions, entry, spectrum, newPrec, NbSpectrum));
                        }
                    }
                    else//*/
                    {
                        nbMissedTrack++;
                        closestTrack = new Track((float)spectrum.PrecursorMZ, (float)spectrum.RetentionTimeInMin, spectrum.PrecursorIntensity,
                                                 (float)(spectrum.RetentionTimeInMin - dbOptions.ComputedRetentionTimeDiff), (float)(spectrum.RetentionTimeInMin + dbOptions.ComputedRetentionTimeDiff),
                                                 true);

                        Precursor prec = new Precursor(closestTrack, spectrum.PrecursorCharge, entry);
                        Tracks.Add(closestTrack, prec);
                        Add(new Query(dbOptions, entry, spectrum, prec, NbSpectrum));
                    }
                }//*/
                
                if (newQueries.Count > 0)
                {
                    //Remove precursors if estimated fragment intensities are too low (based on precursor intensity ratios and isolation window placement)
                    foreach (Query q in newQueries)
                    {
                        //if (q.precursor.Track.INTENSITY > intensityCumul * dbOptions.MinimumPrecursorIntensityRatioInIsolationWindow)//Need to be 5% of all intensity
                        //{
                            this.Add(q);
                        //}
                    }
                }
				//Console.Write("\r{0}%   ", ((100 * NbSpectrum) / spectra.Count));
            }
			//Console.Write("\r{0}%   ", 100);

            //Sort queries to ease search
            this.Sort(AscendingPrecursorMassComparison);
            foreach (Track track in Tracks.Keys)
                if (!Isotopes.ContainsKey(track))
                    Precursors.Add(Tracks[track]);

            //TODO Validate this approach
            //REMOVE QUERIES RELATED TO AN ISOTOPE and Compute the average CoElution 
            Dictionary<ProductSpectrum, double> DicOfSpectrumIntensities = new Dictionary<ProductSpectrum, double>();
            for(int i = 0; i < this.Count; )
			{
				Query query = this[i];
				if(query.precursor.Mass > MaxMass)
					MaxMass = query.precursor.Mass;
                if (!Isotopes.ContainsKey(query.precursor.Track))
                {
                    if (!DicOfSpectrumIntensities.ContainsKey(query.spectrum))
                        DicOfSpectrumIntensities.Add(query.spectrum, query.precursor.Track.INTENSITY);
                    else
                        DicOfSpectrumIntensities[query.spectrum] += query.precursor.Track.INTENSITY;
                    i++;
                }
                else
                    this.RemoveAt(i);
			}
			MaxMass = MassTolerance.MzTop (MaxMass, dbOptions.precursorMassTolerance);
                        
            //REMOVE Queries with Precursor intensities too low
            for (int i = 0; i < this.Count; )
            {
                Query query = this[i];
                if (query.precursor.Track.INTENSITY < DicOfSpectrumIntensities[query.spectrum] * dbOptions.MinimumPrecursorIntensityRatioInIsolationWindow)
                    this.RemoveAt(i);
                else
                    i++;
            }//*/

            Dictionary<ProductSpectrum, int> DicOfSpectrumTracks = new Dictionary<ProductSpectrum, int>();
            for (int i = 0; i < this.Count; )
            {
                Query query = this[i];
                if (!Isotopes.ContainsKey(query.precursor.Track))
                {
                    if (!DicOfSpectrumTracks.ContainsKey(query.spectrum))
                        DicOfSpectrumTracks.Add(query.spectrum, 1);
                    else
                        DicOfSpectrumTracks[query.spectrum]++;
                    i++;
                }
                else
                    this.RemoveAt(i);
            }

            double averageNbPrecursorPerSpectrum = 0;
            int nbSpectrumMatchedToTrack = 0;
            foreach (ProductSpectrum spectrum in DicOfSpectrumTracks.Keys)
            {
                nbSpectrumMatchedToTrack++;
                averageNbPrecursorPerSpectrum += DicOfSpectrumTracks[spectrum];
            }
            dbOptions.ConSole.WriteLine(entry.sSDF + " :" + Precursors.Count + " precursors [" + Isotopes.Count + " isotopes] spreaded in " + Count + " queries [" + nbMissedTrack + " trackless precursors]");
            dbOptions.ConSole.WriteLine("Average Precursors per Spectrum : " + averageNbPrecursorPerSpectrum / (double)nbSpectrumMatchedToTrack);
        }

        public static List<Precursor> GetIsotopes(Track track, DBOptions dbOptions, Tracks listTracks, Sample entry)
        {
            double isotopicMzTolerance = dbOptions.precursorMassTolerance.Value*0.5;
            if (dbOptions.precursorMassTolerance.Units == MassToleranceUnits.ppm)                
                isotopicMzTolerance = (isotopicMzTolerance / 1e6 ) * track.MZ;

            List<Precursor> bestIsotopes = new List<Precursor>();

            for (int charge = dbOptions.MaximumPrecursorChargeState; charge >= dbOptions.MinimumPrecursorChargeState; charge--)
            {
                List<Precursor> isotopes = new List<Precursor>();
                for(int nbIsotope = 1; nbIsotope < 5; nbIsotope++)
                {
                    double bestDeltaMz = isotopicMzTolerance;
                    Track bestTrack = null;
                    double massShift = Numerics.IsotopicMassShift(nbIsotope, charge);
                    double mzIsotope = track.MZ + massShift;
                    
                    foreach (Track trackToTest in listTracks.GetTracksInMzRange(mzIsotope, isotopicMzTolerance))
                    {
                        if(trackToTest.RT >= track.RT_Min &&
                           trackToTest.RT <= track.RT_Max)
                        {
                            double mzDiff = Math.Abs(mzIsotope - trackToTest.MZ);
                            if (mzDiff < bestDeltaMz)//TODO Is the best isotope the most precise one or the most intense?? Use a scoring function!
                            {
                                bestDeltaMz = mzDiff;
                                bestTrack = trackToTest;
                            }
                        }
                    }
                    if (bestTrack != null)
                        isotopes.Add(new Precursor(bestTrack, charge, entry, Constants.C12_C13_MASS_DIFFERENCE * nbIsotope, null));
                    else
                        break;
                }
                if (isotopes.Count > bestIsotopes.Count)//TODO Best way to compare isotope potentials? Number of isotopes? Delta Mass? Intensity ratios?
                    bestIsotopes = isotopes;
            }
            return bestIsotopes;
        }

        public static List<Precursor> GetOtherCharges(Precursor precursor, DBOptions dbOptions, Tracks listTracks, Sample entry)
        {
            List<Precursor> otherPrecursor = new List<Precursor>();

            for (int charge = dbOptions.MaximumPrecursorChargeState; charge >= 1; charge--)
            {
                if(charge != precursor.Charge)
                {
                    double aimedMZ = Numerics.MZFromMass(precursor.Mass, charge);

                    double chargeMzTolerance = dbOptions.precursorMassTolerance.Value*0.5;
                    if (dbOptions.precursorMassTolerance.Units == MassToleranceUnits.ppm)
                        chargeMzTolerance = (chargeMzTolerance / 1e6) * aimedMZ;

                    double bestDeltaMz = chargeMzTolerance;
                    Track bestTrack = null;
                    foreach (Track trackToTest in listTracks.GetTracksInMzRange(aimedMZ, chargeMzTolerance))
                    {
                        if (trackToTest.RT >= precursor.Track.RT_Min &&
                            trackToTest.RT <= precursor.Track.RT_Max)
                        {
                            double mzDiff = Math.Abs(aimedMZ - trackToTest.MZ);
                            if (mzDiff < bestDeltaMz)//TODO Is the best isotope the most precise one or the closest in intensity?? Use a scoring function to get both!
                            {
                                bestDeltaMz = mzDiff;
                                bestTrack = trackToTest;
                            }
                        }
                    }
                    if (bestTrack != null)
                        otherPrecursor.Add(new Precursor(bestTrack, charge, precursor.sample, 0, GetIsotopes(bestTrack, dbOptions, listTracks, entry)));
                }
            }
            return otherPrecursor;
        }

        public static int AscendingRetentionTimeComparison(Query left, Query right)
        {
            return left.spectrum.RetentionTimeInMin.CompareTo(right.spectrum.RetentionTimeInMin);
        }

        public static int AscendingPrecursorMassComparison(Query left, Query right)
        {
            return left.precursor.Mass.CompareTo(right.precursor.Mass);
        }

        public IEnumerable<Query> GetQueryInMassRange(double precursorMass, MassTolerance precursorMassTolerance)
        {
            double minimum_precursor_mass = precursorMass - precursorMassTolerance;
            double maximum_precursor_mass = precursorMass + precursorMassTolerance;
            int low_index = BinarySearch(minimum_precursor_mass);

            if (low_index >= 0 && low_index < Count && this[low_index].precursor.Mass >= minimum_precursor_mass)
                for (int i = low_index; i < Count && this[i].precursor.Mass <= maximum_precursor_mass; i++)
                    yield return this[i];
        }

        public int BinarySearch(double lowestPrecursorMass)
        {
            int low_index = 0;
            int high_index = Count - 1;
            while (low_index <= high_index)
            {
                int mid_index = low_index + ((high_index - low_index) / 2);
                int comparison = this[mid_index].precursor.Mass.CompareTo(lowestPrecursorMass);
                if (comparison == 0)
                {
                    while (mid_index > low_index && this[mid_index - 1].precursor.Mass.CompareTo(lowestPrecursorMass) == 0)
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