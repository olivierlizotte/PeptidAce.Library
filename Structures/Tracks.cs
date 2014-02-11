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
    public partial class Tracks : List<Track>
    {
        public double MaxMZ;
        public double MinMZ;
        public double MaxRt;
        public double MinRt;

        private List<Track> tracksSortedByTime;
        
        public Tracks()
            : base() 
        {
            MaxMZ = double.MinValue;
            MinMZ = double.MaxValue;
            MaxRt = double.MinValue;
            MinRt = double.MaxValue;
        }

        public void Export(string filename)
        {
            vsCSVWriter writer = new vsCSVWriter(filename);
            writer.AddLine(Track.TITLE);

            foreach (Track track in this)
                writer.AddLine(track.ToString());
            writer.WriteToFile();
        }

        public static Tracks Import(string filename, DBOptions dbOptions)
        {
            vsCSV csv = new vsCSV(filename);
            if (csv.LINES_LIST.Count == 0 || csv.LINES_LIST[0].CompareTo(Track.TITLE) != 0)
                return null;
            Tracks tracks = new Tracks();
            for (int i = 1; i < csv.LINES_LIST.Count; i++)
            {
                try
                {
                    string[] splits = csv.LINES_LIST[i].Split(vsCSV._Generic_Separator);
                    tracks.AddTrack(double.Parse(splits[0]), double.Parse(splits[1]), double.Parse(splits[3]), double.Parse(splits[4]), double.Parse(splits[2]));
                }
                catch (Exception)
                {
                    dbOptions.ConSole.WriteLine("Error parsing line : " + csv.LINES_LIST[i]);
                }
            }
            return tracks;
        }

        public void AddTrack(double MZ, double RT, double RT_Min, double RT_Max, double INTENSITY, bool Invented = false)
        {
            Add(new Track(MZ, RT, INTENSITY, RT_Min, RT_Max, Invented));
            if (RT_Min < MinRt)
                MinRt = RT_Min;
            if (RT_Max > MaxRt)
                MaxRt = RT_Max;
            if (MZ < MinMZ)
                MinMZ = MZ;
            if (MZ > MaxMZ)
                MaxMZ = MZ;
        }

        public static int AscendingRetentionTimeComparison(Track left, Track right)
        {
            return left.RT.CompareTo(right.RT);
        }

        public static int AscendingPrecursorMassComparison(Track left, Track right)
        {
            return left.MZ.CompareTo(right.MZ);
        }

        public void PrepareRtSort()
        {
            tracksSortedByTime = new List<Track>(this);            
            tracksSortedByTime.Sort(AscendingRetentionTimeComparison);
        }

        public IEnumerable<Track> GetTracksInMzRange(double precursorMz, double precursorMassTolerance)
        {
            double minimum_precursor_mz = precursorMz - precursorMassTolerance;
            double maximum_precursor_mz = precursorMz + precursorMassTolerance;
            int low_index = BinarySearchMZ(minimum_precursor_mz);

            if (low_index >= 0 && low_index < Count && this[low_index].MZ >= minimum_precursor_mz)
                for (int i = low_index; i < Count && this[i].MZ <= maximum_precursor_mz; i++)
                    yield return this[i];
        }

        private int BinarySearchMZ(double lowestPrecursorMz)
        {
            int low_index = 0;
            int high_index = Count - 1;
            while (low_index <= high_index)
            {
                int mid_index = low_index + ((high_index - low_index) / 2);
                int comparison = this[mid_index].MZ.CompareTo(lowestPrecursorMz);
                if (comparison == 0)
                {
                    while (mid_index > low_index && this[mid_index - 1].MZ.CompareTo(lowestPrecursorMz) == 0)
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

        public IEnumerable<Track> GetTracksInRtRange(double minimum_precursor_rt, double maximum_precursor_rt)
        {
            int low_index = BinarySearchRT(minimum_precursor_rt);

            if (low_index >= 0 && low_index < Count && this[low_index].RT >= minimum_precursor_rt)
                for (int i = low_index; i < Count && this[i].RT <= maximum_precursor_rt; i++)
                    yield return this[i];
        }

        private int BinarySearchRT(double lowestPrecursorRt)
        {
            int low_index = 0;
            int high_index = Count - 1;
            while (low_index <= high_index)
            {
                int mid_index = low_index + ((high_index - low_index) / 2);
                int comparison = this[mid_index].RT.CompareTo(lowestPrecursorRt);
                if (comparison == 0)
                {
                    while (mid_index > low_index && this[mid_index - 1].RT.CompareTo(lowestPrecursorRt) == 0)
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