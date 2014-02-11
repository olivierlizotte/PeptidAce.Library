/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace PeptidAce
{
    public class Track
    {
        public double MZ;
        public double RT;
        public double RT_Min;
        public double RT_Max;
        public double INTENSITY;
        public bool Invented;

        public Track()
        {
        }

        public Track(double mz, double rt, double intensity, double rtMin = -1, double rtMax = -1, bool invented = false)
        {
            MZ = mz;
            RT = rt;
            INTENSITY = intensity;
            if (rtMin > 0)
                RT_Min = rtMin;
            else
                RT_Min = RT;

            if (rtMax > 0)
                RT_Max = rtMax;
            else
                RT_Max = RT;
            this.Invented = invented;
        }

        public Track(Track clone)
        {
            MZ = clone.MZ;
            RT = clone.RT;
            INTENSITY = clone.INTENSITY;
        }

        public static string TITLE
        {
            get { return "MZ,RT,INTENSITY,Min Rt,Max Rt"; }
        }
        public override string ToString()
        {
            return MZ + "," + RT + "," + INTENSITY + "," + RT_Min + "," + RT_Max;
        }

        public static int _CompareMz(Track objA, Track objB)
        {
            return objA.MZ.CompareTo(objB.MZ);            
        }
    }
}
