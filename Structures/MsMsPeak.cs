/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using PeptidAce.Utilities;

namespace PeptidAce
{
    public class MsMsPeak
    {
        public static int Comparer(MsMsPeak p1, MsMsPeak p2)
        {
            return p1.MZ.CompareTo(p2.MZ);
        }

        private double _MZ;
        public double MZ { 
            get{ return _MZ; }  
            set
            {
                _MZ = value;
                if(this.Charge > 0)
                    this._Mass = Numerics.MassFromMZ(this._MZ, this.Charge);
                else
                    this._Mass = Numerics.MassFromMZ(this._MZ, 1);                
            }
        }

        public double Intensity { get; set; }

        public int Charge { get;  set; }

        private double _Mass;
        public double Mass
        {
            get { return _Mass; }
        }

        public MsMsPeak()
        {
        }

        public MsMsPeak(double mz, double intensity, int charge)
        {
            Intensity = intensity;
            Charge = charge;
            MZ = mz;
        }

        public override string ToString()
        {
            if(Charge > 0)
                return MZ + "\t" + Intensity + "\t" + Charge;
            else
                return MZ + "\t" + Intensity;
        }

        public static int DescendingIntensityComparison(MsMsPeak left, MsMsPeak right)
        {
            return -(left.Intensity.CompareTo(right.Intensity));
        }

        public static int AscendingMassComparison(MsMsPeak left, MsMsPeak right)
        {
            return left.Mass.CompareTo(right.Mass);
        }

        public static int AscendingMzComparison(MsMsPeak left, MsMsPeak right)
        {
            return left.MZ.CompareTo(right.MZ);
        }
    }
}