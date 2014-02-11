using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace PeptidAce
{
    public class ProductMatch
    {
        public double theoMz;
        public double obsMz;
        public int charge;
        public double obsIntensity;
        public double normalizedIntensity;
        public double mass_diff;
        public FragmentClass Fragment;
        public int fragmentPos;
        public double weight;
        public ProductMatch()
        {
        }

        public ProductMatch(ProductMatch match)
        {
            this.theoMz = match.theoMz;
            this.obsMz = match.obsMz;
            this.obsIntensity = match.obsIntensity;
            this.normalizedIntensity = match.normalizedIntensity;
            this.mass_diff = match.mass_diff;
            this.charge = match.charge;
            this.Fragment = match.Fragment;
            this.fragmentPos = match.fragmentPos;
            this.weight = match.weight;
        }

        public static int AscendingWeightComparison(ProductMatch left, ProductMatch right)
        {
            if (left.weight == right.weight)
                return left.obsIntensity.CompareTo(right.obsIntensity);
            return left.weight.CompareTo(right.weight);
        }

        public static int DescendingWeightComparison(ProductMatch left, ProductMatch right)
        {
            if (left.weight == right.weight)
                return -left.obsIntensity.CompareTo(right.obsIntensity);
            return -left.weight.CompareTo(right.weight);
        }

        public override string ToString()
        {
            return "[" + obsMz + ";" + obsIntensity + "]";
        }
    }
}
