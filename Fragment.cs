/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using PeptidAce.Utilities;

namespace PeptidAce
{
    /// <summary>
    /// Abstract class for typicall fragment behavior
    /// </summary>
    public abstract class FragmentClass// : Modification
    {
        public abstract string Name { get; }
        public abstract IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul);
        public abstract bool IsReverse { get; }
        public abstract double Distribution { get; }
    }

    public class FragmentGen: FragmentClass
    {
        private string _name;
        public override string Name
        {
            get { return _name; }
        }
        public override double Distribution
        {
            get { return 1; }
        }
        private bool _isReverse;
        public override bool IsReverse
        {
            get { return _isReverse; }
        }
        public override IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul)
        {
            double mass;
            if(_isReverse)
                mass = nTermCumul + addOn;
            else
                mass = cTermCumul + addOn;
            yield return mass;
        }

        public double addOn = 0;
        public FragmentGen(string name, bool reverse, double addedMass)
        {
            _name = name;
            _isReverse = reverse;
            addOn = addedMass;
        }
    }
    /*
    public class FragmentLoss : FragmentClass//Modification, IFragment
    {
        public override string Name { get { return "loss"; } }
        override public IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul)
        {
            yield return 0;
        }
        override public bool IsReverse { get { return false; } }
        public override double Distribution
        {
            get { return 1; }
        }
    }

    public class FragmentA : FragmentClass//Modification, IFragment
    {
        public override string Name { get { return "a"; } }
        override public IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul)
        {
            yield return nTermCumul - 29.002741 + Constants.PROTON_MASS;//-CHO
        }
        override public bool IsReverse { get { return false; } }
        public override double Distribution
        {
            get { return 1; }
        }
    }

    public class FragmentB : FragmentClass//Modification, IFragment
    {
        public override string Name { get { return "b"; } }
        override public IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul)
        {
            yield return nTermCumul;
        }
        override public bool IsReverse { get { return false; } }
        public override double Distribution
        {
            get { return 1; }
        }
    }

    public class FragmentC : FragmentClass//Modification, IFragment
    {
        public override string Name { get { return "c"; } }
        override public IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul)
        {
            yield return nTermCumul + 17.02654915;
        }
        override public bool IsReverse { get { return false; } }
        public override double Distribution
        {
            get { return 1; }
        }
    }

    public class FragmentX : FragmentClass//Modification, IFragment
    {
        public override string Name { get { return "x"; } }
        override public IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul)
        {
            yield return cTermCumul + 43.9898346942;// +15.9949 + 12 + 18.0105646942 - 0.5;
            //yield return cTermCumul + 15.9949 + 12 + 18.0105646942 - 1.007276;
            //yield return cTermCumul + 26.98708959;//CO-H1.007276
        }
        override public bool IsReverse { get { return true; } }
        public override double Distribution
        {
            get { return 1; }
        }
    }

    public class FragmentY : FragmentClass//Modification, IFragment 
    {
        public override string Name { get { return "y"; } }
        override public IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul)
        {
            yield return cTermCumul + Constants.WATER_MONOISOTOPIC_MASS;
        }
        override public bool IsReverse { get { return true; } }
        public override double Distribution
        {
            get { return 1; }
        }
    }

    public class FragmentZ : FragmentClass//Modification, IFragment
    {
        public override string Name { get { return "z"; } }
        override public IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul)
        {
            yield return cTermCumul + 1.991840552567;// -Constants.PROTON_MASS;
        }
        override public bool IsReverse { get { return true; } }
        public override double Distribution
        {
            get { return 1; }
        }
    }//*/

    /// <summary>
    /// This class instantiate the fragment types to be used in a ProPheus search
    /// </summary>
    public class FullFragments
    {
        List<FragmentGen> fragments;
        public FullFragments(bool bNyOnly = false, bool includeLosses = false)
        {
            List<FragmentGen> mainFrags = new List<FragmentGen>();
            
            mainFrags.Add(new FragmentGen("B", false, 0));
            mainFrags.Add(new FragmentGen("X", true, 43.9898346942));

            if (!bNyOnly)
            {
                mainFrags.Add(new FragmentGen("A", false, -29.002741 + Constants.PROTON_MASS));
                mainFrags.Add(new FragmentGen("C", false, 17.02654915));
                mainFrags.Add(new FragmentGen("Y", true, Constants.WATER_MONOISOTOPIC_MASS));
                mainFrags.Add(new FragmentGen("Z", true, 1.991840552567));
            }

            fragments = new List<FragmentGen>();
            foreach (FragmentGen fg in mainFrags)
            {
                fragments.Add(fg);
                if (includeLosses)
                {
                    //Water Loss
                    fragments.Add(new FragmentGen(fg.Name + "h2oLoss", fg.IsReverse, fg.addOn - Constants.WATER_MONOISOTOPIC_MASS));
                    //Amonia Loss
                    fragments.Add(new FragmentGen(fg.Name + "amoniaLoss", fg.IsReverse, fg.addOn - Constants.AMONIA_MASS));
                }
            }
        }

        /// <summary>
        /// Iteratively returns all theoretical fragment masses from a given array of masses (list of Amino Acid mass + modification for each)
        /// </summary>
        /// <param name="masses"></param>
        /// <param name="precursorKnownCharge"></param>
        /// <returns></returns>
        public IEnumerable<ProductMatch> ComputeFragments(double[] masses, int precursorKnownCharge)
        {
            int maxCharge = precursorKnownCharge;
            if (precursorKnownCharge > 1)
                maxCharge = precursorKnownCharge - 1;
            else
                maxCharge = 1;

            ProductMatch match = new ProductMatch();
            double cumulN = 0.0;
            double cumulC = 0.0;
            for (int r = 0; r < masses.Length; r++)
            {
                cumulN += masses[r];
                cumulC += masses[masses.Length - r - 1];

                foreach (FragmentGen fragment in fragments)
                {
                    double product_mass;
                    if (fragment.IsReverse)
                    {
                        match.fragmentPos = masses.Length - r;
                        product_mass = cumulC + fragment.addOn;
                    }
                    else
                    {
                        match.fragmentPos = r + 1;
                        product_mass = cumulN + fragment.addOn;
                    }
                    match.Fragment = fragment;
                    match.weight = fragment.Distribution;

                    for (int c = maxCharge; c > 0; c--)
                    {
                        match.theoMz = Numerics.MZFromMass(product_mass, c);
                        match.charge = c;
                        yield return match;
                    }
                }
            }
        }
    }
    /*
    /// <summary>
    /// Compute theoretical fragments based on a peptide sequence and the precursor charge
    /// </summary>
    public class Fragments : List<FragmentClass>
    {        
        public IEnumerable<ProductMatch> ComputeFragments(double[] masses, string baseSequence, int precursorKnownCharge, DBOptions dbOptions)
        {
            int maxCharge;
            if (precursorKnownCharge > 1)
                maxCharge = precursorKnownCharge - 1;
            else
                maxCharge = 1;

            ProductMatch match = new ProductMatch();
            double cumulN = 0.0;
            double cumulC = 0.0;
            for (int r = 0; r < masses.Length; r++)
            {
                cumulN += masses[r];
                cumulC += masses[masses.Length - r - 1];

                foreach (FragmentClass fragment in this)
                {
                    if (fragment.IsReverse)
                        match.fragmentPos = masses.Length - r;
                    else
                        match.fragmentPos = r + 1;
                    match.Fragment = fragment;

                    foreach (double product_mass in fragment.ComputeFragment(cumulC, cumulN))
                    {
                        match.weight = fragment.Distribution;//TODO Adjust this value by computing overall impact (times 10?)                            
                        
                        for (int c = maxCharge; c > 0; c--)
                        {
                            match.theoMz = Numerics.MZFromMass(product_mass, c);
                            match.charge = c;
                            yield return match;
                        }
                    }
                }
            }
        }
    }//*/
}
