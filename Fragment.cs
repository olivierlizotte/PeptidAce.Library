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
        public static Dictionary<string, FragmentGen> AllFragments = null;
        public List<FragmentGen> fragments;
        public FullFragments(IEnumerable<FragmentGen> frags)
        {
            fragments = new List<FragmentGen>();
            foreach (FragmentGen fg in frags)
                fragments.Add(fg);
        }

        public string GetFragmentSTR(char separator = ';')
        {
            string str = fragments[0].Name;
            for (int i = 1; i < fragments.Count; i++)
                str += separator + fragments[i].Name;
            return str;
        }

        public FullFragments(bool includeCnZ = false, bool includeAnX = false, bool includeLosses = false)
        {
            if(AllFragments == null)
            {
                AllFragments = new Dictionary<string, FragmentGen>();

                AllFragments.Add("b", new FragmentGen("B", false, 0));
                AllFragments.Add("y", new FragmentGen("Y", true, Constants.WATER_MONOISOTOPIC_MASS));
                
                AllFragments.Add("a", new FragmentGen("A", false, -29.002741 + Constants.PROTON_MASS));
                AllFragments.Add("c", new FragmentGen("C", false, 17.02654915));
                AllFragments.Add("x", new FragmentGen("X", true, 43.9898346942));
                AllFragments.Add("z", new FragmentGen("Z", true, 1.991840552567 - Constants.HYDROGEN_MASS));                

                foreach (string fg in AllFragments.Keys.ToArray())
                {
                    //Water Loss
                    AllFragments.Add(fg + "-H2O", new FragmentGen(AllFragments[fg].Name + "h2oLoss", AllFragments[fg].IsReverse, AllFragments[fg].addOn - Constants.WATER_MONOISOTOPIC_MASS));
                    //Amonia Loss
                    AllFragments.Add(fg + "-A", new FragmentGen(AllFragments[fg].Name + "amoniaLoss", AllFragments[fg].IsReverse, AllFragments[fg].addOn - Constants.AMONIA_MASS));                    
                }
            }
            fragments = new List<FragmentGen>();
            fragments.Add(AllFragments["b"]);
            fragments.Add(AllFragments["y"]);
            if(includeLosses)
            {
                fragments.Add(AllFragments["b-H2O"]);
                fragments.Add(AllFragments["y-H2O"]);
                fragments.Add(AllFragments["b-A"]);
                fragments.Add(AllFragments["y-A"]);
            }

            if (includeAnX)
            {
                fragments.Add(AllFragments["a"]);
                fragments.Add(AllFragments["x"]);
                if(includeLosses)
                {
                    fragments.Add(AllFragments["a-H2O"]);
                    fragments.Add(AllFragments["x-H2O"]);
                    fragments.Add(AllFragments["a-A"]);
                    fragments.Add(AllFragments["x-A"]);
                }
            }
            if (includeCnZ)
            {
                fragments.Add(AllFragments["c"]);
                fragments.Add(AllFragments["z"]);
                if (includeLosses)
                {
                    fragments.Add(AllFragments["c-H2O"]);
                    fragments.Add(AllFragments["z-H2O"]);
                    fragments.Add(AllFragments["c-A"]);
                    fragments.Add(AllFragments["z-A"]);
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

        /// <summary>
        /// Iteratively returns all theoretical fragment masses from a given array of masses (list of Amino Acid mass + modification for each)
        /// </summary>
        /// <param name="masses"></param>
        /// <param name="precursorKnownCharge"></param>
        /// <returns></returns>
        public SortedList<double, int> ComputeFragmentFast(double[] masses, int precursorKnownCharge)
        {
            SortedList<double, int> result = new SortedList<double, int>(fragments.Count * masses.Length);
            int maxCharge = precursorKnownCharge;
            if (precursorKnownCharge > 1)
                maxCharge = precursorKnownCharge - 1;
            else
                maxCharge = 1;

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
                        product_mass = cumulC + fragment.addOn;
                    else
                        product_mass = cumulN + fragment.addOn;
                    
                    for (int c = maxCharge; c > 0; c--)
                    {
                        double mz = Numerics.MZFromMass(product_mass, c);
                        if (!result.ContainsKey(mz))
                            result.Add(mz, c);
                        else
                            result[mz] = 0;
                    }
                }
            }
            return result;
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
