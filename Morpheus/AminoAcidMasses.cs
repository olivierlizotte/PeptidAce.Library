/*
 * Wenger CD, Coon JJ. A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra, Journal of Proteome Research, 2013; 12(3): 1377-86
 * http://www.chem.wisc.edu/~coon/software.php#morpheus
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 * Altered by Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 */
using System;
using System.IO;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace PeptidAce
{
    public static class AminoAcidMasses
    {
        private static readonly double[] MONOISOTOPIC_AMINO_ACID_MASSES = new double['Z' - 'A' + 1];
        private static readonly double[] AVERAGE_AMINO_ACID_MASSES = new double['Z' - 'A' + 1];
        public static char[] AMINO_ACIDS = new char[0];
        public static char[] VALID_AMINO_ACIDS = new char[0];

        static AminoAcidMasses()
        {
            try
            {
                List<char> aminoAcids = new List<char>();
                using (StreamReader amino_acids = new StreamReader(Path.Combine(Path.GetDirectoryName(Environment.GetCommandLineArgs()[0]), "Configuration", "amino_acids.tsv")))
                {
                    string header = amino_acids.ReadLine();

                    while (amino_acids.Peek() != -1)
                    {
                        string line = amino_acids.ReadLine();
                        string[] fields = line.Split('\t');

                        char one_letter_code = char.Parse(fields[0]);                        
                        double monoisotopic_mass = double.Parse(fields[1]);
                        MONOISOTOPIC_AMINO_ACID_MASSES[one_letter_code - 'A'] = monoisotopic_mass;
                        double average_mass = double.Parse(fields[2]);
                        AVERAGE_AMINO_ACID_MASSES[one_letter_code - 'A'] = average_mass;
                        aminoAcids.Add(one_letter_code);
                    }
                }
                MONOISOTOPIC_AMINO_ACID_MASSES['X' - 'A'] = 0.0;
                AVERAGE_AMINO_ACID_MASSES['X' - 'A'] = 0.0;
                VALID_AMINO_ACIDS = aminoAcids.ToArray();

                //aminoAcids.Add('X'); 
                AMINO_ACIDS = aminoAcids.ToArray();
            }
            catch (Exception) 
            {
                Console.WriteLine("Could not load configuration file amino_acids.tsv");
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double GetMonoisotopicMass(char aminoAcid)
        {
            return MONOISOTOPIC_AMINO_ACID_MASSES[aminoAcid - 'A'];
        }

        public static double GetAverageMass(char aminoAcid)
        {
            return AVERAGE_AMINO_ACID_MASSES[aminoAcid - 'A'];
        }

        public static double GetHichestMass(int maxPeptideLength)
        {
            double maxMass = 0.0;
            foreach (double mass in MONOISOTOPIC_AMINO_ACID_MASSES)
                if (mass > maxMass)
                    maxMass = mass;
            return maxMass * maxPeptideLength;
        }

        public static double GetLowestMass(int minPeptideLength)
        {
            double minMass = double.MaxValue;
            foreach (double mass in MONOISOTOPIC_AMINO_ACID_MASSES)
                if (mass > 0.0 && mass < minMass)
                    minMass = mass;
            return minMass * minPeptideLength;
        }

        public static double GetMass(string sequence)
        {
            double cumul = PeptidAce.Utilities.Constants.WATER_MONOISOTOPIC_MASS;
            foreach (char c in sequence)
                cumul += AminoAcidMasses.GetMonoisotopicMass(c);
            return cumul;
        }

        public static double GetMass(IList<char> sequence)
        {
            double cumul = PeptidAce.Utilities.Constants.WATER_MONOISOTOPIC_MASS;
            foreach(char c in sequence)
                cumul += AminoAcidMasses.GetMonoisotopicMass(c);
            return cumul;
        }
    }
}