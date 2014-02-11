/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.IO;
using PeptidAce.Utilities;

namespace PeptidAce
{
    /// <summary>
    /// Loads a list of fragments and fragment modifications from a file in the Configuration folder
    /// </summary>
    public class FragmentDictionary// : Dictionary<string, Modification>
    {
        private static Dictionary<string, Modification> instanceFragment = new Dictionary<string, Modification>();
        private static Dictionary<string, Modification> instanceAAFragment = new Dictionary<string, Modification>();

        private static FragmentDictionary Instance = new FragmentDictionary();
        private FragmentDictionary()
            : base()
        {
            using (StreamReader mods = new StreamReader(Path.Combine(Path.GetDirectoryName(Environment.GetCommandLineArgs()[0]), "Configuration", "fragments.tsv")))
            {
                string header = mods.ReadLine();

                while(mods.Peek() != -1)
                {
                    string line = mods.ReadLine();
                    string[] fields = line.Split('\t');

                    string description = fields[0];
                    ModificationType modification_type;
                    switch(fields[1])
                    {
                        case "fragment":
                            modification_type = ModificationType.Fragment;
                            break;
                        case "amino acid residue":
                            modification_type = ModificationType.AminoAcidResidue;
                            break;
                        case "peptide N-terminus":
                            modification_type = ModificationType.PeptideNTerminus;
                            break;
                        case "peptide C-terminus":
                            modification_type = ModificationType.PeptideCTerminus;
                            break;
                        default:
                            modification_type = ModificationType.AminoAcidResidue;
                            break;
                    }
                    char amino_acid = char.MinValue;
                    if(fields[2].Length == 1)
                        amino_acid = char.Parse(fields[2]);

                    double monoisotopic_mass_shift = double.Parse(fields[3]);
                    double probability = double.Parse(fields[4]);
                    string default_mod = fields[5];
                    bool default_fixed = default_mod.ToLower() == "fixed";
                    bool default_variable = default_mod.ToLower() == "variable";
                    bool automatic = default_mod.ToLower() == "auto";
                    
                    if (modification_type == ModificationType.Fragment)
                    {
                        Add(new Modification("Loss : " + description, ModificationType.AminoAcidResidue, amino_acid, -monoisotopic_mass_shift, probability, default_fixed, default_variable, automatic));
                        monoisotopic_mass_shift = Numerics.MZFromMass(monoisotopic_mass_shift, 1);
                        Add(new Modification("Fragment : " + description, modification_type, amino_acid, monoisotopic_mass_shift, probability, default_fixed, default_variable, automatic));
                    }
                    else
                        Add(new Modification(description, modification_type, amino_acid, monoisotopic_mass_shift, probability, default_fixed, default_variable, automatic));
                }
            }
        }

        public static Dictionary<string, Modification> Fragments
        {
            get { return instanceFragment; }
        }

        public static Dictionary<string, Modification> AAFragments
        {
            get { return instanceAAFragment; }
        }

        public void Add(Modification modification)
        {
            if(modification.Type == ModificationType.Fragment)
                instanceAAFragment.Add(modification.Description, modification);
            else
                instanceFragment.Add(modification.Description, modification);
        }
    }
}