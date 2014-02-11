/*
 * Wenger CD, Coon JJ. A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra, Journal of Proteome Research, 2013; 12(3): 1377-86
 * http://www.chem.wisc.edu/~coon/software.php#morpheus
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 * Altered by Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 */
using System;
using System.Collections.Generic;
using System.IO;

namespace PeptidAce
{
    public enum CleavageSpecificity
    {
        None = 0,
        Semi = 1,
        Full = 2
    }

    public class ProteaseDictionary : Dictionary<string, Protease>
    {
        private static readonly ProteaseDictionary instance = new ProteaseDictionary();

        private ProteaseDictionary()
            : base()
        {
            using(StreamReader proteases = new StreamReader(Path.Combine(Path.GetDirectoryName(Environment.GetCommandLineArgs()[0]), "Configuration", "proteases.tsv")))
            {
                string header = proteases.ReadLine();

                while(proteases.Peek() != -1)
                {
                    string line = proteases.ReadLine();
                    string[] fields = line.Split('\t');

                    string name = fields[0];
                    string amino_acids_inducing_cleavage = fields[1];
                    string amino_acids_preventing_cleavage = fields[2];
                    Terminus cleavage_terminus = (Terminus)Enum.Parse(typeof(Terminus), fields[3]);
                    CleavageSpecificity cleavage_specificity = (CleavageSpecificity)Enum.Parse(typeof(CleavageSpecificity), fields[4]);
                    Protease protease = new Protease(name, amino_acids_inducing_cleavage, amino_acids_preventing_cleavage, cleavage_terminus, cleavage_specificity);
                    Add(protease);
                }
            }
        }

        public static ProteaseDictionary Instance
        {
            get { return instance; }
        }

        public void Add(Protease protease)
        {
            base.Add(protease.Name, protease);
        }
    }
}