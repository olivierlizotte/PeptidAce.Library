/*
 * Wenger CD, Coon JJ. A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra, Journal of Proteome Research, 2013; 12(3): 1377-86
 * http://www.chem.wisc.edu/~coon/software.php#morpheus
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 * Altered by Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 */
using System.Collections.Generic;
using System.Text;

namespace PeptidAce
{
    public class Protease
    {
        public string Name { get; set; }

        public Terminus CleavageTerminus { get; set; }

        public string Sense()
        {
            return CleavageTerminus.ToString();
        }

        public string AminoAcidsInducingCleavage { get; set; }

        public string Cut()
        {
            StringBuilder cut = new StringBuilder();

            foreach(char c in AminoAcidsInducingCleavage)
            {
                cut.Append(c);
            }

            return cut.ToString();
        }

        public string AminoAcidsPreventingCleavage { get; set; }

        public string NoCut()
        {
            StringBuilder no_cut = new StringBuilder();

            foreach(char c in AminoAcidsPreventingCleavage)
            {
                no_cut.Append(c);
            }

            return no_cut.ToString();
        }

        public CleavageSpecificity CleavageSpecificity { get; set; }

        public Protease()
        {
        }
        public Protease(string name, string aminoAcidsInducingCleavage, string aminoAcidsPreventingCleavage, Terminus cleavageTerminus, CleavageSpecificity cleavageSpecificity)
        {
            Name = name;
            CleavageTerminus = cleavageTerminus;
            AminoAcidsInducingCleavage = aminoAcidsInducingCleavage;
            AminoAcidsPreventingCleavage = aminoAcidsPreventingCleavage;
            CleavageSpecificity = cleavageSpecificity;
        }

        public override string ToString()
        {
            return Name;
        }

        public List<int> GetDigestionSiteIndices(AminoAcidPolymer aminoAcidPolymer)
        {
            List<int> indices = new List<int>();

            bool inMod = false;
            for(int i = 0; i < aminoAcidPolymer.Length - 1; i++)
            {
                if (aminoAcidPolymer[i] == '(')
                    inMod = true;
                if (aminoAcidPolymer[i] == ')')
                    inMod = false;
                if (!inMod)
                {
                    foreach (char c in AminoAcidsInducingCleavage)
                    {
                        if ((CleavageTerminus != Terminus.N && aminoAcidPolymer[i] == c)
                            || (CleavageTerminus == Terminus.N && i + 1 < aminoAcidPolymer.Length && aminoAcidPolymer[i + 1] == c))
                        {
                            bool cleave = true;
                            foreach (char nc in AminoAcidsPreventingCleavage)
                            {
                                if ((CleavageTerminus != Terminus.N && i + 1 < aminoAcidPolymer.Length && aminoAcidPolymer[i + 1] == nc)
                                    || (CleavageTerminus == Terminus.N && i - 1 >= 0 && aminoAcidPolymer[i - 1] == nc))
                                {
                                    cleave = false;
                                    break;
                                }
                            }
                            if (cleave)
                            {
                                indices.Add(i);
                            }
                        }
                    }
                }
            }

            return indices;
        }
    }
}