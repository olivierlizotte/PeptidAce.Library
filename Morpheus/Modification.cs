/*
 * Wenger CD, Coon JJ. A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra, Journal of Proteome Research, 2013; 12(3): 1377-86
 * http://www.chem.wisc.edu/~coon/software.php#morpheus
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 * Altered by Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 */
namespace PeptidAce
{
    public class Modification
    {
        public string Description { get;  set; }

        public ModificationType Type { get; set; }

        public char AminoAcid { get; set; }

        public double MonoisotopicMassShift { get; set; }
        public double GetMonoisotopicMassShift(char aa)
        {
            if (aa == AminoAcid)
                return MonoisotopicMassShift + AminoAcidMasses.GetMonoisotopicMass(aa);
            else
                return 0;
        }

        public double Probability { get; set; }

        public bool DefaultFixedModification { get; set; }

        public bool DefaultVariableModification { get; set; }

        public bool Automatic { get; set; }

        public Modification()
        {
        }

        public Modification(string description, ModificationType type, char aminoAcid, double monoisotopicMassShift,
            double probability, bool defaultFixedModification, bool defaultVariableModification, bool automatic)
        {
            Description = description;
            Type = type;
            AminoAcid = aminoAcid;
            MonoisotopicMassShift = monoisotopicMassShift;
            Probability = probability;
            DefaultFixedModification = defaultFixedModification;
            DefaultVariableModification = defaultVariableModification;
            this.Automatic = automatic;
        }

        public override string ToString()
        {
            return Description;
        }
    }
}