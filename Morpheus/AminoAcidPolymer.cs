/*
 * Wenger CD, Coon JJ. A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra, Journal of Proteome Research, 2013; 12(3): 1377-86
 * http://www.chem.wisc.edu/~coon/software.php#morpheus
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 * Altered by Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 */
using System;
using System.Collections.Generic;
using System.Text;
using System.Text.RegularExpressions;
using PeptidAce.Utilities;

namespace PeptidAce
{
    public abstract class AminoAcidPolymer
    {
        private static MassType productMassType = MassType.Monoisotopic;

        private string baseSequence;

        public string BaseSequence
        {
            get
            {
                return baseSequence;
            }
            set
            {
                baseSequence = value;
                Length = value.Length;
            }
        }

        public char this[int index]
        {
            get
            {
                return baseSequence[index];
            }
        }

        public int Length { get; private set; }

        public static double ComputeMonoisotopicMass(string pepSeqWithMods)
        {
            double mass = Constants.WATER_MONOISOTOPIC_MASS;

            bool cumulMod = false;
            string strMod = "";
            for(int i = 0; i < pepSeqWithMods.Length; i++)
            {
                if(pepSeqWithMods[i] == '(')
                {
                    strMod = "";
                    cumulMod = true;
                }
                else if(pepSeqWithMods[i] == ')')
                {
                    cumulMod = false;
                    mass += ModificationDictionary.Instance[strMod].MonoisotopicMassShift;
                } else if(cumulMod)
                    strMod += pepSeqWithMods[i];
                else
                    mass += AminoAcidMasses.GetMonoisotopicMass(pepSeqWithMods[i]);
            }
            return mass;
        }

        private double _preCalcMonoMass = 0.0;
        public double MonoisotopicMass
        {
            get
            {
                if (_preCalcMonoMass > 0)
                    return _preCalcMonoMass;
                else
                {
                    _preCalcMonoMass = ComputeMonoisotopicMass(baseSequence);

                    if (fixedModifications != null)
                        foreach (List<Modification> fixed_modifications in fixedModifications.Values)
                            foreach (Modification fixed_modification in fixed_modifications)
                                _preCalcMonoMass += fixed_modification.MonoisotopicMassShift;

                    if (variableModifications != null)
                        foreach (Modification variable_modification in variableModifications.Values)
                            _preCalcMonoMass += variable_modification.MonoisotopicMassShift;

                    return _preCalcMonoMass;
                }
            }
        }

        public string BaseLeucineSequence
        {
            get { return baseSequence.Replace('I', 'L'); }
        }

        protected string _preCompSeq = null;
        public string Sequence
        {
            get
            {
                if(_preCompSeq == null)
                {
                    StringBuilder sequence = new StringBuilder();

                    // fixed modifications on protein N-terminus
                    if(fixedModifications != null)
                    {
                        List<Modification> prot_n_term_fixed_mods;
                        if(fixedModifications.TryGetValue(0, out prot_n_term_fixed_mods))
                        {
                            foreach(Modification fixed_modification in prot_n_term_fixed_mods)
                            {
                                sequence.Append('[' + fixed_modification.Description.ToString() + ']');
                            }
                        }
                    }
                    // variable modification on protein N-terminus
                    if(variableModifications != null)
                    {
                        Modification prot_n_term_variable_mod;
                        if(variableModifications.TryGetValue(0, out prot_n_term_variable_mod))
                        {
                            sequence.Append('(' + prot_n_term_variable_mod.Description.ToString() + ')');
                        }
                    }

                    // fixed modifications on peptide N-terminus
                    if(fixedModifications != null)
                    {
                        List<Modification> pep_n_term_fixed_mods;
                        if(fixedModifications.TryGetValue(1, out pep_n_term_fixed_mods))
                        {
                            foreach(Modification fixed_modification in pep_n_term_fixed_mods)
                            {
                                sequence.Append('[' + fixed_modification.Description.ToString() + ']');
                            }
                        }
                    }
                    // variable modification on peptide N-terminus
                    if(variableModifications != null)
                    {
                        Modification pep_n_term_variable_mod;
                        if(variableModifications.TryGetValue(1, out pep_n_term_variable_mod))
                        {
                            sequence.Append('(' + pep_n_term_variable_mod.Description.ToString() + ')');
                        }
                    }

                    for(int r = 0; r < Length; r++)
                    {
                        sequence.Append(this[r]);
                        // fixed modifications on this residue
                        if(fixedModifications != null)
                        {
                            List<Modification> residue_fixed_mods;
                            if(fixedModifications.TryGetValue(r + 2, out residue_fixed_mods))
                            {
                                foreach(Modification fixed_modification in residue_fixed_mods)
                                {
                                    sequence.Append('[' + fixed_modification.Description.ToString() + ']');
                                }
                            }
                        }
                        // variable modification on this residue
                        if(variableModifications != null)
                        {
                            Modification residue_variable_mod;
                            if(variableModifications.TryGetValue(r + 2, out residue_variable_mod))
                            {
                                sequence.Append('(' + residue_variable_mod.Description.ToString() + ')');
                            }
                        }
                    }

                    // fixed modifications on peptide C-terminus
                    if(fixedModifications != null)
                    {
                        List<Modification> pep_c_term_fixed_mods;
                        if(fixedModifications.TryGetValue(Length + 2, out pep_c_term_fixed_mods))
                        {
                            foreach(Modification fixed_modification in pep_c_term_fixed_mods)
                            {
                                sequence.Append('[' + fixed_modification.Description.ToString() + ']');
                            }
                        }
                    }
                    // variable modification on peptide C-terminus
                    if(variableModifications != null)
                    {
                        Modification pep_c_term_variable_mod;
                        if(variableModifications.TryGetValue(Length + 2, out pep_c_term_variable_mod))
                        {
                            sequence.Append('(' + pep_c_term_variable_mod.Description.ToString() + ')');
                        }
                    }

                    // fixed modifications on protein C-terminus
                    if(fixedModifications != null)
                    {
                        List<Modification> prot_c_term_fixed_mods;
                        if(fixedModifications.TryGetValue(Length + 3, out prot_c_term_fixed_mods))
                        {
                            foreach(Modification fixed_modification in prot_c_term_fixed_mods)
                            {
                                sequence.Append('[' + fixed_modification.Description.ToString() + ']');
                            }
                        }
                    }
                    // variable modification on protein C-terminus
                    if(variableModifications != null)
                    {
                        Modification prot_c_term_variable_mod;
                        if(variableModifications.TryGetValue(Length + 3, out prot_c_term_variable_mod))
                        {
                            sequence.Append('(' + prot_c_term_variable_mod.Description.ToString() + ')');
                        }
                    }
                    _preCompSeq = sequence.ToString();                    
                }
                return _preCompSeq;
            }
        }

        public string LeucineSequence
        {
            get { return Sequence.Replace('I', 'L'); }
        }

        private static readonly Regex INVALID_AMINO_ACIDS = new Regex("[^ACDEFGHIKLMNPQRSTVWY]");

        public AminoAcidPolymer()
        {
            BaseSequence = "";
            variableModifications = new Dictionary<int, Modification>();
            fixedModifications = new Dictionary<int, List<Modification>>();
        }
        protected AminoAcidPolymer(string baseSequence, bool checkMods)
        {            
            string[] seqSplits = baseSequence.Split(new char[]{ '(', ')' });            
            string seq = "";
            for (int i = 0; i < seqSplits.Length; i += 2)
                seq += seqSplits[i];
            BaseSequence = INVALID_AMINO_ACIDS.Replace(seq, string.Empty);
        }

        protected AminoAcidPolymer(string baseSequence)
        {
            BaseSequence = baseSequence;
        }

        public override string ToString()
        {
            return Sequence;
        }
        
        private bool initializeProductArrays = true;

        protected Dictionary<int, List<Modification>> fixedModifications;

        public Dictionary<int, List<Modification>> FixedModifications
        {
            get { return fixedModifications; }
        }
        public void SetFixedModifications(Dictionary<int, List<Modification>> value)
        {
            fixedModifications = value;
            initializeProductArrays = true;
        }
        private static string StringifyListMods(Dictionary<int, List<Modification>> mods)
        {
            if (mods == null)
                return null;
            string str = "";
            foreach (int key in mods.Keys)
            {
                str += "|" + key;
                foreach (Modification mod in mods[key])
                    str += "," + mod.Description;
            }
            if (string.IsNullOrEmpty(str))
                return "";
            else
                return str.Substring(1);
        }

        private static Dictionary<int, List<Modification>> UnStringifyListMods(string str)
        {
            Dictionary<int, List<Modification>> dic = new Dictionary<int, List<Modification>>();
            string[] mods = str.Split('|');
            foreach (string mod in mods)
            {
                string[] pair = mod.Split(',');
                List<Modification> newMods = new List<Modification>();
                for(int i = 1; i < pair.Length; i++)
                    newMods.Add(ModificationDictionary.Instance[pair[i]]);
                dic.Add(int.Parse(pair[0]), newMods);
            }
            return dic;
        }
        public string FixedModificationsInString
        {
            get { return StringifyListMods(fixedModifications); }
            set
            {
                SetFixedModifications(UnStringifyListMods(value));
            }
        }

        private Dictionary<int, Modification> variableModifications;

        public Dictionary<int, Modification> VariableModifications
        {
            get { return variableModifications; }
        }

        public bool IsPionylated()
        {
            return (variableModifications != null && variableModifications.ContainsValue(ModificationDictionary.Pionylation));
        }

        public void SetVariableModifications(Dictionary<int, Modification> value)
        {
            variableModifications = new Dictionary<int,Modification>(value);
            initializeProductArrays = true;        
        }
        private static string StringifyMods(Dictionary<int, Modification> mods)
        {
            if (mods == null)
                return null;
            string str = "";
            foreach (int key in mods.Keys)
                str += "|" + key + "," + mods[key].Description;
            if (string.IsNullOrEmpty(str))
                return "";
            else
                return str.Substring(1);
        }

        private static Dictionary<int, Modification> UnStringifyMods(string str)
        {
            Dictionary<int, Modification> dic = new Dictionary<int,Modification>();
            string[] mods = str.Split('|');
            foreach (string mod in mods)
            {
                string[] pair = mod.Split(',');
                dic.Add(int.Parse(pair[0]), ModificationDictionary.Instance[pair[1]]);
            }
            return dic;
        }
        public string VariableModificationsInString
        {
            get { return StringifyMods(variableModifications); }
            set
            {
                SetVariableModifications(UnStringifyMods(value));
            }
        }

        public void SetFixedModifications(IEnumerable<Modification> fixedModifications)
        {
            this.fixedModifications = new Dictionary<int, List<Modification>>(Length + 4);

            foreach(Modification fixed_modification in fixedModifications)
            {
                if(fixed_modification.Type == ModificationType.ProteinNTerminus && (this is Protein ||
                    (this is Peptide && (((Peptide)this).StartResidueNumber == 1 || (((Peptide)this).StartResidueNumber == 2 && ((Peptide)this).Parent[0] == 'M')))))
                {
                    List<Modification> prot_n_term_fixed_mods;
                    if(!this.fixedModifications.TryGetValue(0, out prot_n_term_fixed_mods))
                    {
                        prot_n_term_fixed_mods = new List<Modification>();
                        prot_n_term_fixed_mods.Add(fixed_modification);
                        this.fixedModifications.Add(0, prot_n_term_fixed_mods);
                    }
                    else
                    {
                        prot_n_term_fixed_mods.Add(fixed_modification);
                    }
                }

                if(fixed_modification.Type == ModificationType.PeptideNTerminus)
                {
                    List<Modification> pep_n_term_fixed_mods;
                    if(!this.fixedModifications.TryGetValue(1, out pep_n_term_fixed_mods))
                    {
                        pep_n_term_fixed_mods = new List<Modification>();
                        pep_n_term_fixed_mods.Add(fixed_modification);
                        this.fixedModifications.Add(1, pep_n_term_fixed_mods);
                    }
                    else
                    {
                        pep_n_term_fixed_mods.Add(fixed_modification);
                    }
                }

                for(int r = 0; r < Length; r++)
                {
                    if(fixed_modification.Type == ModificationType.AminoAcidResidue && this[r] == fixed_modification.AminoAcid && (variableModifications == null || !variableModifications.ContainsKey(r + 2)))
                    {
                        List<Modification> residue_fixed_mods;
                        if(!this.fixedModifications.TryGetValue(r + 2, out residue_fixed_mods))
                        {
                            residue_fixed_mods = new List<Modification>();
                            residue_fixed_mods.Add(fixed_modification);
                            this.fixedModifications.Add(r + 2, residue_fixed_mods);
                        }
                        else
                        {
                            residue_fixed_mods.Add(fixed_modification);
                        }
                    }
                }

                if(fixed_modification.Type == ModificationType.PeptideCTerminus)
                {
                    List<Modification> pep_c_term_fixed_mods;
                    if(!this.fixedModifications.TryGetValue(Length + 2, out pep_c_term_fixed_mods))
                    {
                        pep_c_term_fixed_mods = new List<Modification>();
                        pep_c_term_fixed_mods.Add(fixed_modification);
                        this.fixedModifications.Add(Length + 2, pep_c_term_fixed_mods);
                    }
                    else
                    {
                        pep_c_term_fixed_mods.Add(fixed_modification);
                    }
                }

                if(fixed_modification.Type == ModificationType.ProteinCTerminus && (this is Protein || (this is Peptide && ((Peptide)this).EndResidueNumber == ((Peptide)this).Parent.Length - 1)))
                {
                    List<Modification> prot_c_term_fixed_mods;
                    if(!this.fixedModifications.TryGetValue(Length + 3, out prot_c_term_fixed_mods))
                    {
                        prot_c_term_fixed_mods = new List<Modification>();
                        prot_c_term_fixed_mods.Add(fixed_modification);
                        this.fixedModifications.Add(Length + 3, prot_c_term_fixed_mods);
                    }
                    else
                    {
                        prot_c_term_fixed_mods.Add(fixed_modification);
                    }
                }
            }

            if(this.fixedModifications.Count == 0)
            {
                this.fixedModifications = null;
            }

            initializeProductArrays = true;
        }

        private double[] cumulativeNTerminalMass;
        private double[] cumulativeCTerminalMass;

        private void InitializeProductArrays()
        {
            double mass_shift;

            cumulativeNTerminalMass = new double[Length];

            mass_shift = 0.0;
            // fixed modifications on protein N-terminus
            if(fixedModifications != null)
            {
                List<Modification> prot_n_term_fixed_mods;
                if(fixedModifications.TryGetValue(0, out prot_n_term_fixed_mods))
                {
                    foreach(Modification fixed_modification in prot_n_term_fixed_mods)
                    {
                        mass_shift += fixed_modification.MonoisotopicMassShift;
                    }
                }
            }
            // variable modification on the protein N-terminus
            if(variableModifications != null)
            {
                Modification protein_n_term_variable_mod;
                if(variableModifications.TryGetValue(0, out protein_n_term_variable_mod))
                {
                    mass_shift += protein_n_term_variable_mod.MonoisotopicMassShift;
                }
            }
            // fixed modifications on peptide N-terminus
            if(fixedModifications != null)
            {
                List<Modification> pep_n_term_fixed_mods;
                if(fixedModifications.TryGetValue(1, out pep_n_term_fixed_mods))
                {
                    foreach(Modification fixed_modification in pep_n_term_fixed_mods)
                    {
                        mass_shift += fixed_modification.MonoisotopicMassShift;
                    }
                }
            }
            // variable modification on peptide N-terminus
            if(variableModifications != null)
            {
                Modification pep_n_term_variable_mod;
                if(variableModifications.TryGetValue(1, out pep_n_term_variable_mod))
                {
                    mass_shift += pep_n_term_variable_mod.MonoisotopicMassShift;
                }
            }
            cumulativeNTerminalMass[0] = mass_shift;

            for(int r = 1; r < Length; r++)
            {
                mass_shift = 0.0;
                // fixed modifications on this residue
                if(fixedModifications != null)
                {
                    List<Modification> residue_fixed_mods;
                    if(fixedModifications.TryGetValue(r + 1, out residue_fixed_mods))
                    {
                        foreach(Modification fixed_modification in residue_fixed_mods)
                        {
                            mass_shift += fixed_modification.MonoisotopicMassShift;
                        }
                    }
                }
                // variable modification on this residue
                if(variableModifications != null)
                {
                    Modification residue_variable_mod;
                    if(variableModifications.TryGetValue(r + 1, out residue_variable_mod))
                    {
                        mass_shift += residue_variable_mod.MonoisotopicMassShift;
                    }
                }
                cumulativeNTerminalMass[r] = cumulativeNTerminalMass[r - 1] + (productMassType == MassType.Average ? AminoAcidMasses.GetAverageMass(this[r - 1]) : AminoAcidMasses.GetMonoisotopicMass(this[r - 1])) + mass_shift;
            }

            cumulativeCTerminalMass = new double[Length];

            mass_shift = 0.0;
            // fixed modifications on protein C-terminus
            if(fixedModifications != null)
            {
                List<Modification> prot_c_term_fixed_mods;
                if(fixedModifications.TryGetValue(Length + 3, out prot_c_term_fixed_mods))
                {
                    foreach(Modification fixed_modification in prot_c_term_fixed_mods)
                    {
                        mass_shift += fixed_modification.MonoisotopicMassShift;
                    }
                }
            }
            // variable modification on protein C-terminus
            if(variableModifications != null)
            {
                Modification prot_c_term_variable_mod;
                if(variableModifications.TryGetValue(Length + 3, out prot_c_term_variable_mod))
                {
                    mass_shift += prot_c_term_variable_mod.MonoisotopicMassShift;
                }
            }
            // fixed modifications on peptide C-terminus
            if(fixedModifications != null)
            {
                List<Modification> pep_c_term_fixed_mods;
                if(fixedModifications.TryGetValue(Length + 2, out pep_c_term_fixed_mods))
                {
                    foreach(Modification fixed_modification in pep_c_term_fixed_mods)
                    {
                        mass_shift += fixed_modification.MonoisotopicMassShift;
                    }
                }
            }
            // variable modification on peptide C-terminus
            if(variableModifications != null)
            {
                Modification pep_c_term_variable_mod;
                if(variableModifications.TryGetValue(Length + 2, out pep_c_term_variable_mod))
                {
                    mass_shift += pep_c_term_variable_mod.MonoisotopicMassShift;
                }
            }
            cumulativeCTerminalMass[0] = mass_shift;

            for(int r = 1; r < Length; r++)
            {
                mass_shift = 0.0;
                // fixed modifications on this residue
                if(fixedModifications != null)
                {
                    List<Modification> residue_fixed_mods;
                    if(fixedModifications.TryGetValue(Length - r + 2, out residue_fixed_mods))
                    {
                        foreach(Modification fixed_modification in residue_fixed_mods)
                        {
                            mass_shift += fixed_modification.MonoisotopicMassShift;
                        }
                    }
                }
                // variable modification on this residue
                if(variableModifications != null)
                {
                    Modification residue_variable_mod;
                    if(variableModifications.TryGetValue(Length - r + 2, out residue_variable_mod))
                    {
                        mass_shift += residue_variable_mod.MonoisotopicMassShift;
                    }
                }

                cumulativeCTerminalMass[r] = cumulativeCTerminalMass[r - 1] + (productMassType == MassType.Average ? AminoAcidMasses.GetAverageMass(this[Length - r]) : AminoAcidMasses.GetMonoisotopicMass(this[Length - r])) + mass_shift;
            }

            initializeProductArrays = false;
        }

        private static readonly ProductCaps PRODUCT_CAPS = ProductCaps.Instance;

        public double CalculateProductMass(ProductType productType, int productNumber)
        {
            if(initializeProductArrays)
            {
                InitializeProductArrays();
            }

            switch(productType)
            {
                case ProductType.b:
                case ProductType.c:
                    return cumulativeNTerminalMass[productNumber] + PRODUCT_CAPS[productType, productMassType];
                case ProductType.y:
                case ProductType.zdot:
                    return cumulativeCTerminalMass[productNumber] + PRODUCT_CAPS[productType, productMassType];
                default:
                    return double.NaN;
            }
        }

        public Product CalculateProduct(ProductType productType, int productNumber)
        {
            if(initializeProductArrays)
            {
                InitializeProductArrays();
            }

            switch(productType)
            {
                case ProductType.b:
                case ProductType.c:
                    return new Product(productType, productNumber, cumulativeNTerminalMass[productNumber] + PRODUCT_CAPS[productType, productMassType]);
                case ProductType.y:
                case ProductType.zdot:
                    return new Product(productType, productNumber, cumulativeCTerminalMass[productNumber] + PRODUCT_CAPS[productType, productMassType]);
                default:
                    return null;
            }
        }

        protected IEnumerable<Dictionary<int, Modification>> GetModificationDic(Dictionary<int, List<Modification>> possibleMods, int maxModPerPeptide, Dictionary<int, Modification> dic)
        {
            if (dic.Count < maxModPerPeptide)
            {
                foreach (int key in possibleMods.Keys)
                {
                    if(!dic.ContainsKey(key))
                    {
                        List<Modification> list = possibleMods[key];
                        for (int i = 0; i < list.Count; i++)
                        {
                            dic.Add(key, list[i]);
                            yield return dic;
                            foreach (Dictionary<int, Modification> otherDic in GetModificationDic(possibleMods, maxModPerPeptide, dic))
                                yield return otherDic;
                            dic.Remove(key);
                        }
                    }
                }
            }
        }

        protected IEnumerable<Dictionary<int, Modification>> GetVariableModificationPatterns(Dictionary<int, List<Modification>> possibleVariableModifications, int maxModPerPeptide)
        {
            if(possibleVariableModifications.Count > 0 && maxModPerPeptide > 0)
            {
                foreach (Dictionary<int, Modification> dic in GetModificationDic(possibleVariableModifications, maxModPerPeptide, new Dictionary<int, Modification>()))
                    yield return dic;
                
                //List<KeyValuePair<int, List<Modification>>> possible_variable_modifications = new List<KeyValuePair<int, List<Modification>>>(possibleVariableModifications);
                //int[] base_variable_modification_pattern = new int[Length + 4];
                //for (int variable_modifications = 0; variable_modifications <= possibleVariableModifications.Count; variable_modifications++)
                //{
                //    foreach (int[] variable_modification_pattern in GetVariableModificationPatterns(possibleVariableModifications, possibleVariableModifications.Count - variable_modifications, base_variable_modification_pattern, 0, maxModPerPeptide))
                //    {
                   //     Dictionary<int, Modification> dic = GetVariableModificationPattern(variable_modification_pattern, possibleVariableModifications);
                     //   if (dic.Count > 0)
                       //     yield return dic;
                //    }
                //}
            }
        }
        /*
        protected IEnumerable<Dictionary<int, Modification>> GetVariableModificationPatterns(Dictionary<int, List<Modification>> possibleVariableModifications, int maxModPerPeptide)
        {
            if(possibleVariableModifications.Count > 0)
            {
                List<KeyValuePair<int, List<Modification>>> possible_variable_modifications = new List<KeyValuePair<int, List<Modification>>>(possibleVariableModifications);
                int[] base_variable_modification_pattern = new int[Length + 4];
                for(int variable_modifications = 0; variable_modifications <= possible_variable_modifications.Count; variable_modifications++)
                {
                    foreach (int[] variable_modification_pattern in GetVariableModificationPatterns(possible_variable_modifications, possible_variable_modifications.Count - variable_modifications, base_variable_modification_pattern, 0, maxModPerPeptide))
                    {
                        Dictionary<int, Modification> dic = GetVariableModificationPattern(variable_modification_pattern, possibleVariableModifications);
                        if (dic.Count > 0)
                            yield return dic;
                    }
                }
            }
        }//*/
        private static IEnumerable<int[]> GetVariableModificationPatterns(List<KeyValuePair<int, List<Modification>>> possibleVariableModifications, int unmodifiedResiduesDesired, int[] variableModificationPattern, int index, int maxModPerPeptide)
        {
            if (index < possibleVariableModifications.Count - 1 && index < maxModPerPeptide)
            {
                if(unmodifiedResiduesDesired > 0)
                {
                    variableModificationPattern[possibleVariableModifications[index].Key] = 0;
                    foreach (int[] new_variable_modification_pattern in GetVariableModificationPatterns(possibleVariableModifications, unmodifiedResiduesDesired - 1, variableModificationPattern, index + 1, maxModPerPeptide))
                    {
                        yield return new_variable_modification_pattern;
                    }
                }
                if(unmodifiedResiduesDesired < possibleVariableModifications.Count - index)
                {
                    for(int i = 1; i <= possibleVariableModifications[index].Value.Count; i++)
                    {
                        variableModificationPattern[possibleVariableModifications[index].Key] = i;
                        foreach (int[] new_variable_modification_pattern in GetVariableModificationPatterns(possibleVariableModifications, unmodifiedResiduesDesired, variableModificationPattern, index + 1, maxModPerPeptide))
                        {
                            yield return new_variable_modification_pattern;
                        }
                    }
                }
            }
            else
            {
                if(unmodifiedResiduesDesired > 0)
                {
                    variableModificationPattern[possibleVariableModifications[index].Key] = 0;
                    yield return variableModificationPattern;
                }
                else
                {
                    for(int i = 1; i <= possibleVariableModifications[index].Value.Count; i++)
                    {
                        variableModificationPattern[possibleVariableModifications[index].Key] = i;
                        yield return variableModificationPattern;
                    }
                }
            }
        }

        private static Dictionary<int, Modification> GetVariableModificationPattern(int[] variableModificationArray, Dictionary<int, List<Modification>> possibleVariableModifications)
        {
            Dictionary<int, Modification> modification_pattern = new Dictionary<int, Modification>();

            foreach(KeyValuePair<int, List<Modification>> kvp in possibleVariableModifications)
            {
                if(variableModificationArray[kvp.Key] > 0)
                {
                    modification_pattern.Add(kvp.Key, kvp.Value[variableModificationArray[kvp.Key] - 1]);
                }
            }

            return modification_pattern;
        }
    }
}