/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System.Collections.Generic;
using System;

namespace PeptidAce
{
	public class PeptideHit
	{
		public double[] masses;
		public int startIndex;
		public int firstQueryIndex;
		public Protein protein;
		public double mass;
		public int missedCleavage;
		public PeptideHit (Protein p, int proteinIndexStart, double[] aaMasses, double pMass, int index, int pMissedCleavage)
		{
			masses = aaMasses;
			protein = p;
			startIndex = proteinIndexStart;
			firstQueryIndex = index;
			mass = pMass;
			missedCleavage = pMissedCleavage;
		}
	}

    public class stPeptide
    {
        public string Sequence;
        public List<string> Proteins;
    }

    /// <summary>
    /// Peptide class, mostly as its defined in Morpheus, with a few twists
    /// </summary>
    public class Peptide : AminoAcidPolymer
    {
        public Protein Parent { get; set; }

        public int StartResidueNumber { get; set; }

        public int EndResidueNumber { get; set; }

        public int MissedCleavages { get; set; }

        public char PreviousAminoAcid { get; set; }

        public char NextAminoAcid { get; set; }

        public bool Decoy{get; set;}
        /*{
            get { return Parent.Decoy; }
            private set { }
        }//*/

        public bool Target
        {
            get { return !Decoy; }
            private set { Decoy = !value; }
        }

        public string ExtendedSequence()
        {
                return PreviousAminoAcid.ToString() + '.' + Sequence + '.' + NextAminoAcid.ToString();
        }

        public string ExtendedLeucineSequence()
        {
                return ExtendedSequence().Replace('I', 'L');
        }

        public string ModificationString()
        {
            string cumul = "";
            foreach (int key in VariableModifications.Keys)
                cumul += key + "[" + VariableModifications[key].Description + "]";
            return cumul;
        }

        public bool IsSamePeptide(Peptide peptide, bool checkMods)
        {
            if (BaseSequence.Equals(peptide.BaseSequence))
            {
                if (checkMods)
                {
                    if (VariableModifications != null && peptide.VariableModifications != null)
                        foreach (int key in VariableModifications.Keys)
                            if (!peptide.VariableModifications.ContainsKey(key) || VariableModifications[key] != peptide.VariableModifications[key])
                                return false;
                }
                return true;
            }
            else
                return false;
        }

        public Peptide()
        {
        }

        public static string GetSequence(string protSequence, int start, int stop)
        {
            if(start < stop)
                return protSequence.Substring(start, stop - start + 1);
            else
            {
                char[] sequence_array = protSequence.Substring(stop, start - stop + 1).ToCharArray();
                Array.Reverse(sequence_array);

                return new string(sequence_array);
            }
        }

        public Peptide Reverse()
        {
            return new Peptide(Parent, EndResidueNumber-1, StartResidueNumber-1, MissedCleavages);
        }

        public Peptide(Protein parent, int startResidueNumber, int endResidueNumber, int missedCleavages)
            : base(GetSequence(parent.BaseSequence, startResidueNumber, endResidueNumber))
        {
            Parent = parent;
            Decoy = parent.Decoy;

            StartResidueNumber = startResidueNumber + 1;
            EndResidueNumber = endResidueNumber + 1;
            MissedCleavages = missedCleavages;
            if (startResidueNumber < endResidueNumber)
            {
                if (startResidueNumber - 1 >= 0)
                    PreviousAminoAcid = parent[startResidueNumber - 1];
                else
                    PreviousAminoAcid = '-';
                if (endResidueNumber + 1 < parent.Length)
                    NextAminoAcid = parent[endResidueNumber + 1];
                else
                    NextAminoAcid = '-';
            }
            else
            {
                //Reverse sequences are decoy
                Decoy = true;
                if (endResidueNumber - 1 >= 0)
                    PreviousAminoAcid = parent[endResidueNumber - 1];
                else
                    PreviousAminoAcid = '-';
                if (startResidueNumber + 1 < parent.Length)
                    NextAminoAcid = parent[startResidueNumber + 1];
                else
                    NextAminoAcid = '-';
            }
        }

		public Peptide(PeptideHit hit)
			: base(GetSequence(hit.protein.BaseSequence, hit.startIndex, hit.startIndex + hit.masses.Length))
		{
			Parent = hit.protein;
			Decoy = Parent.Decoy;

			StartResidueNumber = hit.startIndex + 1;
			EndResidueNumber = hit.startIndex + 1 + hit.masses.Length;
			MissedCleavages = hit.missedCleavage;
			if (StartResidueNumber < EndResidueNumber)
			{
				if (StartResidueNumber - 1 >= 0)
					PreviousAminoAcid = Parent[StartResidueNumber - 1];
				else
					PreviousAminoAcid = '-';
				if (EndResidueNumber + 1 < Parent.Length)
					NextAminoAcid = Parent[EndResidueNumber + 1];
				else
					NextAminoAcid = '-';
			}
			else
			{
				//Reverse sequences are decoy
				Decoy = true;
				if (EndResidueNumber - 1 >= 0)
					PreviousAminoAcid = Parent[EndResidueNumber - 1];
				else
					PreviousAminoAcid = '-';
				if (StartResidueNumber + 1 < Parent.Length)
					NextAminoAcid = Parent[StartResidueNumber + 1];
				else
					NextAminoAcid = '-';
			}
			//Add modifications

		}

        public double[] GetMasses()
        {
            double cumul = Utilities.Constants.WATER_MONOISOTOPIC_MASS;
            double[] array = new double[BaseSequence.Length];
            for (int r = 1; r <= BaseSequence.Length; r++)
            {
                double tmp = 0;
                if (FixedModifications != null && FixedModifications.ContainsKey(r + 1))
                    foreach (Modification mod in FixedModifications[r + 1])
                        tmp += mod.MonoisotopicMassShift;
                if (VariableModifications != null && VariableModifications.ContainsKey(r + 1))
                    tmp += VariableModifications[r + 1].MonoisotopicMassShift;

                array[r - 1] = AminoAcidMasses.GetMonoisotopicMass(BaseSequence[r - 1]) + tmp;
                cumul += array[r - 1];
            }
            return array;
        }

        private Peptide(Peptide peptide) : this(peptide.Parent, peptide.StartResidueNumber-1, peptide.EndResidueNumber-1, peptide.MissedCleavages) { }

        public IEnumerable<Peptide> GetVariablyModifiedPeptides(IEnumerable<Modification> variableModifications, int maximumVariableModificationIsoforms)
        {
            Dictionary<int, List<Modification>> possible_modifications = new Dictionary<int, List<Modification>>(Length + 4);

            foreach(Modification variable_modification in variableModifications)
            {
                if(variable_modification.Type == ModificationType.ProteinNTerminus && (StartResidueNumber == 1 || (StartResidueNumber == 2 && Parent[0] == 'M')))
                {
                    List<Modification> prot_n_term_variable_mods;
                    if(!possible_modifications.TryGetValue(0, out prot_n_term_variable_mods))
                    {
                        prot_n_term_variable_mods = new List<Modification>();
                        prot_n_term_variable_mods.Add(variable_modification);
                        possible_modifications.Add(0, prot_n_term_variable_mods);
                    }
                    else
                    {
                        prot_n_term_variable_mods.Add(variable_modification);
                    }
                }

                if(variable_modification.Type == ModificationType.PeptideNTerminus)
                {
                    List<Modification> pep_n_term_variable_mods;
                    if(!possible_modifications.TryGetValue(1, out pep_n_term_variable_mods))
                    {
                        pep_n_term_variable_mods = new List<Modification>();
                        pep_n_term_variable_mods.Add(variable_modification);
                        possible_modifications.Add(1, pep_n_term_variable_mods);
                    }
                    else
                    {
                        pep_n_term_variable_mods.Add(variable_modification);
                    }
                }

                for(int r = 0; r < Length; r++)
                {
                    if(variable_modification.Type == ModificationType.AminoAcidResidue && this[r] == variable_modification.AminoAcid)
                    {
                        List<Modification> residue_variable_mods;
                        if(!possible_modifications.TryGetValue(r + 2, out residue_variable_mods))
                        {
                            residue_variable_mods = new List<Modification>();
                            residue_variable_mods.Add(variable_modification);
                            possible_modifications.Add(r + 2, residue_variable_mods);
                        }
                        else
                        {
                            residue_variable_mods.Add(variable_modification);
                        }
                    }
                }

                if(variable_modification.Type == ModificationType.PeptideCTerminus)
                {
                    List<Modification> pep_c_term_variable_mods;
                    if(!possible_modifications.TryGetValue(Length + 2, out pep_c_term_variable_mods))
                    {
                        pep_c_term_variable_mods = new List<Modification>();
                        pep_c_term_variable_mods.Add(variable_modification);
                        possible_modifications.Add(Length + 2, pep_c_term_variable_mods);
                    }
                    else
                    {
                        pep_c_term_variable_mods.Add(variable_modification);
                    }
                }

                if(variable_modification.Type == ModificationType.ProteinCTerminus && (EndResidueNumber == Parent.Length - 1))
                {
                    List<Modification> prot_c_term_variable_mods;
                    if(!possible_modifications.TryGetValue(Length + 3, out prot_c_term_variable_mods))
                    {
                        prot_c_term_variable_mods = new List<Modification>();
                        prot_c_term_variable_mods.Add(variable_modification);
                        possible_modifications.Add(Length + 3, prot_c_term_variable_mods);
                    }
                    else
                    {
                        prot_c_term_variable_mods.Add(variable_modification);
                    }
                }
            }

            int variable_modification_isoforms = 0;
            foreach(Dictionary<int, Modification> kvp in GetVariableModificationPatterns(possible_modifications))
            {
                Peptide peptide = new Peptide(this);
                peptide.SetFixedModifications(FixedModifications);
                peptide.SetVariableModifications(kvp);
                yield return peptide;
                variable_modification_isoforms++;
                if(variable_modification_isoforms >= maximumVariableModificationIsoforms)
                {
                    yield break;
                }
            }
        }

        public IEnumerable<Peptide> GetSNPsdPeptides()
        {
            string strCumul = "";
            for (int i = 0; i < BaseSequence.Length; i++)
            {
                foreach (char letter in AminoAcidMasses.AMINO_ACIDS)
                {
                    if (BaseSequence[i] != letter)
                    {
                        Peptide newPeptide = new Peptide(this);
                        if (i + 1 < BaseSequence.Length)
                            newPeptide.BaseSequence = strCumul + letter + BaseSequence.Substring(i + 1);
                        else
                            newPeptide.BaseSequence = strCumul + letter;
                        yield return newPeptide;
                    }
                }
            }
        }
    }
}