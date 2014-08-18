/*
 * Wenger CD, Coon JJ. A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra, Journal of Proteome Research, 2013; 12(3): 1377-86
 * http://www.chem.wisc.edu/~coon/software.php#morpheus
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 * Altered by Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 */
using System;
using System.Collections.Generic;
using PeptidAce.Utilities;

namespace PeptidAce
{
    public class Protein : AminoAcidPolymer
    {
        public string Description { get; set; }

        public Protein()
        {
            Decoy = false;
        }

        public Protein(string sequence, string description, bool decoy)
            : base(sequence, true)
        {
            this.Description = description;
            this.Decoy = decoy;// Description.Contains("DECOY");
        }

        public bool Decoy;
        /*
        public static int TargetDecoyComparison(Protein left, Protein right)
        {
            if (left.Target)
                if (right.Target)
                    return 0;
                else
                    return -1;
            else
                if (right.Target)
                    return 1;
                else
                    return 0;
        }//*/

        public IEnumerable<Peptide> Digest(Protease protease, int maximumMissedCleavages, InitiatorMethionineBehavior initiatorMethionineBehavior,
            int? minimumPeptideLength, int? maximumPeptideLength, bool onTheFlyDecoy)
        {
            List<int> indices = protease.GetDigestionSiteIndices(this);
            indices.Insert(0, -1);
            indices.Add(Length - 1);

            for(int missed_cleavages = 0; missed_cleavages <= maximumMissedCleavages; missed_cleavages++)
            {
                for(int i = 0; i < indices.Count - missed_cleavages - 1; i++)
                {
                    if(initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || indices[i] + 1 != 0 || this[0] != 'M')
                    {
                        Peptide peptide = new Peptide(this, indices[i] + 1, indices[i + missed_cleavages + 1], missed_cleavages);

                        if((!minimumPeptideLength.HasValue || peptide.Length >= minimumPeptideLength.Value)
                            && (!maximumPeptideLength.HasValue || peptide.Length <= maximumPeptideLength.Value))
                        {
                            yield return peptide;
                            
                            if(onTheFlyDecoy)
                                yield return new Peptide(this, indices[i + missed_cleavages + 1], indices[i] + 1, missed_cleavages);
                        }
                    }

                    if(initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && indices[i] + 1 == 0 && this[0] == 'M')
                    {
                        Peptide peptide_without_initiator_methionine = new Peptide(this, indices[i] + 1 + 1, indices[i + missed_cleavages + 1], missed_cleavages);

                        if((!minimumPeptideLength.HasValue || peptide_without_initiator_methionine.Length >= minimumPeptideLength.Value)
                            && (!maximumPeptideLength.HasValue || peptide_without_initiator_methionine.Length <= maximumPeptideLength.Value))
                        {
                            yield return peptide_without_initiator_methionine;

                            if(onTheFlyDecoy)
                                yield return new Peptide(this, indices[i + missed_cleavages + 1], indices[i] + 1 + 1, missed_cleavages);
                        }
                    }
                }
            }
        }

        public IEnumerable<Protein> GetVariablyModifiedProteins(IEnumerable<Modification> variableModifications, int maximumVariableModificationIsoforms)
        {
            Dictionary<int, List<Modification>> possible_modifications = new Dictionary<int, List<Modification>>();

            foreach(Modification variable_modification in variableModifications)
            {
                if(variable_modification.Type == ModificationType.ProteinNTerminus)
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

                if(variable_modification.Type == ModificationType.ProteinCTerminus)
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
            foreach(Dictionary<int, Modification> kvp in GetVariableModificationPatterns(possible_modifications, maximumVariableModificationIsoforms))
            {
                Protein protein = new Protein(BaseSequence, Description, Decoy);
                protein.SetFixedModifications(FixedModifications);
                protein.SetVariableModifications(kvp);
                yield return protein;
                variable_modification_isoforms++;
                if(variable_modification_isoforms == maximumVariableModificationIsoforms)
                {
                    yield break;
                }
            }
        }

        public double CalculateSequenceCoverage(IEnumerable<PeptideSpectrumMatch> peptideSpectrumMatches)
        {
            HashSet<int> covered_residues = new HashSet<int>();
            foreach(PeptideSpectrumMatch psm in peptideSpectrumMatches)
            {
                for(int r = psm.Peptide.StartResidueNumber; r <= psm.Peptide.EndResidueNumber; r++)
                {
                    covered_residues.Add(r);
                }
            }
            return (double)covered_residues.Count / Length;
        }
    }
}