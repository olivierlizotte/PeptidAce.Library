/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using PeptidAce.Utilities;
using PeptidAce.Utilities.Methods;

namespace PeptidAce
{
    /// <summary>
    /// list of Groups of similar protein sequences
    /// </summary>
    public class ProteinGroupMatches : List<ProteinGroupMatch>
    {
        public ProteinGroupMatches(IEnumerable<ProteinGroupMatch> list) : base(list) { }
        public ProteinGroupMatches() { }
        public List<ProteinGroupMatch> ComputeAtFDR(double desired_fdr)
        {
            this.Sort(CompareProbabilityScore);
            int index = FDRizer.Extract(this, desired_fdr);
            if (index >= 0)
                return this.GetRange(0, index + 1);
            else
                return new List<ProteinGroupMatch>();
        }
        public static int CompareProbabilityScore(ProteinGroupMatch left, ProteinGroupMatch right)
        {
            return -left.ProbabilityScore().CompareTo(right.ProbabilityScore());
        }
        public static int ComparePeptideMatchesCount(ProteinGroupMatch left, ProteinGroupMatch right)
        {
            return -left.PeptideMatches.Count.CompareTo(right.PeptideMatches.Count);
        }
    }

    /// <summary>
    /// Group of similar protein sequences, that this experiment cannot differentiate
    /// </summary>
    public class ProteinGroupMatch : ITargetDecoy
    {
        public List<Protein> Proteins;
        public ProteinGroupMatch()
            : base()
        {
            PeptideMatches = new List<PeptideMatch>();
            Proteins = new List<Protein>();
        }

        public bool Target { get { return !Decoy; } }

        public bool Decoy
        {
            get
            {
                bool onlyDecoy = true;
                foreach (PeptideMatch match in PeptideMatches)//Is this method adequate? SummedScore comparison?
                    if (match.Target)
                        onlyDecoy = false;
                return onlyDecoy;
            }
        }

        public List<PeptideMatch> PeptideMatches { get; set; }
        
        public HashSet<string> BaseLeucinePeptideSequences
        {
            get
            {
                HashSet<string> base_leucine_peptide_sequences = new HashSet<string>();

                foreach (PeptideMatch match in PeptideMatches)
                    base_leucine_peptide_sequences.Add(match.peptide.BaseSequence);

                return base_leucine_peptide_sequences;
            }
        }
                
        public double ProbabilityScore()
        {
            double score = 0;
            
            foreach (PeptideMatch match in PeptideMatches)  // need option to score based on all PSMs rather than unique peptide PSMs?
                score += /*(1 - score) * */match.ProbabilityScore();
            return score;
        }
        
        public static int DescendingProbabilityScore(ProteinGroupMatch left, ProteinGroupMatch right)
        {
            int comparison = -(left.ProbabilityScore().CompareTo(right.ProbabilityScore()));
            if (comparison != 0)
                return comparison;
            else
                return left.Target.CompareTo(right.Target);
        }

        public static int AscendingProbabilityScore(ProteinGroupMatch left, ProteinGroupMatch right)
        {
            return left.ProbabilityScore().CompareTo(right.ProbabilityScore());
        }

        public static readonly string Header = "Protein Description,Protein Sequence,Protein Length,Number of Proteins in Group,Number of Peptide-Spectrum Matches,Protein Sequence Coverage (%),Summed Morpheus Score,Decoy";

        public override string ToString()
        {
            StringBuilder description = new StringBuilder();
            StringBuilder sequence = new StringBuilder();
            StringBuilder length = new StringBuilder();
            if (Proteins.Count > 0)
            {
                foreach (Protein protein in this.Proteins)
                {
                    description.Append(protein.Description.Replace(',', ' ') + " / ");
                    sequence.Append(protein.Sequence + '/');
                    length.Append(protein.Sequence.Length.ToString() + '/');
                }
                description = description.Remove(description.Length - 3, 3);
                sequence = sequence.Remove(sequence.Length - 1, 1);
                length = length.Remove(length.Length - 1, 1);
            }
            StringBuilder sb = new StringBuilder();

            sb.Append(description.ToString() + ',');
            sb.Append(sequence.ToString() + ',');
            sb.Append(length.ToString() + ',');
            sb.Append(Proteins.Count.ToString() + ',');
            sb.Append(PeptideMatches.Count.ToString() + ',');
            StringBuilder sequence_coverage = new StringBuilder();
            if (Proteins.Count > 0)
            {
                foreach (Protein protein in this.Proteins)
                    sequence_coverage.Append((CalculateSequenceCoverage(PeptideMatches, protein) * 100.0).ToString() + '/');
                sequence_coverage = sequence_coverage.Remove(sequence_coverage.Length - 1, 1);
            }
            sb.Append(sequence_coverage.ToString() + ',');
            sb.Append(ProbabilityScore().ToString() + ',');
            sb.Append(Decoy.ToString());

            return sb.ToString();
        }

        public double CalculateSequenceCoverage(List<PeptideMatch> matches, Protein protein)
        {
            HashSet<int> covered_residues = new HashSet<int>();
            foreach (PeptideMatch match in matches)
            {
                for (int r = match.peptide.StartResidueNumber; r <= match.peptide.EndResidueNumber; r++)
                {
                    covered_residues.Add(r);
                }
            }
            return (double)covered_residues.Count / protein.Length;
        }
    }

    /// <summary>
    /// Methods to search for protein discovered, based on the list of peptides and clusters found
    /// TODO aggregate peptides into proteins based on the quantification (profile)
    /// </summary>
    public class ProteinSearcher
    {
        DBOptions options;

        public ProteinSearcher(DBOptions options)
        {
            this.options = options;
        }

        public ProteinGroupMatches SearchLatest(List<PeptideMatch> peptides)//, IDictionary<string, List<Protein>> dicOfProteins)//, List<Protein> AllProteins)// Dictionary<string, List<Protein>> dicOfPeptides)
        {
            Dictionary<string, PeptideMatch> dicOfPeptideMatches = new Dictionary<string,PeptideMatch>();
            foreach(PeptideMatch match in peptides)
                dicOfPeptideMatches.Add(match.peptide.BaseSequence, match);

            Dictionary<Protein, List<string>> allPossibleProteins = new Dictionary<Protein, List<string>>();
            foreach (PeptideMatch match in peptides)
            //foreach (string peptideSequence in dicOfProteins.Keys)
            {
                Protein protein = match.peptide.Parent;

                //foreach (Protein protein in dicOfProteins[peptideSequence])
                {
                    if (!allPossibleProteins.ContainsKey(protein))
                        allPossibleProteins.Add(protein, new List<string>());
                    allPossibleProteins[protein].Add(match.peptide.Sequence);
                }
            }

            options.ConSole.WriteLine("Merging undistinguishable proteins...");

            foreach(Protein protein in allPossibleProteins.Keys)
                allPossibleProteins[protein].Sort();

            //Merge only undistinguishable protein sequences
            //Dictionary<Protein, List<Protein>> dicOfSimilarProteins = new Dictionary<Protein,List<Protein>>();
            ProteinGroupMatches proteinMatches = new ProteinGroupMatches();
            List<Protein> proteins = new List<Protein>(allPossibleProteins.Keys);
            for(int i = 0; i < proteins.Count; i++)
            {
                ProteinGroupMatch groupMatch = new ProteinGroupMatch();
                groupMatch.Proteins.Add(proteins[i]);
                foreach(string sequence in allPossibleProteins[proteins[i]])
                    if(dicOfPeptideMatches.ContainsKey(sequence))
                        groupMatch.PeptideMatches.Add(dicOfPeptideMatches[sequence]);

                for(int j = i + 1; j < proteins.Count;)
                {
                    bool isSame = false;
                    if(allPossibleProteins[proteins[i]].Count == allPossibleProteins[proteins[j]].Count)
                    {
                        isSame = true;
                        for(int k = 0; k < allPossibleProteins[proteins[i]].Count; k++)
                        {
                            if(allPossibleProteins[proteins[i]][k].CompareTo(allPossibleProteins[proteins[j]][k]) != 0)
                                isSame = false;
                        }
                    }
                    if(isSame)
                    {
                        proteins.RemoveAt(j);
                        groupMatch.Proteins.Add(proteins[j]);
                    }
                    else
                        j++;
                }
                proteinMatches.Add(groupMatch);
            }
            options.ConSole.WriteLine("Found " + proteinMatches.Count + " protein groups.");
            return proteinMatches;
        }

        public static void Export(string filename, List<ProteinGroupMatch> proteins)
        {
            vsCSVWriter writer = new vsCSVWriter(filename);
            writer.AddLine(ProteinGroupMatch.Header);
            foreach (ProteinGroupMatch group in proteins)
                writer.AddLine(group.ToString());
            writer.WriteToFile();
        }

        public static IEnumerable<Peptide> ProteinDigest(DBOptions options, List<Protein> Proteins, bool allowSNP)
        {
            int nbProteins = 0;
            bool tooLong = false;
            foreach (Protein protein in Proteins)
            {
                nbProteins++;
                if (tooLong)
                    break;
                List<int> indices = options.DigestionEnzyme.GetDigestionSiteIndices(protein);
                indices.Insert(0, -1);
                indices.Add(protein.Length - 1);

                for (int missed_cleavages = 0; missed_cleavages <= options.ToleratedMissedCleavages; missed_cleavages++)
                {
                    for (int i = 0; i < indices.Count - missed_cleavages - 1; i++)
                    {
                        if (indices[i + missed_cleavages + 1] + 1 - (indices[i] + 1 + 1) + 1 >= options.MinimumPeptideLength)
                        {
                            if (options.initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || indices[i] + 1 != 0 || protein[0] != 'M')
                            {
                                Peptide newPep = new Peptide(protein, indices[i] + 1, indices[i + missed_cleavages + 1], missed_cleavages);
                                yield return newPep;

                                if (allowSNP)
                                    foreach (Peptide possibleSnp in newPep.GetSNPsdPeptides())
                                        yield return possibleSnp;
                            }

                            if (options.initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && indices[i] + 1 == 0 && protein[0] == 'M')
                            {
                                Peptide newPep = new Peptide(protein, indices[i] + 1 + 1, indices[i + missed_cleavages + 1], missed_cleavages);
                                yield return newPep;

                                if (allowSNP)
                                    foreach (Peptide possibleSnp in newPep.GetSNPsdPeptides())
                                        yield return possibleSnp;
                            }
                        }
                    }
                }
                Console.Write("\r{0}%   ", ((100 * nbProteins) / Proteins.Count));
            }
            Console.Write("\r{0}%   ", 100);
        }
    }
}