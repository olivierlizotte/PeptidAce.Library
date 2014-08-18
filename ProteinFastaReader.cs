/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace PeptidAce
{
    public class FastProtein
    {
        public char[][] Sequence;
        public Modification[][][] Modifications;
        public double[][] Masses;
        public FastProtein(char[][] sequence, double[][] masses, Modification[][][] modifications)
        {
            Sequence = sequence;
            Masses = masses;
            Modifications = modifications;
        }
    }

    public class FastGene
    {
        public FastProtein[] Proteins;
        public FastGene(string superSequence, DBOptions options)
        {
            List<FastProtein> proteins = new List<FastProtein>();
            string[] splits = superSequence.Split('*');
            foreach (string subSeq in splits)
            {
                List<List<char>> list = new List<List<char>>();
                int indexSeq = 0;
                while (indexSeq < subSeq.Length)
                {
                    if (subSeq[indexSeq] == '/')
                    {
                        indexSeq++;
                        list[list.Count - 1].Add(subSeq[indexSeq]);
                    }
                    else list.Add(new List<char>() { subSeq[indexSeq] });
                    indexSeq++;
                }
                char[][] array = new char[list.Count][];
                double[][] arrayMass = new double[list.Count][];
                Modification[][][] arrayMod = new Modification[list.Count][][];
                for (int i = 0; i < list.Count; i++)
                {
                    array[i] = list[i].ToArray();
                    arrayMass[i] = new double[array[i].Length];
                    arrayMod[i] = new Modification[array[i].Length][];
                    for (int j = 0; j < array[i].Length; j++)
                    {
                        arrayMass[i][j] = AminoAcidMasses.GetMonoisotopicMass(array[i][j]);
                        foreach (Modification mod in options.fixedModifications)
                            if (mod.AminoAcid == array[i][j])
                                arrayMass[i][j] += mod.MonoisotopicMassShift;

                        List<Modification> mods = new List<Modification>();
                        foreach (Modification mod in options.variableModifications)
                            if (mod.AminoAcid == array[i][j])
                                mods.Add(mod);
                        if (mods.Count > 0)
                            arrayMod[i][j] = mods.ToArray();
                    }
                }
                proteins.Add(new FastProtein(array, arrayMass, arrayMod));
            }
            Proteins = proteins.ToArray();
        }
    }

    /// <summary>
    /// Reads a fasta file
    /// </summary>
    public static class ProteinFastaReader
    {
        
        private static List<char[][]> BuildProteinAAArray(string superSequence)
        {
            List<char[][]> results = new List<char[][]>();
            string[] splits = superSequence.Split('*');
            foreach (string sequence in splits)
            {
                List<List<char>> list = new List<List<char>>();
                int indexSeq = 0;
                while (indexSeq < sequence.Length)
                {
                    if (sequence[indexSeq] == '/')
                    {
                        indexSeq++;
                        list[list.Count - 1].Add(sequence[indexSeq]);
                    }
                    else list.Add(new List<char>() { sequence[indexSeq] });
                    indexSeq++;
                }
                char[][] array = new char[list.Count][];
                for (int i = 0; i < list.Count; i++)
                    array[i] = list[i].ToArray();
                results.Add(array);
            }
            return results;
        }//*/        

        
        private static List<char[][]> BuildProteinAAArrayWithMods(string superSequence)
        {
            List<char[][]> results = new List<char[][]>();
            string[] splits = superSequence.Split('*');
            foreach (string sequence in splits)
            {
                List<List<char>> list = new List<List<char>>();
                int indexSeq = 0;
                while (indexSeq < sequence.Length)
                {
                    if (sequence[indexSeq] == '/')
                    {
                        indexSeq++;
                        list[list.Count - 1].Add(sequence[indexSeq]);
                    }
                    else list.Add(new List<char>() { sequence[indexSeq] });
                    indexSeq++;
                }
                char[][] array = new char[list.Count][];
                for (int i = 0; i < list.Count; i++)
                    array[i] = list[i].ToArray();
                results.Add(array);
            }
            return results;
        }//*/

        public static Dictionary<string, FastGene> ReadProteinChars(string superFastaFile, DBOptions options)
        {
            //superFasta format:
            //>protein_id|other info|incluing Gene symbol
            //ASEQUENCEWITHALTERNATE/SAAFORSomePOSITIONSAND*FORSTOPCODONS
            //This version does not support shifts (indels and inserts)

            //Dic of proteins <Description,ListOf*EndingSequences<[AA Position,AA possibilities]>>
            //
            Dictionary<string, FastGene> genes = new Dictionary<string, FastGene>();
            FileStream proteinFastaDatabase = new FileStream(superFastaFile, FileMode.Open, FileAccess.Read, FileShare.Read);
            using (StreamReader fasta = new StreamReader(proteinFastaDatabase))
            {
                string description = null;
                string sequence = "";

                while (true)
                {
                    string line = fasta.ReadLine();

                    if (line.StartsWith(">"))
                        description = line.Substring(1);
                    else
                        sequence += line.Trim();

                    if (fasta.Peek() == '>' || fasta.Peek() == -1)
                    {
                        if (!string.IsNullOrEmpty(description) && sequence.Length > 0)     
                            genes.Add(description, new FastGene(sequence, options));

                        description = null;
                        sequence = "";

                        if (fasta.Peek() == -1)
                        {
                            break; ;
                        }
                    }
                }
                fasta.Close();
            }
            proteinFastaDatabase.Close();
            return genes;
        }
        /*
        public static Dictionary<string, List<char[][]>> ReadProteinChars(string superFastaFile)
        {
            //superFasta format:
            //>protein_id|other info|incluing Gene symbol
            //ASEQUENCEWITHALTERNATE/SAAFORSomePOSITIONSAND*FORSTOPCODONS
            //This version does not support shifts (indels and inserts)

            //Dic of proteins <Description,ListOf*EndingSequences<[AA Position,AA possibilities]>>
            //
            Dictionary<string, List<char[][]>> proteins = new Dictionary<string, List<char[][]>>();
            FileStream proteinFastaDatabase = new FileStream(superFastaFile, FileMode.Open, FileAccess.Read, FileShare.Read);
            using (StreamReader fasta = new StreamReader(proteinFastaDatabase))
            {
                string description = null;
                string sequence = "";

                while (true)
                {
                    string line = fasta.ReadLine();

                    if (line.StartsWith(">"))
                        description = line.Substring(1);
                    else
                        sequence += line.Trim();

                    if (fasta.Peek() == '>' || fasta.Peek() == -1)
                    {
                        if (!string.IsNullOrEmpty(description) && sequence.Length > 0)                         
                            proteins.Add(description, BuildProteinAAArray(sequence));

                        description = null;
                        sequence = "";

                        if (fasta.Peek() == -1)
                        {
                            break; ;
                        }
                    }
                }
                fasta.Close();
            }
            proteinFastaDatabase.Close();
            return proteins;
        }
        //*/
        public static Dictionary<string, string> ReadProteins(string[] fastaFiles)
        {
            Dictionary<string, string> allProts = new Dictionary<string, string>();

            foreach (string fasta in fastaFiles)
                foreach (KeyValuePair<string, string> pair in ReadProteins(fasta))
                    allProts.Add(pair.Key, pair.Value);
            return allProts;
        }

        public static Dictionary<string, string> ReadProteins(string fastafile)
        {
            Dictionary<string, string> proteins = new Dictionary<string, string>();
            FileStream proteinFastaDatabase = new FileStream(fastafile, FileMode.Open, FileAccess.Read, FileShare.Read);
            using (StreamReader fasta = new StreamReader(proteinFastaDatabase))
            {

                string description = null;
                string sequence = null;

                while (true)
                {
                    string line = fasta.ReadLine();

                    if (line.StartsWith(">"))
                        description = line.Substring(1);
                    else
                        sequence += line.Trim();

                    if (fasta.Peek() == '>' || fasta.Peek() == -1)
                    {
                        if (!string.IsNullOrEmpty(description) && !string.IsNullOrEmpty(sequence))
                            proteins.Add(description, sequence);

                        description = null;
                        sequence = null;

                        if (fasta.Peek() == -1)
                        {
                            break; ;
                        }
                    }
                }
                fasta.Close();
            }
            proteinFastaDatabase.Close();
            return proteins;
        }

        public static IEnumerable<Protein> ReadProteins(FileStream proteinFastaDatabase, bool onTheFlyDecoys)
        {
            StreamReader fasta = new StreamReader(proteinFastaDatabase);

            string description = null;
            string sequence = null;

            while(true)
            {
                string line = fasta.ReadLine();

                if(line.StartsWith(">"))
                {
                    description = line.Substring(1);
                }
                else
                {
                    sequence += line.Trim();
                }

                if(fasta.Peek() == '>' || fasta.Peek() == -1)
                {
                    Protein protein = new Protein(sequence, description, false);
                    
                    yield return protein;

                    if(onTheFlyDecoys)
                    {
                        if(protein.Decoy)
                        {
                            throw new ArgumentException(proteinFastaDatabase.Name + " contains decoy proteins; database should not contain decoy proteins when \"create target–decoy database on the fly\" option is enabled");
                        }
                    
                        char[] sequence_array = sequence.ToCharArray();
                        if(sequence.StartsWith("M"))
                        {
                            Array.Reverse(sequence_array, 1, sequence.Length - 1);
                        }
                        else
                        {
                            Array.Reverse(sequence_array);
                        }
                        string reversed_sequence = new string(sequence_array);
                        Protein decoy_protein = new Protein(reversed_sequence, (description.Length > 2 && description[2] == '|') ? description.Insert(3, "DECOY_") : "DECOY_" + description, true);
                        yield return decoy_protein;//*/
                    }

                    description = null;
                    sequence = null;

                    if(fasta.Peek() == -1)
                    {
                        break;
                    }
                }
            }

            proteinFastaDatabase.Seek(0, SeekOrigin.Begin);
        }
    }
}