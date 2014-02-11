/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.IO;

namespace PeptidAce
{
    /// <summary>
    /// Reads a fasta file
    /// </summary>
    public static class ProteinFastaReader
    {
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