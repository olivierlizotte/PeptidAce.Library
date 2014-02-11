using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using PeptidAce;

namespace PeptidAce.Utilities.Fasta
{

    public enum ProteinIdType { Forward, Reverse, Unknown };

    public static class FastaRead
    {
        public static string Reverse(string str)
        {
            char[] sequence_array = str.ToCharArray();
            Array.Reverse(sequence_array);
            return new string(sequence_array);
        }

        public static IEnumerable<string[]> GetSequences(string fastaFileIn)
        {
            FileStream fs;
            try
            {
                fs = new FileStream(fastaFileIn, FileMode.Open, FileAccess.Read, FileShare.Read);
            }
            catch (System.Exception)
            {
                fs = new FileStream(fastaFileIn, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
            }
                
            using (StreamReader sr = new StreamReader(fs))
            {
                string line;
                string header = null;
                string sequence = "";
                while ((line = sr.ReadLine()) != null)
                {
                    if (line.StartsWith(">"))
                    {
                        if(!string.IsNullOrEmpty(header))
                            yield return new string[]{header, sequence};
                        header = line;
                        sequence = "";
                    }
                    else
                        sequence += line;
                }
                if (!string.IsNullOrEmpty(header))
                    yield return new string[] { header, sequence };
            }
            fs.Close();
        }
    }
}
