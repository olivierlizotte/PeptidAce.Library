using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using PeptidAce.Utilities;

namespace CMinor.Custom
{
    public static class ERKMotifFindor
    {
        public static void Find(string fastaFile, string peptidesFile)
        {
// TO Put in db one day!
//DEF motif list
//$def_listarray("F.[FY]P"); // paper ERK
 
// D motif list
//$d_list = array("[KR]..[KR].{1,5}[LI].[LI]", "[KR]{2,5}.{1,6}[LIV].[LIV]"); // paper ERK ;

            Dictionary<string, Regex> dicOfMotifs = new Dictionary<string, Regex>();
            dicOfMotifs.Add("Def", new Regex("F.[FY]P"));
            dicOfMotifs.Add("D 1", new Regex("[KR]..[KR].{1,5}[LI].[LI]"));
            dicOfMotifs.Add("D 2", new Regex("[KR]{2,5}.{1,6}[LIV].[LIV]"));

            vsCSV csvPep = new vsCSV(peptidesFile);
            int indexPepSeq = csvPep.GetColumnIndex("Peptide Sequence");
            Dictionary<string, List<string[]>> dicOfPeptides = new Dictionary<string, List<string[]>>();
            for(int i = 1; i < csvPep.LINES_LIST.Count; i++)
            {
                string[] splits = vsCSV.Split_HandleQuotes(csvPep.LINES_LIST[i]);
                if (!dicOfPeptides.ContainsKey(splits[indexPepSeq]))
                    dicOfPeptides.Add(splits[indexPepSeq], new List<string[]>());

                dicOfPeptides[splits[indexPepSeq]].Add(splits);
            }
            
            Dictionary<string, string> proteins = PeptidAce.ProteinFastaReader.ReadProteins(fastaFile);
            PeptidAce.Utilities.vsCSVWriter writer = new vsCSVWriter(vsCSV.GetFolder(peptidesFile) + vsCSV.GetFileName_NoExtension(peptidesFile) + "_ERK_Motif.csv");
            writer.AddLine(csvPep.LINES_LIST[0] + ",Start index in Protein,Index of Ph in protein,Motif Found,Index of Motif in protein,Motif sequence,Distance from peptide,SP Motif");
            foreach (KeyValuePair<string, string> protPair in proteins)
            {
                string protein = protPair.Value;
                foreach (string peptide in dicOfPeptides.Keys)
                {
                    foreach(string[] occurance in dicOfPeptides[peptide])
                    { 
                        if((!string.IsNullOrEmpty(occurance[0]) && protPair.Key.Contains(occurance[0])) || 
                            (!string.IsNullOrEmpty(occurance[1]) && protPair.Key.Contains(occurance[1])))
                        {
                            int index = protein.IndexOf(peptide);
                            while (index < protein.Length && index >= 0)
                            {                                
                                if (index > 0)
                                {
                                    int startIndex = index - 50;
                                    if (startIndex < 0) startIndex = 0;
                                    int length = peptide.Length + 40;
                                    if (startIndex + length >= protein.Length)
                                        length = protein.Length - startIndex;
                                    //int endIndex   = index + peptide.Length + 40;
                                    //if (endIndex >= protein.Length) endIndex = protein.Length - 1;
                                    string subProt = protein.Substring(startIndex, length);
                                    foreach (KeyValuePair<string, Regex> pair in dicOfMotifs)
                                        foreach (Match match in pair.Value.Matches(subProt))
                                        {
                                            string line = "";
                                            int otherIndexes = 0;
                                            int bestDistance = int.MaxValue;
                                            while (otherIndexes >= 0)
                                            {
                                                otherIndexes++;
                                                bool hasSP = false;
                                                bool hasPH = false;
                                                int indexOfPh = occurance[3].IndexOf("S(ph)P", otherIndexes);
                                                if (indexOfPh < 0)
                                                    indexOfPh = occurance[3].IndexOf("(ph)", otherIndexes);
                                                else
                                                    hasSP = true;
                                                otherIndexes = indexOfPh;
                                                if (indexOfPh < 0)
                                                    indexOfPh = peptide.Length / 2;
                                                else
                                                    hasPH = true;

                                                indexOfPh += index + 1;
                                                int indexOfMotif = match.Index + startIndex + 1;
                                                int distance = Math.Abs(indexOfMotif - indexOfPh);
                                                if(distance < bestDistance)
                                                {
                                                    bestDistance = distance;
                                                    line = vsCSV.GetLineFromColumns(occurance, 0) + "," + (index + 1) + "," + indexOfPh + "," + pair.Key + "," + indexOfMotif + "," + subProt.Substring(match.Index, 4) + "," + distance + "," + hasSP;
                                                }
                                            }
                                            writer.AddLine(line);
                                            
                                        }

                                } 
                                index = protein.IndexOf(peptide, index + 1);
                            }
                        }
                    }
                }
                
            }
            writer.WriteToFile();
        }
    }
}
