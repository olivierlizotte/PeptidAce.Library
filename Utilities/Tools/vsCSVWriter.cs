/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;


// This class WILL BE used to create pre-visualized-and-then-saved CSV files
namespace PeptidAce.Utilities
{
    /// <summary>
    /// Eases exporting text to comma separated value files. Just don't forget to call 
    /// WriteToFile, as everything is stored in memory until you call that method!
    /// </summary>
    public class vsCSVWriter
    {
        public static char[] _Generic_Separator = { ',' };
        
        public void Replace(string oldString, string newString)
        {
            for (int i = 0; i < m_strLines.Count; i++)
                m_strLines[i] = m_strLines[i].Replace(oldString, newString);
        }
        
        public long CountLength()
        {
            long cptChar = 0;
            foreach (string str in m_strLines)
                cptChar += str.Length + 2;
            return cptChar;
        }

        public void AddLine(string line)
        {
            m_strLines.Add(line);
        }

        public void AddLines(string[] lines)
        {
            m_strLines.AddRange(lines);
        }

        public void InsertLine(string line, int index)
        {
            m_strLines.Insert(index, line);
        }
        public string m_fileName;
        public List<string> m_strLines;
        public string m_strTitle = "";

        //Constructor still available if no visual is required
        public vsCSVWriter(string fileOutput)
        {
            m_strLines = new List<string>();
            m_fileName = fileOutput;
        }
        
        public bool WriteToFile()
        {
            int fileTries = 20;
            while (fileTries > 0)
            {
                try
                {
                    if (!System.IO.Directory.Exists(vsCSV.GetFolder(m_fileName)))
                        System.IO.Directory.CreateDirectory(vsCSV.GetFolder(m_fileName));

                    using (StreamWriter sw = new StreamWriter(m_fileName))
                    {
                        if (!string.IsNullOrEmpty(m_strTitle))
                            sw.WriteLine(m_strTitle);
                        foreach (string line in m_strLines)
                            sw.WriteLine(line);
                        sw.Close();
                    }
                    return true;
                }
                catch (System.Exception)
                {
                    m_fileName = vsCSV.GetFolder(m_fileName) + vsCSV.GetFileName_NoExtension(m_fileName) + "_" + (new Random()).Next() + ".csv";                    
                    fileTries--;
                }
            }

            return false;
        }
    }
}
