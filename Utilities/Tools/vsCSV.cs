/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace PeptidAce.Utilities
{
    /// <summary>
    /// This class makes reading comma separated value files a little less painful 
    /// </summary>
    public class vsCSV
    {
        public static char[] _Generic_Separator = { ',' };
        
        public List<string> LINES_LIST
        {
            get { return m_strLines; }
        }

        public static void Export(string filename, params List<double>[] values)
        {
            vsCSVWriter writer = new vsCSVWriter(filename);
            bool keepGoing = true;
            int lineNumber = 0;
            while(keepGoing)
            {
                keepGoing = false;
                string line = "";
                for (int i = 0; i < values.Length; i++)
                    if (lineNumber < values[i].Count)
                    {
                        keepGoing = true;
                        line += values[i][lineNumber] + ",";
                    }
                    else
                        line += ",";

                lineNumber++;
                writer.AddLine(line);
            }
            writer.WriteToFile();
        }

        private static int Count(string text, string pattern)
        {
            // Loop through all instances of the string 'text'.
            int count = 0;
            int i = 0;
            while ((i = text.IndexOf(pattern, i)) != -1)
            {
                i += pattern.Length;
                count++;
            }
            return count;
        }
        
        public static string[] Split_HandleQuotes(string line, long nbColMax = int.MaxValue)
        {
            char[] separator = { ',' };
            string[] array = line.Split(_Generic_Separator);
            string[] columns = new string[nbColMax];
            bool inQuotes = false;
            int iter = -1;
            for(int i = 0; i < array.Length && iter + 1 < nbColMax; i++)
            {
                if (inQuotes)
                    columns[iter] += "," + array[i];
                else
                {
                    iter++;
                    columns[iter] = array[i];
                }

                if (Count(array[i], "\"") % 2 != 0)
                    inQuotes = !inQuotes;
            }
            return columns;
        }
        
        //Necessary variables 
        protected string m_strFileName = "";
        protected List<String> m_strLines;
        
        public String m_strFirstLine = "";
        protected int m_validLineNumber = 0;

        public string m_name;

        //Constructor with 6 values
        public vsCSV(string csvFileName)
        {
            m_strFileName = csvFileName;
            m_name = getNameFromFilename(m_strFileName);
            loadCSVFile(csvFileName);
        }
        
        //Retrieves the name of the file by removing the folder names
        public static String getNameFromFilename(String fileName)
        {
            char[] splitter = { '\\' };
            string[] splits = fileName.Split(splitter);
            return splits[splits.Length - 1];
        }

        //Retrieves the first line of the csv file (usually the header)
        public String getFirstLine()
        {
            return m_strFirstLine;
        }

        //Open CSV File 
        protected void loadCSVFile(string CSVFileName)
        {
            int line_number = 0;
            try
            {
                m_strLines = ReadCSVFile(CSVFileName);

                if (m_strLines != null)
                {
                    line_number = m_strLines.Count;
                    if (m_strLines.Count > 0)
                        m_strFirstLine = m_strLines[0];
                }
            }
            catch (Exception e)
            {
                throw new Exception("Error in  loadCSVFile! " + e.Message + " at line " + line_number);
            }
        }

        public int GetColumnIndex(string title)
        {
            if (string.IsNullOrEmpty(m_strFirstLine))
                m_strFirstLine = m_strLines[0];
            string[] splits = m_strFirstLine.Split(_Generic_Separator);
            for (int i = 0; i < splits.Length; i++)
            {
                if (string.Compare(title, splits[i]) == 0)
                    return i;
            }
            return -1;
        }

        public int GetColumnIndex(string title, int nbLine)
        {
            string[] splits = m_strLines[nbLine].Split(_Generic_Separator);
            for (int i = 0; i < splits.Length; i++)
            {
                if (string.Compare(title, splits[i]) == 0)
                    return i;
            }
            return -1;
        }
        
        public static List<string> ReadCSVFile(string csvFile)
        {
            List<string> strLines = new List<string>();
            try
            {
                FileStream fs;
                try
                {
                    fs = new FileStream(csvFile, FileMode.Open, FileAccess.Read, FileShare.Read);
                }
                catch (System.Exception)
                {
                    fs = new FileStream(csvFile, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
                }
                using (StreamReader sr = new StreamReader(fs))
                {
                    string line;

                    while ((line = sr.ReadLine()) != null)
                    {
                        strLines.Add(line);
                    }
                }
                fs.Close();
            }
            catch (System.Exception ex)
            {
                //Well, the file didn't open.... what to do...?
                Console.WriteLine("Could not open file :" + ex.Message);
                //strLines = null;
            }
            return strLines;
        }

        public static int GetColumnIndex(List<string> titleLine, string columnTitle)
        {
            //Lets try to find it
            for (int i = 0; i < titleLine.Count; i++)
                if (titleLine[i].Equals(columnTitle))
                    return i;

            //Could not find exact match, try lower case
            string lowerCase = columnTitle.ToLower();
            for (int i = 0; i < titleLine.Count; i++)
                if (titleLine[i].ToLower().Equals(lowerCase))
                    return i;

            return -1;
        }

        public static int GetColumnIndex(string[] titleLine, string columnTitle)
        {
            //Lets try to find it
            for (int i = 0; i < titleLine.Length; i++)
                if (titleLine[i].Equals(columnTitle))
                    return i;

            //Could not find exact match, try lower case
            string lowerCase = columnTitle.ToLower();
            for (int i = 0; i < titleLine.Length; i++)
                if (titleLine[i].ToLower().Equals(lowerCase))
                    return i;

            return -1;
        }

        public static List<string> GetColumn(List<string> lines, int index, string skip, char separator)
        {
            List<string> sdfs = new List<string>();
            
            foreach (string line in lines)
            {
                string[] splits = line.Split(separator);
                if (string.Compare(splits[index], skip, true) != 0)//skip this pattern
                    sdfs.Add(splits[index]);
            }
            return sdfs;
        }
        
        public static string GetFolder(string fileName)
        {
            int index = fileName.LastIndexOf("\\");
            if(index == -1)
                index = fileName.LastIndexOf("/");
            return fileName.Substring(0, index + 1);
        }

        public static string GetFileName(string folderAndFileName)
        {
            int index = folderAndFileName.LastIndexOf("\\");
            if (index == -1)
                index = folderAndFileName.LastIndexOf("/");
            if (index != -1)
                return folderAndFileName.Substring(index + 1, folderAndFileName.Length - (index + 1));
            else
                return folderAndFileName;
        }

        public static string GetFileName_NoExtension(string folderAndFilename)
        {
            string result = GetFileName(folderAndFilename);
            int index = result.LastIndexOf(".");
            if (index != -1)
                return result.Remove(index);
            else
                return result;
        }

        public static string GetLineFromColumns(string[] columns, int indexStart)
        {
            return GetLineFromColumns(columns, indexStart, columns.Length);
        }

        public static string GetLineFromColumns(string[] columns, int indexStart, int maxCloumns)
        {
            string result = "";
            if (columns.Length > indexStart)
            {
                result = columns[indexStart];
                if (maxCloumns > columns.Length)
                    maxCloumns = columns.Length;

                for (int i = indexStart + 1; i < maxCloumns; i++)
                    result += "," + columns[i];
            }
            return result;
        }

        public static string Concatenate(string[] splits, string separator)
        {
            string results = splits[0];
            for (int i = 1; i < splits.Length; i++)
                results += separator + splits[i];
            return results;
        }
        
        public static bool StartsWithNumber(string line)
        {
            if (!string.IsNullOrEmpty(line) && line.Length > 0)
                return line[0] >= '0' && line[0] <= '9';
            
            return false;
        }

        public static bool StartsWithNumberOrLetter(string line)
        {
            if (!string.IsNullOrEmpty(line) && line.Length > 0)
                return (line[0] >= '0' && line[0] <= '9') || (line[0] >= 'A' && line[0] <= 'Z') || (line[0] >= 'a' && line[0] <= 'z');

            return false;
        }
    }
}
