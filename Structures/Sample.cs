/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using PeptidAce.Utilities;

namespace PeptidAce
{
    public class Samples : List<Sample>
    {
        public string FileName = "";
        public DBOptions dbOptions;
        public Samples()
        {
        }

        public Samples(DBOptions options)
        {
            dbOptions = options;
        }

        public Samples(string projectFileName, int maxFractionSpreading, DBOptions options)
        {
            loadProjectFile(projectFileName, maxFractionSpreading);
            Sort(Compare);
            FileName = projectFileName;
            dbOptions = options;
        }
        /*
        public static vsSDF LoadSDF(Sample project)
        {
            try
            {
                //if (project.m_vsSDF == null || string.Compare(project.m_vsSDF.m_name, project.sSDF) != 0)
                //{
                    vsSDF result = null;
                    result = vsSDF.Load(project.sSDF, null);
                    result.MZ_GAPS = Proteomics.Cluster.Scoring.MzTol(result.MAX_MZ_VALUE);
                    project.sSDF = result.m_name;//Change the name in case there was misspelling fixed by the SDF.Load method
                    //project.m_vsSDF = result;
                //}
                    return result;// project.m_vsSDF;
            }
            catch (Exception ex)
            {
                Sol.CONSOLE.OutputLine("Error jk302me :" + ex.Message + "    " + ex.StackTrace);
                return null;
            }
        }//*/
        
        private void loadProjectFile(string project_file, int maxFractionSpreading)
        {
            // current line number in the project csv file
            int lineNum = 0;
            //vsList<stCondition> iterativeHelper = null;
            try
            {
                this.Clear();

                // project csv file is comma delimited
                char[] splitter = { ',' };
                
                vsCSV csvProject = new vsCSV(project_file);
                for (int i = 0; i < csvProject.LINES_LIST.Count; i++)
                {
                    try
                    {
                        string line = csvProject.LINES_LIST[i];
                        lineNum++;

                        string[] lineParts = line.Split(splitter);

                        // skip header line
                        if (lineNum == 1)
                        {
                            continue;
                        }

                        // project csv file format: REPLICATE,REPLICATE,FRACTION,peptide map location, sdf file location
                        // ex. 1,1,1,c:\work\adhoc\Promix_070706\Promix_10_070706_Out_peptides.csv,c:\work\adhoc\Promix_070706\Promix_10_070706_Out.sdf
                        int condition = int.Parse(lineParts[0]);
                        int replicate = int.Parse(lineParts[1]);
                        int fraction = int.Parse(lineParts[2]);
                        string peptideMap = lineParts[3];
                        string sdf = lineParts[4];
                        string nameCol = "";
                        if (lineParts.Length >= 6)
                            nameCol = lineParts[5];

                        // KE Aug 7 2008 - lower case
                        //OLI Why lower case? Files does not exist under linux!
                        this.Add(new Sample(condition, replicate, fraction, peptideMap, sdf, maxFractionSpreading, nameCol));
                    }
                    catch (System.Exception ex)
                    {
                        dbOptions.ConSole.WriteLine("Error sfj4aau34 : LineNum = " + lineNum + "   " + ex.Message + "    " + ex.StackTrace);
                    }
                }
            }
            catch (Exception ex)
            {
                dbOptions.ConSole.WriteLine("Error sfj4u34 : " + ex.Message + "  \n  " + ex.StackTrace);
            }
            //ITERATIVE_HELPER = iterativeHelper;
//            GenerateIterativeHelper();
  //          return iterativeHelper;
        }

        private static int Compare(Sample objA, Sample objB)
        {
            int compare = objA.PROJECT.CONDITION.CompareTo(objB.PROJECT.CONDITION);
            if (compare == 0)
            {
                compare = objA.PROJECT.REPLICATE.CompareTo(objB.PROJECT.REPLICATE);
                if (compare == 0)
                    compare = objA.PROJECT.FRACTION.CompareTo(objB.PROJECT.FRACTION);
            }
            return compare;
        }
    }

    public class Sample
    {
        //public static object m_userControl = null;
        public SampleEntry PROJECT;
		public string sPeptideMap;
		public string sSDF;

        public string Name
        {
            get {
                if (string.IsNullOrEmpty(nameColumn))
                    return vsCSV.GetFileName_NoExtension(sSDF);
                else
                    return nameColumn;
            }
        }
        //private vsSDF m_vsSDF;
        //public vsSDF m_vsSDF;
        /*public vsSDF GetSDF()
        {
            if (m_vsSDF == null)
                m_vsSDF = Samples.LoadSDF(this);
            return m_vsSDF;
        }//*/
        public string nameColumn = "";

        public Sample()
        {
        }

        public Sample(int ConditionNum, int ReplicateNum, int FractionNum, string PeptideMap, string SDF, int maxFractionSpreading, string nameCol)
		{
            this.PROJECT = new SampleEntry(ConditionNum, ReplicateNum, FractionNum, maxFractionSpreading);
			this.sPeptideMap = PeptideMap;
			this.sSDF = SDF;
            //this.m_vsSDF = null;
            this.nameColumn = nameCol;
		}
    }

    public class SampleEntry
    {    
        public int CONDITION = -1;
        public int REPLICATE = -1;
        public int FRACTION = -1;
        public int MIN_FRACTION = -65536;
        public int MAX_FRACTION = 65536;
        public long ID_PROJECT_ENTRY = 0;

        public SampleEntry()
        {
        }

        public SampleEntry(int condition, int replicate, int fraction, int ecart)
        {
            CONDITION = condition;
            REPLICATE = replicate;
            FRACTION = fraction;
            if (FRACTION != -1)
            {
                MIN_FRACTION = (int)(FRACTION - ecart);
                MAX_FRACTION = (int)(FRACTION + ecart);
            }
        }
    }
}
