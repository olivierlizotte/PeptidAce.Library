﻿/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using PeptidAce.Utilities;

namespace PeptidAce
{
    /// <summary>
    /// Propheus is the main search engine for PeptidAce (previously known as Trinity). At its most basic, 
    /// it is a modification of the Morpheus 
    /// peptide identification scheme, tailored for no enzyme searches
    /// </summary>
    public class Ace
    {
        public static Result Start(DBOptions dbOptions, string fasta, Samples project, bool loadMS1 = true, bool filterMS2 = true, Queries queries = null)
        {
            Ace pAce = new Ace(dbOptions, project);
            pAce.Preload(fasta, loadMS1, filterMS2);
            if (queries == null)
                pAce.PrepareQueries();
            else
                pAce.AllQueries = queries;

            if(pAce.AllQueries.Count > 0)
                return pAce.LaunchSearch(pAce.AllQueries);
            return null;
        }

        public DBOptions dbOptions;
        public Samples Project;

        //Preloaded lists
        public List<Protein> AllProteins;
        public Dictionary<Sample, Spectra> AllSpectras;
        public Queries AllQueries;
        
        public Ace(DBOptions dbOptions, Samples Project)
        {
            this.dbOptions = dbOptions;
            this.Project = Project; 
        }

        /// <summary>
        /// Builds the list of proteins to digest, the spectrum to match them to
        /// Also creates the list of queries (some spectrum are used more than once, when multiple 
        /// precursors are found in the mass range of the fragmentation window
        /// </summary>
        public void Preload(string fasta, bool loadMS1, bool filterMS2 = true)
        {
            AllProteins             = ReadProteomeFromFasta(fasta, !dbOptions.DecoyFusion);
            AllSpectras             = LoadSpectras(loadMS1, filterMS2);
        }

        /// <summary>
        /// Builds the list of proteins to digest, the spectrum to match them to
        /// Also creates the list of queries (some spectrum are used more than once, when multiple 
        /// precursors are found in the mass range of the fragmentation window
        /// </summary>
        public void PrepareQueries()
        {
            AllQueries = CreateQueries(AllSpectras, dbOptions);
        }

        /// <summary>
        /// Load a Propheus searh object from a Result object
        /// Can be used to load a save state
        /// </summary>
        /// <param name="tmp"></param>
        public void Load(Result tmp, string fasta)
        {
            AllProteins = ReadProteomeFromFasta(fasta, !dbOptions.DecoyFusion);
            AllQueries = tmp.queries;
            
            AllSpectras = new Dictionary<Sample, Spectra>();
            foreach(Query query in tmp.queries)
            {
                if(!AllSpectras.ContainsKey(query.sample))
                    AllSpectras.Add(query.sample, null);
            }
        }

        /// <summary>
        /// Reads the proteins from a fasta file format
        /// </summary>
        /// <param name="fileName"></param>
        /// <param name="onTheFlyDecoys"> Adds reverse sequnce proteins </param>
        /// <returns></returns>
        public static List<Protein> ReadProteomeFromFasta(string fileName, bool addReverseProteins)
        {
            List<Protein> AllProteins = new List<Protein>();
            try
            {
                //dbOptions.ConSole.WriteLine("Reading FASTA file " + fileName + " ... ");
                //Extract Proteins from Fasta file
                FileStream protein_fasta_database = new FileStream(fileName, FileMode.Open, FileAccess.Read, FileShare.Read);
                foreach (Protein protein in ProteinFastaReader.ReadProteins(protein_fasta_database, addReverseProteins))
                {
                    AllProteins.Add(protein);
                }
                //AllProteins.Sort(Protein.TargetDecoyComparison);
                protein_fasta_database.Close();
                //dbOptions.ConSole.WriteLine("Proteins in fasta file : " + AllProteins.Count);
            }
            catch(Exception e)
            {
                //dbOptions.ConSole.WriteLine("Error reading fasta file : " + fileName);
                Console.WriteLine(e.StackTrace);
            }
            return AllProteins;
        }        
        
        public static Spectra LoadSpectra(string file, DBOptions dbOptions, bool loadMS = true, bool filterMS2 = true)
        {
            string trackFile = vsCSV.GetFolder(file) + vsCSV.GetFileName_NoExtension(file) + "_Tracks.csv";
            string msmsIonFile = vsCSV.GetFolder(file) + vsCSV.GetFileName_NoExtension(file) + "_MSMSIons.csv";
            if (dbOptions.LoadSpectraIfFound && System.IO.File.Exists(trackFile)
                                            && System.IO.File.Exists(msmsIonFile))
            {
                dbOptions.ConSole.WriteLine("Loading Sectra from " + trackFile + " AND " + msmsIonFile);
                if (loadMS)
                    return Spectra.Import(msmsIonFile, trackFile, dbOptions);
                else
                    return Spectra.Import(msmsIonFile, null, dbOptions);
            }
            else
            {
                dbOptions.ConSole.WriteLine("Loading Sectra " + file);

                pwiz.CLI.msdata.MSDataFile msFile = new pwiz.CLI.msdata.MSDataFile(file);
                Spectra spectra = Spectra.Load(msFile, dbOptions, file, loadMS, filterMS2);
                spectra.Sort(ProductSpectrum.AscendingPrecursorMassComparison);

                dbOptions.ConSole.WriteLine(file + " [" + spectra.Count + " msms scans]");
                if (dbOptions.SaveMS1Peaks)
                    spectra.ExportTracks(trackFile);

                if (dbOptions.SaveMSMSPeaks)
                    spectra.ExportMSMS(msmsIonFile);

                return spectra;
            }
        }

        /// <summary>
        /// Loads spectra from Raw files
        /// </summary>
        /// <returns></returns>
        public Dictionary<Sample, Spectra> LoadSpectras(bool loadMS = true, bool filterMS2 = true)
        {
            //TODO test compatibility with mzML, .d folders ... other known formats compatible with ProteoWizard
            AllSpectras = new Dictionary<Sample, Spectra>();
            for (int i = 0; i < Project.Count; i++)
            {
                Sample sample = Project[i];
                AllSpectras.Add(sample, LoadSpectra(sample.sSDF, dbOptions, loadMS, filterMS2));
            }
            List<Sample> keys = new List<Sample>(AllSpectras.Keys);
            foreach (Sample s in keys)
            {
                Spectra toKeep = new Spectra();
                foreach (ProductSpectrum ps in AllSpectras[s])
                    if (ps.PrecursorIntensity > 0)
                        toKeep.Add(ps);
                AllSpectras[s] = toKeep;
            }
            return AllSpectras;            
        }
        
        /// <summary>
        /// Creates the list of queries
        /// </summary>
        /// <param name="spectras"></param>
        /// <returns></returns>
        public static Queries CreateQueries(Dictionary<Sample, Spectra> spectras, DBOptions dbOptions)
        {
            Queries AllQueries = new Queries(dbOptions);
            foreach(Sample entry in spectras.Keys)
                AllQueries.GenerateQueries(entry, spectras[entry], spectras[entry].tracks);
                        
            return AllQueries;
        }

        /// <summary>
        /// Launches ProPheus to match queries to peptide sequences
        /// </summary>
        /// <param name="queries"></param>
        /// <returns></returns>
        public Result LaunchSearch(Queries queries)
        {
            Result result = new Result();
            result.queries = queries;
            result.dbOptions = dbOptions;
            result.samples = Project;

            ProPheus dbSearcher = new ProPheus(dbOptions);
            Digestion ps = new Digestion(dbOptions);

            result.SetPrecursors(dbSearcher.Search(queries, ps.DigestProteomeOnTheFly(AllProteins, false, queries)));

            dbOptions.ConSole.WriteLine(result.precursors.Count + " precursors matched !");
            
            foreach (Precursor precursor in result.precursors)
                foreach (PeptideSpectrumMatch psm in precursor.psms_AllPossibilities)
                        precursor.psms.Add(psm);

            long nbTargets = result.SetPrecursors(result.precursors);
            dbOptions.ConSole.WriteLine("Targets : " + nbTargets + " [" + result.matchedPrecursors.Count + "]");
            dbOptions.ConSole.WriteLine(result.matchedPrecursors.Count + " precursors remaining after ProPheus Search!");
            return result;
        }
    }
}
