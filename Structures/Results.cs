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
using PeptidAce.Utilities;

namespace PeptidAce
{
    public class Result
    {
        public DBOptions dbOptions;
        public Samples samples;
        public Queries queries;
        public Precursors precursors;
        public Precursors matchedPrecursors { get; private set; }
        public List<Cluster> clusters;
        public PeptideMatches peptides;
        public PeptideMatches peptideSequences;
        public ProteinGroupMatches proteins;

        public Result()
        {
            queries = new Queries();
            precursors = new Precursors();
            matchedPrecursors = new Precursors();
            clusters = new List<Cluster>();
            peptides = new PeptideMatches(new PeptideMatch[0]);
            peptideSequences = new PeptideMatches(new PeptideMatch[0]);
            proteins = new ProteinGroupMatches();
            dbOptions = new DBOptions("");
            samples = new Samples();
        }
        /*
        public Result(Precursors precursors, List<Cluster> clusters, PeptideMatches peptides, PeptideMatches peptideSequences, ProteinGroupMatches proteins, Queries queries, DBOptions dbOptions)
        {
            this.queries = queries;
            this.precursors = precursors;
            this.SetPrecursors(precursors);
            this.clusters = clusters;
            this.peptides = peptides;
            this.peptideSequences = peptideSequences;
            this.proteins = proteins;
            this.dbOptions = dbOptions;
        }//*/

        public long SetPrecursors(Precursors precursors)
        {
            this.precursors = precursors;
            this.matchedPrecursors = new Precursors();
            long nbTargets = 0;
            foreach (Precursor precursor in precursors)
                if (precursor.psms != null && precursor.psms.Count > 0)
                {
                    if (precursor.Target)
                        nbTargets++;
                    matchedPrecursors.Add(precursor);
                }
            return nbTargets;
        }
        /*
        public List<PeptideSpectrumMatch> GetPSMs(double fdr, double decoyOverTargetRatio = 1, int nbMaxPsm = 1)
        {
            return FDR.PSMs(clusters, fdr, decoyOverTargetRatio, nbMaxPsm);
        }//*/
        /*
        public List<PeptideMatch> GetPeptides(double fdr)
        {
            if (peptides != null)
                return FDR.Peptides(peptides, fdr);
            else
                return null;
        }

        public List<ProteinGroupMatch> GetProteins(double fdr)
        {
            if (proteins != null)
                return FDR.Proteins(proteins, fdr);
            else
                return null;
        }//*/
        
        public static int DescendingProteinScoreComparison(Precursor left, Precursor right)
        {
            PeptideSpectrumMatch psmLeft = left.OptimizedBestPsm();
            if (psmLeft != null)
            {
                PeptideSpectrumMatch psmRight = right.OptimizedBestPsm();
                if (psmRight != null)
                    return -(psmLeft.ProteinScore.CompareTo(psmRight.ProteinScore));
                else
                    return -1;
            }
            else return 1;
        }

        public static int DescendingIntensityFractionComparison(Precursor left, Precursor right)
        {
            return -(left.OptimizedBestPsm().MatchingIntensityFraction.CompareTo(right.OptimizedBestPsm().MatchingIntensityFraction));
        }

        public static int CountTargets(List<ITargetDecoy> elems)
        {
            int nbTarget = 0;
            foreach (ITargetDecoy t in elems)
                if (t.Target)
                    nbTarget++;
            return nbTarget;
        }

        public void WriteInfoToCsv(bool light = false)
        {
            int target = 0;
            foreach (Precursor precursor in matchedPrecursors)
                if (precursor.Target)
                    target++;
            dbOptions.ConSole.WriteLine("  ---  Number of precursors          : " + target + " targets [" + matchedPrecursors.Count + "]" + "  ---  ");

            target = 0;
            foreach (PeptideMatch peptide in peptides)
                if (peptide.Target)
                    target++;
            dbOptions.ConSole.WriteLine("  ---  Number of peptides            : " + target + " targets [" + peptides.Count + "]" + "  ---  ");

            target = 0;
            foreach (PeptideMatch peptide in peptideSequences)
                if (peptide.Target)
                    target++;
            dbOptions.ConSole.WriteLine("  ---  Number of peptide sequences   : " + target + " targets [" + peptideSequences.Count + "]" + "  ---  ");

            target = 0;
            foreach (ProteinGroupMatch protein in proteins)
                if (protein.Target)
                    target++;
            dbOptions.ConSole.WriteLine("  ---  Number of proteins            : " + target + " targets [" + proteins.Count + "]" + "  ---  ");
            if (!light)
            {
                //WriteFragmentation(true);
                //WriteFragmentation(false);
            }
        }
        /*
        public void WriteFragmentation(bool target)
        {
            vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + "FragmentStats_" + (target ? "Targets" : "Decoy") + ".csv");
            writer.AddLine("  === Fragmentation of " + (target ? "Targets" : "Decoys") + " ===");
            foreach (FragmentClass fragment in dbOptions.fragments)
            {
                double cumulIntensity = 0;
                int nbFrag = 0;
                Dictionary<int, int> positions = new Dictionary<int, int>();
                foreach (Precursor precursor in matchedPrecursors)
                {
                    PeptideSpectrumMatch psm = precursor.OptimizedBestPsm();
                    if (psm.Target == target)
                    {
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment == match.Fragment)
                            {
                                nbFrag++;
                                if (!positions.ContainsKey(match.fragmentPos))
                                    positions.Add(match.fragmentPos, 1);
                                else
                                    positions[match.fragmentPos]++;
                                cumulIntensity += match.obsIntensity;
                            }
                    }
                }
                string strPos = "";
                if (positions.Count > 0)
                    foreach (int key in positions.Keys)
                        strPos += "|" + key + ":" + positions[key];
                else
                    strPos += ",";
                writer.AddLine("    " + fragment.Name + ", Number of fragments = , " + nbFrag + ",   Intensity = ," + cumulIntensity + ", fragment matched [" + strPos.Substring(1) + "]");
            }
            foreach (FragmentClass fragment in dbOptions.fragments)
            //foreach (string fragment in FragmentDictionary.Fragments.Keys)
            {
                double cumulIntensity = 0;
                int nbFrag = 0;
                foreach (Precursor precursor in matchedPrecursors)
                {
                    PeptideSpectrumMatch psm = precursor.OptimizedBestPsm();
                    if (psm.Target == target)
                    {
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment == match.Fragment)
                            {
                                nbFrag++;
                                cumulIntensity += match.obsIntensity;
                            }
                    }
                }
                writer.AddLine("    " + fragment + ", Number of fragments = ," + nbFrag + ",   Intensity = ," + cumulIntensity);
            }
            foreach (FragmentClass fragment in dbOptions.fragments)
            //foreach (string fragment in FragmentDictionary.AAFragments.Keys)
            {
                double cumulIntensity = 0;
                int nbFrag = 0;
                foreach (Precursor precursor in matchedPrecursors)
                {
                    PeptideSpectrumMatch psm = precursor.OptimizedBestPsm();
                    if (psm.Target == target)
                    {
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment == match.Fragment)
                            {
                                nbFrag++;
                                cumulIntensity += match.obsIntensity;
                            }
                    }
                }
                writer.AddLine("    " + fragment + ", Number of fragments = ," + nbFrag + ",   Intensity = ," + cumulIntensity);
            }
            writer.WriteToFile();
        }//*/
        /*
        public void ExportFragments(PeptideSpectrumMatch psm)
        {
            vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + psm.Peptide.Sequence + "_" + vsCSV.GetFileName_NoExtension(psm.Query.sample.sSDF) + "_" + psm.Query.precursor.Track.RT + ".csv");
            List<FragmentClass> fragments = new List<FragmentClass>();
            foreach (FragmentClass fragment in dbOptions.fragments)
            //foreach (string fragment in FragmentDictionary.Fragments.Keys)
            {
                bool found = false;
                foreach (ProductMatch match in dbOptions.fragments.ComputeFragments(psm.Peptide.GetMasses(), psm.Peptide.BaseSequence, psm.Query.precursor.Charge, dbOptions))
                {
                    if (fragment == match.Fragment)
                    {
                        found = true;
                        break;
                    }
                }
                if (found)
                    fragments.Add(fragment);
            }

            string title = "Theoretical Fragments";
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (FragmentClass fragment in dbOptions.fragments)
                    title += "," + fragment.Name + " ^" + charge;
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (FragmentClass fragment in fragments)
                    title += "," + fragment.Name + " ^" + charge;
            writer.AddLine(title);

            for (int i = 1; i <= psm.Peptide.Length; i++)
            {
                string line = i.ToString();
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (FragmentClass fragment in dbOptions.fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in dbOptions.fragments.ComputeFragments(psm.Peptide.GetMasses(), psm.Peptide.BaseSequence, psm.Query.precursor.Charge, dbOptions))
                        {
                            if (fragment == match.Fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.theoMz;
                                found = true;
                                break;
                            }
                        }
                        if (!found)
                            line += ",";
                    }
                }
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (FragmentClass fragment in fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in dbOptions.fragments.ComputeFragments(psm.Peptide.GetMasses(), psm.Peptide.BaseSequence, psm.Query.precursor.Charge, dbOptions))
                        {
                            if (fragment == match.Fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.theoMz;
                                found = true;
                                break;
                            }
                        }
                        if (!found)
                            line += ",";
                    }
                }
                writer.AddLine(line);
            }

            title = "Observed Fragments Intensities";
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (FragmentClass fragment in dbOptions.fragments)
                    title += "," + fragment.Name + " ^" + charge;
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (FragmentClass fragment in fragments)
                    title += "," + fragment.Name + " ^" + charge; 
            writer.AddLine(title);

            for(int i = 1; i <= psm.Peptide.Length; i++)
            {
                string line = i.ToString();
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (FragmentClass fragment in dbOptions.fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment == match.Fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.obsIntensity;
                                found = true;
                                break;
                            }
                        if (!found)
                            line += ",";
                    }
                }
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (FragmentClass fragment in fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment == match.Fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.obsIntensity;
                                found = true;
                                break;
                            }
                        if (!found)
                            line += ",";
                    }
                }
                writer.AddLine(line);
            }

            title = "Observed Fragments Mz";
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (FragmentClass fragment in dbOptions.fragments)
                    title += "," + fragment.Name + " ^" + charge;
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (FragmentClass fragment in fragments)
                    title += "," + fragment.Name + " ^" + charge;
            writer.AddLine(title);

            for (int i = 1; i <= psm.Peptide.Length; i++)
            {
                string line = i.ToString();
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (FragmentClass fragment in dbOptions.fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment == match.Fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.obsMz;
                                found = true;
                                break;
                            }
                        if (!found)
                            line += ",";
                    }
                }
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (FragmentClass fragment in fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment == match.Fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.obsMz;
                                found = true;
                                break;
                            }
                        if (!found)
                            line += ",";
                    }
                }
                writer.AddLine(line);
            }

            title = "Error on Fragments";
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (FragmentClass fragment in dbOptions.fragments)
                    title += "," + fragment.Name + " ^" + charge;
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (FragmentClass fragment in fragments)
                    title += "," + fragment.Name + " ^" + charge;
            writer.AddLine(title);

            for (int i = 1; i <= psm.Peptide.Length; i++)
            {
                string line = i.ToString();
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (FragmentClass fragment in dbOptions.fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment == match.Fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.mass_diff;
                                found = true;
                                break;
                            }
                        if (!found)
                            line += ",";
                    }
                }
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (FragmentClass fragment in fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment == match.Fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.mass_diff;
                                found = true;
                                break;
                            }
                        if (!found)
                            line += ",";
                    }
                }
                writer.AddLine(line);
            }
            writer.WriteToFile();
        }//*/
        /*
        public void ExportFragmentIntensities(List<PeptideSpectrumMatch> psms, Peptide peptide, int psmCharge, string fileName)
        {
            vsCSVWriter writer = new vsCSVWriter(fileName);
            List<FragmentClass> fragments = new List<FragmentClass>();
            foreach (FragmentClass fragment in dbOptions.fragments)
            //foreach (string fragment in FragmentDictionary.Fragments.Keys)
            {
                bool found = false;
                foreach (ProductMatch match in dbOptions.fragments.ComputeFragments(peptide.GetMasses(), peptide.BaseSequence, psmCharge, dbOptions))
                {
                    if (fragment == match.Fragment)
                    {
                        found = true;
                        break;
                    }
                }
                if (found)
                    fragments.Add(fragment);
            }

            string title = "Cumulated Product Intensities";
            for (int charge = 1; charge <= psmCharge; charge++)
                foreach (FragmentClass fragment in dbOptions.fragments)
                    title += "," + fragment.Name + " ^" + charge;
            for (int charge = 1; charge <= psmCharge; charge++)
                foreach (FragmentClass fragment in fragments)
                    title += "," + fragment.Name + " ^" + charge;
            writer.AddLine(title);

            for (int i = 1; i <= peptide.Length; i++)
            {
                string line = i.ToString();
                for (int charge = 1; charge <= psmCharge; charge++)
                {
                    foreach (FragmentClass fragment in dbOptions.fragments)
                    {
                        double cumul = 0.0;
                        foreach(PeptideSpectrumMatch psm in psms)
                            foreach (ProductMatch match in psm.AllProductMatches)
                                if (fragment == match.Fragment && match.fragmentPos == i && match.charge == charge)
                                    cumul += match.obsIntensity;
                        line += "," + cumul;
                    }
                }
                for (int charge = 1; charge <= psmCharge; charge++)
                {
                    foreach (FragmentClass fragment in fragments)
                    {
                        double cumul = 0.0;
                        foreach (PeptideSpectrumMatch psm in psms)
                            foreach (ProductMatch match in psm.AllProductMatches)
                                if (fragment == match.Fragment && match.fragmentPos == i && match.charge == charge)
                                    cumul += match.obsIntensity;
                        line += "," + cumul;
                    }
                }
                writer.AddLine(line);
            }
            writer.WriteToFile();
        }
        //*/
        /*
        public void ExportFragmentIntensitiesForAllPSM(List<PeptideSpectrumMatch> psms, Peptide peptide, int psmCharge, string fileName)
        {
            vsCSVWriter writer = new vsCSVWriter(fileName);
            string title = "Retention Time";
            for (int i = 1; i <= peptide.Length; i++)
                for (int charge = 1; charge <= psmCharge; charge++)
                    foreach (FragmentClass fragment in dbOptions.fragments)
                        title += "," + i + fragment.Name + " ^" + charge;
            writer.AddLine(title);

            foreach (PeptideSpectrumMatch psm in psms)
            {
                string line = psm.Query.spectrum.RetentionTimeInMin.ToString();
                for (int i = 1; i <= peptide.Length; i++)
                {
                    for (int charge = 1; charge <= psmCharge; charge++)
                    {
                        foreach (FragmentClass fragment in dbOptions.fragments)
                        {
                            double cumul = 0.0;
                            foreach (ProductMatch match in psm.AllProductMatches)
                                if (fragment == match.Fragment && match.fragmentPos == i && match.charge == charge)
                                    cumul += match.obsIntensity;
                            line += "," + cumul;
                        }
                    }
                }
                writer.AddLine(line);
            }
            writer.WriteToFile();
        }//*/

        public void Export(double fdr, string keyword = "", bool onlyPrecursors = false)
        {
            dbOptions.ConSole.WriteLine("Exporting at " + (fdr * 100) + "% FDR (Decoy/Target)...");

            List<Precursor> prec = null;
            if (precursors != null)
            {
                if (matchedPrecursors == null)
                {
                    this.matchedPrecursors = new Precursors();
                    foreach (Precursor precursor in precursors)
                        if (precursor.psms.Count > 0)
                            matchedPrecursors.Add(precursor);
                }/*
                List<Precursor> prec = FDR.PrecursorsV2(precursors, fdr, 1);
                Sol.CONSOLE.OutputLine(">   " + prec.Count + " Precursors");
                MSSearcher.Export(dbOptions.outputFolder + keyword + "precursors.csv", prec);//*/
                /*
                prec = Optimizer.PrecursorOptimizer(matchedPrecursors, fdr);
                Sol.CONSOLE.OutputLine(">   " + prec.Count + " Optimized Precursors");
                MSSearcher.Export(dbOptions.outputFolder + keyword + "Optimized_precursors.csv", prec);//*/

                prec = matchedPrecursors.ComputeAtFDR(fdr);
                dbOptions.ConSole.WriteLine(">   " + prec.Count + " Uptimized V5 Precursors");
                MSSearcher.Export(dbOptions.OutputFolder + keyword + "UptimizedV5_precursors.csv", prec);
            }

            if (!onlyPrecursors)
            {
                if (queries != null)
                {
                    List<Query> qs = queries.ComputeAtFDR(fdr);
                    dbOptions.ConSole.WriteLine(">   " + qs.Count + " PSMs (Top 10)");
                    MSSearcher.Export(dbOptions.OutputFolder + keyword + "queries.csv", qs);
                }
                if (clusters != null)
                {
                    //List<PeptideSpectrumMatch> psms = FDR.PSMs(clusters, fdr, 1, 10);
                    //dbOptions.ConSole.WriteLine(">   " + psms.Count + " PSMs (Top 10)");
                    //MSSearcher.Export(dbOptions.outputFolder + keyword + "psms_Top10.csv", psms);

                    //psms = FDR.PSMs(clusters, fdr, 1, 1);
                    //dbOptions.ConSole.WriteLine(">   " + psms.Count + " PSMs");
                    //MSSearcher.Export(dbOptions.outputFolder + keyword + "psms_Best.csv", psms);
                }

                if (peptides != null)
                {
                    List<PeptideMatch> pep = peptideSequences.ComputeAtFDR(fdr);
                    dbOptions.ConSole.WriteLine(">   " + pep.Count + " Peptides Sequences (Version 5)");
                    PeptideSearcher.Export(dbOptions.OutputFolder + keyword + "peptideSequencesV5_.csv", pep);

                    PeptideSearcher sr = new PeptideSearcher(dbOptions);
                    PeptideMatches seqs = sr.Search(clusters, prec, false);
                    pep = seqs.ComputeAtFDR(fdr);
                    dbOptions.ConSole.WriteLine(">   " + pep.Count + " Peptides Sequences (Version 5b)");
                    PeptideSearcher.Export(dbOptions.OutputFolder + keyword + "peptideSequencesV5b_PrecursorFDRed.csv", pep);

                    pep = peptides.ComputeAtFDR(fdr);
                    dbOptions.ConSole.WriteLine(">   " + pep.Count + " Peptides (Version 5)");
                    PeptideSearcher.Export(dbOptions.OutputFolder + keyword + "peptidesV5_.csv", pep);

                }
                if (proteins != null)
                {
                    List<ProteinGroupMatch> prots = proteins.ComputeAtFDR(fdr);
                    dbOptions.ConSole.WriteLine(">   " + prots.Count + " Proteins");
                    ProteinSearcher.Export(dbOptions.OutputFolder + keyword + "proteins_.csv", prots);
                }
            }
        }

        public void ExportPSMs(int nbPsmPerQuery, string filename)
        {
            //List<Query> qs = queries.ComputeAtFDR(1);
            MSSearcher.Export(filename, queries);
        }
    }    
}
