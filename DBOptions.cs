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
using PeptidAce.Utilities.Interfaces;

namespace PeptidAce
{
    /// <summary>
    /// Stores options related to a sample file or group of files
    /// Values are adjusted during calculations to reflec detected thresholds
    /// </summary>
    public class DBOptions
    {
        public string DBName;

        public int MinimumPrecursorChargeState = 1;
        public int MaximumPrecursorChargeState = 8;
        public int MaximumNumberOfFragmentsPerSpectrum = 100;
        public double MaximumPeptideMass;
        public bool DecoyFusion;

        public Protease DigestionEnzyme;
        public int ToleratedMissedCleavages;

        public double PSMFalseDiscoveryRate = 0.0504;

        public int NbPSMToKeep = 8;
        public string OutputFolder;
        public string TempFolder;
        public double MinimumPSMScore;
        public int MinimumPeptideLength = 5;
        public int MaximumPeptideLength = 100;

        public FullFragments fullFragment = new FullFragments(new char[]{'b','y'}, false, false);//false, false, false);
        public FullFragments fragmentErrorFragment = new FullFragments(new char[] { 'b', 'y' }, false, false);
        public FullFragments everyFragments = new FullFragments(new char[] { 'a', 'b', 'c', 'x', 'y', 'z' }, true, true);

        //public Fragments fragments;
        public int NbMinProducts = 4;

        public bool addFragmentLoss = false;
        public bool addFragmentMods = false;
        public bool WriteMaxQuantPeakFile = false;

        public bool SaveMS1Peaks = true;
        public bool SaveMSMSPeaks = true;

        public bool LoadSpectraIfFound = true;

        public double ComputedRetentionTimeDiff = 1.0;//TODO compute this after alignment step, based on common identifications
        public double EffectiveIsolationWindowRatio = 0.9;//TODO optimize this value 
        public double MinimumPrecursorIntensityRatioInIsolationWindow = 0.05;
        
        public double dProduct = 0.0;
        public double dPrecursor = 0.12;
        public double dMatchingProductFraction = 0.45;
        public double dMatchingProduct = 0;
        public double dIntensityFraction = 0.13;
        public double dIntensity = 0;
        public double dProtein = 0;
        public double dPeptideScore = 0.3;
        public double dFragmentScore = 0.0;

        public Dictionary<string, double> dicOfMinPrecursorThr;
        public Dictionary<string, double> dicOfMaxPrecursorThr;
        public Dictionary<string, double> dicOfMinFragmentThr;
        public Dictionary<string, double> dicOfMaxFragmentThr;

        ///Values kept from the original Morpheus source code
        public InitiatorMethionineBehavior initiatorMethionineBehavior;
        public List<Modification> fixedModifications;
        public List<Modification> variableModifications;
        public int maximumVariableModificationIsoforms;
        public MassTolerance precursorMassTolerance = new MassTolerance(0.005, MassToleranceUnits.Da);
        public MassTolerance productMassTolerance = new MassTolerance(0.005, MassToleranceUnits.Da);

        public string[] FastaFiles;
        public string[] FastaFilesMinors;
        public string[] MS2s;

        public IConSol ConSole;
        public DBOptions(IConSol console = null)
        {
            if (console == null)
                ConSole = new ConSolCommandLine();
            else
                ConSole = console;
            //Create with default values
            this.DecoyFusion = true;
            this.MaximumPeptideMass = 10000;
            ProteaseDictionary proteases = ProteaseDictionary.Instance;
            this.DigestionEnzyme = proteases["no enzyme"];// proteases["trypsin (no proline rule)"];
            this.ToleratedMissedCleavages = 100;// 3;//determines the length of peptides with no-enzyme option
            this.initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            this.fixedModifications = new List<Modification>();
            this.variableModifications = new List<Modification>();
            this.maximumVariableModificationIsoforms = 1024;
            
            this.MinimumPrecursorChargeState = 1;
            this.MaximumPrecursorChargeState = 4;
            this.MaximumNumberOfFragmentsPerSpectrum = 400;
            //TODO Add precision to the precursor by reading MS part of file
            this.precursorMassTolerance = new MassTolerance(0.05, MassToleranceUnits.Da);//2.1
            //TODO Add precision to the product masses by reading corresponding MS part of raw file
            this.productMassTolerance = new MassTolerance(0.05, MassToleranceUnits.Da);

            this.PSMFalseDiscoveryRate = 0.0504;

            //this.OutputFolder = @"C:\_IRIC\DATA\Test2";//C:\Documents and Settings\ProteoAdmin\Desktop\AEffacer\Morpheus\Output";
            this.MinimumPSMScore = 0.0001;
            this.TempFolder = System.IO.Path.GetTempPath();
        }

        public DBOptions Clone()
        {
            DBOptions clone = new DBOptions(ConSole);
            clone.DBName = DBName;
            clone.MinimumPrecursorChargeState = MinimumPrecursorChargeState;
            clone.MaximumPrecursorChargeState = MaximumPrecursorChargeState;
            clone.MaximumNumberOfFragmentsPerSpectrum = MaximumNumberOfFragmentsPerSpectrum;
            clone.MaximumPeptideMass = MaximumPeptideMass;
            clone.DecoyFusion = DecoyFusion;
            clone.DigestionEnzyme = DigestionEnzyme;
            clone.ToleratedMissedCleavages = ToleratedMissedCleavages;
            clone.PSMFalseDiscoveryRate = PSMFalseDiscoveryRate;

            clone.NbPSMToKeep = NbPSMToKeep;
            clone.OutputFolder = OutputFolder;
            clone.TempFolder = TempFolder;
            clone.MinimumPSMScore = MinimumPSMScore;
            clone.MinimumPeptideLength = MinimumPeptideLength;
            clone.MaximumPeptideLength = MaximumPeptideLength;
            
            clone.fullFragment = fullFragment;
            clone.everyFragments = everyFragments;

            clone.NbMinProducts = NbMinProducts;

            clone.addFragmentLoss = addFragmentLoss;
            clone.addFragmentMods = addFragmentMods;
            clone.WriteMaxQuantPeakFile = WriteMaxQuantPeakFile;

            clone.SaveMS1Peaks = SaveMS1Peaks;
            clone.SaveMSMSPeaks = SaveMSMSPeaks;

            clone.LoadSpectraIfFound = LoadSpectraIfFound;

            clone.ComputedRetentionTimeDiff = ComputedRetentionTimeDiff;
            clone.EffectiveIsolationWindowRatio = EffectiveIsolationWindowRatio;
            clone.MinimumPrecursorIntensityRatioInIsolationWindow = MinimumPrecursorIntensityRatioInIsolationWindow;


            clone.dProduct = dProduct;
            clone.dPrecursor = dPrecursor;
            clone.dMatchingProductFraction = dMatchingProductFraction;
            clone.dMatchingProduct = dMatchingProduct;
            clone.dIntensityFraction = dIntensityFraction;
            clone.dIntensity = dIntensity;
            clone.dProtein = dProtein;
            clone.dPeptideScore = dPeptideScore;
            clone.dFragmentScore = dFragmentScore;

            clone.dicOfMinPrecursorThr = null;
            clone.dicOfMaxPrecursorThr = null;
            clone.dicOfMinFragmentThr = null;
            clone.dicOfMaxFragmentThr = null;

            clone.initiatorMethionineBehavior = initiatorMethionineBehavior;
            clone.fixedModifications = new List<Modification>(fixedModifications);
            clone.variableModifications = new List<Modification>(variableModifications);
            clone.maximumVariableModificationIsoforms = maximumVariableModificationIsoforms;
            clone.precursorMassTolerance = precursorMassTolerance;
            clone.productMassTolerance = productMassTolerance;

            if(FastaFiles != null)
                clone.FastaFiles = new List<string>(FastaFiles).ToArray();
            if (FastaFilesMinors != null)
                clone.FastaFilesMinors = new List<string>(FastaFilesMinors).ToArray();
            if (MS2s != null)
                clone.MS2s = new List<string>(MS2s).ToArray();

            return clone;
        }

        /*
        private Dictionary<double, double> dicOfScore = null;
        public void UptimizeParams(bool currentTrendIsGood)
        {
            if (dicOfScore == null)
            {
                dicOfScore = new Dictionary<double, double>();
                dicOfScore.Add(1.0, dProduct);
                dicOfScore.Add(2.0, dPrecursor);
                dicOfScore.Add(3.0, dMatchingProductFraction);
                dicOfScore.Add(4.0, dMatchingProduct);
                dicOfScore.Add(5.0, dIntensityFraction);
                dicOfScore.Add(6.0, dIntensity);
                dicOfScore.Add(7.0, dProtein);
                dicOfScore.Add(8.0, dPeptideScore);
                dicOfScore.Add(9.0, dFragmentScore);
            }

            Utilities.Methods.GradientDescent.SolveMaxFlowStyle(unitSpectrum, mixedSpectrum, out solution, out tmpUnderflow, ConSole, stepSize);
                
            Random r = new Random();

            dProduct = r.NextDouble();
            dPrecursor = r.NextDouble();
            dMatchingProductFraction = r.NextDouble();
            dMatchingProduct = 0;// r.NextDouble();
            dIntensityFraction = r.NextDouble();
            dIntensity = 0;// r.NextDouble() * 0.0001;
            dProtein = r.NextDouble();
            dPeptideScore = r.NextDouble();
            dFragmentScore = r.NextDouble();
        }//*/

        public void RandomizeParams()
        {
            Random r = new Random();

            dProduct = r.NextDouble() - 0.5;
            dPrecursor = r.NextDouble() - 0.5;
            dMatchingProductFraction = r.NextDouble() - 0.5;
            dMatchingProduct = (r.NextDouble() - 0.5) * 0.01;
            dIntensityFraction = r.NextDouble() - 0.5;
            dIntensity = (r.NextDouble() - 0.5) * 0.0001;
            dProtein = r.NextDouble() - 0.5;
            dPeptideScore = r.NextDouble() - 0.5;
            dFragmentScore = r.NextDouble() - 0.5;
        }

        public static DBOptions Load(string fileName, IConSol consol)
        {
            DBOptions options = new DBOptions(consol);
            vsCSV csv = new vsCSV(fileName);
            foreach (string line in csv.LINES_LIST)
            {
                string[] splits = line.Split(',');
                switch (splits[0])
                {
                    case "MinimumPrecursorChargeState":
                        options.MinimumPrecursorChargeState = int.Parse(splits[1]); break;
                    case "MaximumPrecursorChargeState":
                        options.MaximumPrecursorChargeState = int.Parse(splits[1]); break;
                    case "MaximumNumberOfFragmentsPerSpectrum":
                        options.MaximumNumberOfFragmentsPerSpectrum = int.Parse(splits[1]); break;
                    case "MaximumPeptideMass": options.MaximumPeptideMass = double.Parse(splits[1]); break;
                    case "DecoyFusion": options.DecoyFusion = bool.Parse(splits[1]); break;
                    //case "DigestionEnzyme":                    options.DigestionEnzyme = double.Parse(splits[1]); break; 
                    //case "ToleratedMissedCleavages":                    options.ToleratedMissedCleavages = double.Parse(splits[1]); break; 
                    case "PSMFalseDiscoveryRate": options.PSMFalseDiscoveryRate = double.Parse(splits[1]); break;

                    case "NbPSMToKeep": options.NbPSMToKeep = int.Parse(splits[1]); break;
                    case "OutputFolder": options.OutputFolder = splits[1]; break;
                    case "MinimumPSMScore": options.MinimumPSMScore = double.Parse(splits[1]); break;
                    case "MinimumPeptideLength": options.MinimumPeptideLength = int.Parse(splits[1]); break;
                    case "MaximumPeptideLength": options.MaximumPeptideLength = int.Parse(splits[1]); break;
                    case "fullFragment":
                        List<FragmentGen> frags = new List<FragmentGen>();
                        foreach (string frag in splits[1].Split(';'))
                            if(!string.IsNullOrEmpty(frag))
                                frags.Add(FullFragments.AllFragments[frag]);
                        options.fullFragment = new FullFragments(frags);
                        break;
                    //public Fragments fragments = double.Parse(splits[1]); break; 
                    case "NbMinProducts": options.NbMinProducts = int.Parse(splits[1]); break;

                    case "addFragmentLoss": options.addFragmentLoss = bool.Parse(splits[1]); break;
                    case "addFragmentMods": options.addFragmentMods = bool.Parse(splits[1]); break;
                    case "WriteMaxQuantPeakFile": options.WriteMaxQuantPeakFile = bool.Parse(splits[1]); break;

                    case "SaveMS1Peaks": options.SaveMS1Peaks = bool.Parse(splits[1]); break;
                    case "SaveMSMSPeaks": options.SaveMSMSPeaks = bool.Parse(splits[1]); break;

                    case "LoadSpectraIfFound": options.LoadSpectraIfFound = bool.Parse(splits[1]); break;

                    case "ComputedRetentionTimeDiff": options.ComputedRetentionTimeDiff = double.Parse(splits[1]); break; //TODO compute this after alignment step, based on common identifications
                    case "EffectiveIsolationWindowRatio": options.EffectiveIsolationWindowRatio = double.Parse(splits[1]); break; //TODO optimize this value 
                    case "MinimumPrecursorIntensityRatioInIsolationWindow": options.MinimumPrecursorIntensityRatioInIsolationWindow = double.Parse(splits[1]); break;


                    case "dProduct": options.dProduct = double.Parse(splits[1]); break;
                    case "dPrecursor": options.dPrecursor = double.Parse(splits[1]); break;
                    case "dMatchingProductFraction": options.dMatchingProductFraction = double.Parse(splits[1]); break;
                    case "dMatchingProduct": options.dMatchingProduct = double.Parse(splits[1]); break;
                    case "dIntensityFraction": options.dIntensityFraction = double.Parse(splits[1]); break;
                    case "dIntensity": options.dIntensity = double.Parse(splits[1]); break;
                    case "dProtein": options.dProtein = double.Parse(splits[1]); break;
                    case "dPeptideScore": options.dPeptideScore = double.Parse(splits[1]); break;
                    case "dFragmentScore": options.dFragmentScore = double.Parse(splits[1]); break;

                    ///Values kept from the original Morpheus source code
                    //case "initiatorMethionineBehavior ":                    options.initiatorMethionineBehavior = double.Parse(splits[1]); break; 
                    case "fixedModifications":
                        options.fixedModifications = new List<Modification>();
                        foreach (string mod in splits[1].Split(';'))
                            if (!string.IsNullOrEmpty(mod))
                                options.fixedModifications.Add(ModificationDictionary.Instance[mod]);
                        break;
                    case "variableModifications":
                        options.variableModifications = new List<Modification>();
                        foreach (string mod in splits[1].Split(';'))
                            if (!string.IsNullOrEmpty(mod))
                                options.fixedModifications.Add(ModificationDictionary.Instance[mod]);
                        break;
                    case "maximumVariableModificationIsoforms": options.maximumVariableModificationIsoforms = int.Parse(splits[1]); break;
                    case "precursorMassTolerance": options.precursorMassTolerance = new MassTolerance(double.Parse(splits[1]), MassToleranceUnits.ppm); break;
                    case "productMassTolerance": options.productMassTolerance = new MassTolerance(double.Parse(splits[1]), MassToleranceUnits.Da); break;
                    case "dicOfMinPrecursorThr":
                        options.dicOfMinPrecursorThr = StrToDic(splits[1]); break;
                    case "dicOfMaxPrecursorThr":
                        options.dicOfMaxPrecursorThr = StrToDic(splits[1]); break;
                    case "dicOfMinFragmentThr":
                        options.dicOfMinFragmentThr = StrToDic(splits[1]); break;
                    case "dicOfMaxFragmentThr":
                        options.dicOfMaxFragmentThr = StrToDic(splits[1]); break;
                    case "FastaFiles":
                        options.FastaFiles = TxtToArray(splits[1]);       break;
                    case "FastaFilesMinors":
                        options.FastaFilesMinors = TxtToArray(splits[1]); break;
                    case "MS2s":
                        options.MS2s = TxtToArray(splits[1]); break;
                }
            }
            return options;
        }

        public void Save(string fileName)
        {
            vsCSVWriter writer = new vsCSVWriter(fileName);
            writer.AddLine( "MinimumPrecursorChargeState," +  MinimumPrecursorChargeState);
            writer.AddLine( "MaximumPrecursorChargeState,"+ MaximumPrecursorChargeState);
            writer.AddLine( "MaximumNumberOfFragmentsPerSpectrum,"+ MaximumNumberOfFragmentsPerSpectrum);
            writer.AddLine( "MaximumPeptideMass,"+ MaximumPeptideMass);
            writer.AddLine( "DecoyFusion,"+ DecoyFusion);

            writer.AddLine( "DigestionEnzyme,"+ DigestionEnzyme);
            writer.AddLine( "ToleratedMissedCleavages,"+ ToleratedMissedCleavages);

            writer.AddLine( "PSMFalseDiscoveryRate,"+ PSMFalseDiscoveryRate);

            writer.AddLine( "NbPSMToKeep,"+ NbPSMToKeep);
            writer.AddLine( "OutputFolder,"+ OutputFolder);
            writer.AddLine( "MinimumPSMScore,"+ MinimumPSMScore);
            writer.AddLine( "MinimumPeptideLength,"+ MinimumPeptideLength);
            writer.AddLine( "MaximumPeptideLength,"+ MaximumPeptideLength);

            writer.AddLine( "fullFragment,"+ fullFragment.GetFragmentSTR());
            //public Fragments fragments;
            writer.AddLine( "NbMinProducts,"+ NbMinProducts);

            writer.AddLine( "addFragmentLoss,"+ addFragmentLoss);
            writer.AddLine( "addFragmentMods,"+ addFragmentMods);
            writer.AddLine( "WriteMaxQuantPeakFile,"+ WriteMaxQuantPeakFile);

            writer.AddLine( "SaveMS1Peaks,"+ SaveMS1Peaks);
            writer.AddLine( "SaveMSMSPeaks,"+ SaveMSMSPeaks);

            writer.AddLine( "LoadSpectraIfFound,"+ LoadSpectraIfFound);

            writer.AddLine( "ComputedRetentionTimeDiff,"+ ComputedRetentionTimeDiff);//TODO compute this after alignment step, based on common identifications
            writer.AddLine( "EffectiveIsolationWindowRatio,"+ EffectiveIsolationWindowRatio);//TODO optimize this value 
            writer.AddLine( "MinimumPrecursorIntensityRatioInIsolationWindow,"+ MinimumPrecursorIntensityRatioInIsolationWindow);


            writer.AddLine("dProduct," + dProduct);
            writer.AddLine("dPrecursor," + dPrecursor);
            writer.AddLine("dMatchingProductFraction," + dMatchingProductFraction);
            writer.AddLine("dMatchingProduct," + dMatchingProduct);
            writer.AddLine("dIntensityFraction," + dIntensityFraction);
            writer.AddLine("dIntensity," + dIntensity);
            writer.AddLine("dProtein," + dProtein);
            writer.AddLine("dPeptideScore," + dPeptideScore);
            writer.AddLine("dFragmentScore," + dFragmentScore);

            ///Values kept from the original Morpheus source code
            writer.AddLine(  "initiatorMethionineBehavior ," +initiatorMethionineBehavior);
            writer.AddLine("fixedModifications," + ModsToStr(fixedModifications));
            writer.AddLine(  "variableModifications ," + ModsToStr(variableModifications));
            writer.AddLine(  "maximumVariableModificationIsoforms," + maximumVariableModificationIsoforms);
            writer.AddLine(  "precursorMassTolerance," + precursorMassTolerance.Value);
            writer.AddLine(  "productMassTolerance," +productMassTolerance.Value);
            writer.AddLine(  "dicOfMinPrecursorThr,"+ThrDicToStr(dicOfMinPrecursorThr));
            writer.AddLine(  "dicOfMaxPrecursorThr,"+ThrDicToStr(dicOfMaxPrecursorThr));
            writer.AddLine(  "dicOfMinFragmentThr,"+ThrDicToStr(dicOfMinFragmentThr));
            writer.AddLine(  "dicOfMaxFragmentThr,"+ThrDicToStr(dicOfMaxFragmentThr));

            writer.AddLine("FastaFiles," + ThrArrayToStr(FastaFiles));
            writer.AddLine("FastaFilesMinors," + ThrArrayToStr(FastaFilesMinors));
            writer.AddLine("MS2s," + ThrArrayToStr(MS2s));
            writer.WriteToFile();
        }
        private static Dictionary<string, double> StrToDic(string linePart)
        {
            Dictionary<string, double> dic = new Dictionary<string,double>();
            string[] items = linePart.Split(';');
            foreach(string item in items)
            {
                if (!string.IsNullOrEmpty(item))
                {
                    string[] splits = item.Split('|');
                    dic.Add(splits[0], double.Parse(splits[1]));
                }
            }
            return dic;
        }
        private static string ModsToStr(IEnumerable<Modification> mods)
        {
            string line = "";
            foreach (Modification mod in mods)
                line += mod.Description + ";";
            return line;
        }
        private static string ThrDicToStr(Dictionary<string, double> dic)
        {
            string line = "";
            foreach (string key in dic.Keys)
                line += key + "|" + dic[key] + ";";
            return line;
        }
        private static string ThrArrayToStr(string[] array)
        {
            string line = "";
            if(array != null)
                foreach (string key in array)
                    line += key + ";";
            return line;
        }
        private static string[] TxtToArray(string line)
        {
            List<string> list1 = new List<string>();
            foreach (string file in line.Split(';'))
                if (!string.IsNullOrEmpty(file))
                    list1.Add(file);
            return list1.ToArray();
        }
    }
}
