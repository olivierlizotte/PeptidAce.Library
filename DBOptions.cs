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
        public int MinimumPrecursorChargeState;
        public int MaximumPrecursorChargeState;
        public int MaximumNumberOfFragmentsPerSpectrum;
        public string FastaDatabaseFilepath;
        public double MaximumPeptideMass;
        public bool DecoyFusion;

        public Protease DigestionEnzyme;
        public int ToleratedMissedCleavages;

        public double PSMFalseDiscoveryRate = 0.05;

        public int NbPSMToKeep = 32;
        public string OutputFolder;
        public double MinimumPSMScore;
        public int MinimumPeptideLength = 5;
        public int MaximumPeptideLength = 100;

        public FullFragments fullFragment = new FullFragments(true);
        //public Fragments fragments;
        public int NbMinProducts = 4;

        public bool addFragmentLoss = false;
        public bool addFragmentMods = false;
        public bool NoEnzymeSearch = false;
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

        ///Values kept from the original Morpheus source code
        public InitiatorMethionineBehavior initiatorMethionineBehavior;
        public List<Modification> fixedModifications;
        public List<Modification> variableModifications;
        public int maximumVariableModificationIsoforms;
        public MassTolerance precursorMassTolerance;
        public MassTolerance productMassTolerance;

        public IConSol ConSole;
        /// <summary>
        /// Parameter less constructor used for save states
        /// </summary>
        public DBOptions(IConSol console = null)
        {
            fixedModifications = new List<Modification>();
            variableModifications = new List<Modification>();
            //fragments = new Fragments();
            if (console == null)
                ConSole = new ConSolCommandLine();
            else
                ConSole = console;
        }

        public DBOptions(string fasta, IConSol console = null)
        {
            if (console == null)
                ConSole = new ConSolCommandLine();
            else
                ConSole = console;
            //Create with default values
            this.DecoyFusion = true;
            this.FastaDatabaseFilepath = fasta;
            this.MaximumPeptideMass = 10000;
            ProteaseDictionary proteases = ProteaseDictionary.Instance;
            this.DigestionEnzyme = proteases["no enzyme"];// proteases["trypsin (no proline rule)"];
            this.NoEnzymeSearch = true;
            this.ToleratedMissedCleavages = 100;// 3;//determines the length of peptides with no-enzyme option
            this.initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            this.fixedModifications = new List<Modification>();
            this.variableModifications = new List<Modification>();
            this.maximumVariableModificationIsoforms = 1024;
            
            this.MinimumPrecursorChargeState = 1;
            this.MaximumPrecursorChargeState = 4;
            this.MaximumNumberOfFragmentsPerSpectrum = 400;
            //TODO Add precision to the precursor by reading MS part of file
            this.precursorMassTolerance = new MassTolerance(0.005, MassToleranceUnits.Da);//2.1
            //TODO Add precision to the product masses by reading corresponding MS part of raw file
            this.productMassTolerance = new MassTolerance(0.005, MassToleranceUnits.Da);

            this.PSMFalseDiscoveryRate = 0.25;// 0.05;

            this.OutputFolder = @"C:\_IRIC\DATA\Test2";//C:\Documents and Settings\ProteoAdmin\Desktop\AEffacer\Morpheus\Output";
            this.MinimumPSMScore = 0.0001;         
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

            dProduct = r.NextDouble();
            dPrecursor = r.NextDouble();
            dMatchingProductFraction = r.NextDouble();
            dMatchingProduct = 0;// r.NextDouble();
            dIntensityFraction = r.NextDouble();
            dIntensity = 0;// r.NextDouble() * 0.0001;
            dProtein = r.NextDouble();
            dPeptideScore = r.NextDouble();
            dFragmentScore = r.NextDouble();
        }

        public void Save(string fileName)
        {
            vsCSVWriter writer = new vsCSVWriter(fileName);
            writer.AddLine( "MinimumPrecursorChargeState = " +  MinimumPrecursorChargeState);
            writer.AddLine( "MaximumPrecursorChargeState"+ MaximumPrecursorChargeState);
            writer.AddLine( "MaximumNumberOfFragmentsPerSpectrum"+ MaximumNumberOfFragmentsPerSpectrum);
            writer.AddLine( "FastaDatabaseFilepath"+ FastaDatabaseFilepath);
            writer.AddLine( "MaximumPeptideMass"+ MaximumPeptideMass);
            writer.AddLine( "DecoyFusion"+ DecoyFusion);

            writer.AddLine( "DigestionEnzyme"+ DigestionEnzyme);
            writer.AddLine( "ToleratedMissedCleavages"+ ToleratedMissedCleavages);

            writer.AddLine( "PSMFalseDiscoveryRate"+ PSMFalseDiscoveryRate);

            writer.AddLine( "NbPSMToKeep"+ NbPSMToKeep);
            writer.AddLine( "OutputFolder"+ OutputFolder);
            writer.AddLine( "MinimumPSMScore"+ MinimumPSMScore);
            writer.AddLine( "MinimumPeptideLength"+ MinimumPeptideLength);
            writer.AddLine( "MaximumPeptideLength"+ MaximumPeptideLength);

            writer.AddLine( "fullFragment"+ fullFragment);
            //public Fragments fragments;
            writer.AddLine( "NbMinProducts"+ NbMinProducts);

            writer.AddLine( "addFragmentLoss"+ addFragmentLoss);
            writer.AddLine( "addFragmentMods"+ addFragmentMods);
            writer.AddLine( "NoEnzymeSearch"+ NoEnzymeSearch);
            writer.AddLine( "WriteMaxQuantPeakFile"+ WriteMaxQuantPeakFile);

            writer.AddLine( "SaveMS1Peaks"+ SaveMS1Peaks);
            writer.AddLine( "SaveMSMSPeaks"+ SaveMSMSPeaks);

            writer.AddLine( "LoadSpectraIfFound"+ LoadSpectraIfFound);

            writer.AddLine( "ComputedRetentionTimeDiff"+ ComputedRetentionTimeDiff);//TODO compute this after alignment step, based on common identifications
            writer.AddLine( "EffectiveIsolationWindowRatio"+ EffectiveIsolationWindowRatio);//TODO optimize this value 
            writer.AddLine( "MinimumPrecursorIntensityRatioInIsolationWindow"+ MinimumPrecursorIntensityRatioInIsolationWindow);


            writer.AddLine("dProduct =" + dProduct + ";");
            writer.AddLine("dPrecursor =" + dPrecursor + ";");
            writer.AddLine("dMatchingProductFraction =" + dMatchingProductFraction + ";");
            writer.AddLine("dMatchingProduct =" + dMatchingProduct + ";");
            writer.AddLine("dIntensityFraction =" + dIntensityFraction + ";");
            writer.AddLine("dIntensity =" + dIntensity + ";");
            writer.AddLine("dProtein =" + dProtein + ";");
            writer.AddLine("dPeptideScore =" + dPeptideScore + ";");
            writer.AddLine("dFragmentScore =" + dFragmentScore + ";");

            ///Values kept from the original Morpheus source code
            writer.AddLine(  "initiatorMethionineBehavior ," +initiatorMethionineBehavior);
            writer.AddLine(  "fixedModifications ," +fixedModifications);
            writer.AddLine(  "variableModifications ," +variableModifications);
            writer.AddLine(  "maximumVariableModificationIsoforms," + maximumVariableModificationIsoforms);
            writer.AddLine(  "precursorMassTolerance,," + precursorMassTolerance);
            writer.AddLine(  "productMassTolerance," +productMassTolerance);
            writer.WriteToFile();
        }
    }
}
