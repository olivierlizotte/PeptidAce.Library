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

        public FullFragments fullFragment = new FullFragments();
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
    }
}
