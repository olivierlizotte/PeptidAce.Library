using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace PeptidAce
{
    public class AnnotatedSpectrumModification
    {
        public int index;
        public double modMass;
        public string aminoAcid;
        public AnnotatedSpectrumModification(char AA, double massShift, int position)
        {
            index = position;
            modMass = massShift;
            aminoAcid = AA.ToString();
        }
    }
    public class AnnotatedSpectrum
    {
        public string sequence;
        public double[][] peaks;
        public int scanNum;
        public int charge;
        public double precursorMz;
        public string fileName;
        public string Type = "AnnotatedSpectrum";

        public double ntermMod = 0;
        public double ctermMod = 0;
        public AnnotatedSpectrumModification[] staticMods;
        public AnnotatedSpectrumModification[] variableMods;

        public AnnotatedSpectrum(Sample sample, ProductSpectrum spectrum, Peptide peptide)
        {
            sequence = peptide.BaseSequence;
            peaks = new double[spectrum.Peaks.Count][];
            for (int i = 0; i < peaks.Length; i++)
            {
                peaks[i] = new double[2];
                peaks[i][0] = spectrum.Peaks[i].MZ;
                peaks[i][1] = spectrum.Peaks[i].Intensity;
            }
            scanNum = spectrum.ScanNumber;
            charge = spectrum.PrecursorCharge;
            precursorMz = spectrum.PrecursorMZ;
            fileName = sample.sSDF;

            //Fixed modifications
            List<AnnotatedSpectrumModification> fixedMods = new List<AnnotatedSpectrumModification>();
            if(peptide.FixedModifications != null)
                foreach(int position in peptide.FixedModifications.Keys)
                {
                    if(position == 1)
                        foreach(Modification mod in peptide.FixedModifications[position])
                            ntermMod += mod.MonoisotopicMassShift;
                    else
                        if(position == sequence.Length + 2)
                            foreach(Modification mod in peptide.FixedModifications[position])
                                ctermMod += mod.MonoisotopicMassShift;
                        else
                            foreach(Modification mod in peptide.FixedModifications[position])
                                fixedMods.Add(new AnnotatedSpectrumModification(mod.AminoAcid, mod.MonoisotopicMassShift, position-1));
                }
            staticMods = fixedMods.ToArray();

            //Variable modifications
            List<AnnotatedSpectrumModification> varMods = new List<AnnotatedSpectrumModification>();
            if(peptide.VariableModifications != null)
                foreach(int position in peptide.VariableModifications.Keys)
                {
                    if(position == 1)
                        ntermMod += peptide.VariableModifications[position].MonoisotopicMassShift;
                    else
                        if(position == sequence.Length + 2)
                            ctermMod += peptide.VariableModifications[position].MonoisotopicMassShift;
                        else
                            varMods.Add(new AnnotatedSpectrumModification(peptide.VariableModifications[position].AminoAcid, peptide.VariableModifications[position].MonoisotopicMassShift, position-1));
                }
            variableMods = varMods.ToArray();

            /*
    var sequence = "FDSFGDLSSASAIMGNPK";
    var varMods = [];
    // modification index = 14; modification mass = 16.0; modified residue = 'M'
    varMods[0] = { index: 14, modMass: 16.0, aminoAcid: 'M' };
    // mass to be added to the N-terminus
    var ntermMod = 164.07;//*/
        }
    }
}
