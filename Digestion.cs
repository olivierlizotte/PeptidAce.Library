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
using System.Threading.Tasks;

namespace PeptidAce
{
    /// <summary>
    /// Methods to digest proteins into peptides, either with an enzyme, or with the noEnzyme, GPU extension
    /// </summary>
    public class Digestion
    {
        private DBOptions options;
        public Digestion(DBOptions dbOptions)
        {
            this.options = dbOptions;
        }

        public static double[] GetMasses(string sequence)
        {
            double[] proteinMasses = new double[sequence.Length];
            for (int i = 0; i < sequence.Length; i++)
                proteinMasses[i] = AminoAcidMasses.GetMonoisotopicMass(sequence[i]);
            return proteinMasses;
        }

        private bool IsDone = true;
        private System.Collections.Concurrent.ConcurrentQueue<Tuple<Peptide, int>> CumulMatch;
        /// <summary>
        /// Runs the protein digestion on the GPU. Does not support modifications for now.
        /// </summary>
        /// <param name="options"></param>
        /// <param name="Proteins"></param>
        /// <param name="listOfQueries"></param>
        private void RunGPUProteinDigest(DBOptions options, List<Protein> Proteins, Queries listOfQueries)
        {
            List<double> precursorMasses = new List<double>(listOfQueries.Count);
            //double[] precursors = new double[listOfQueries.Count];
            for (int i = 0; i < listOfQueries.Count; i++)
                precursorMasses.Add(listOfQueries[i].precursor.Mass);//TODO Why waterLoss here? - Constants.WATER_MONOISOTOPIC_MASS);
            //precursors[i] = listOfQueries[i].precursor.Mass - Constants.WATER_MONOISOTOPIC_MASS;
            precursorMasses.Sort();

            //for each protein, build matrix of mass
            int nbProteins = 0;
            GPU.ProteinDigest pg = new GPU.ProteinDigest(precursorMasses.ToArray(), options.MaximumPeptideLength, options.MinimumPeptideLength);
            foreach (Protein protein in Proteins)
            {
                double[] proteinMasses = GetMasses(protein.BaseSequence);

                foreach (GPU.ProteinPrecursorMatch match in pg.Execute(proteinMasses, options.precursorMassTolerance.Value, options.MaximumPeptideMass))//TODO compute correct tolerance window
                {
                    CumulMatch.Enqueue(new Tuple<Peptide, int>(new Peptide(protein, match.proteinStartPos, match.proteinEndPos, 0), match.firstQueryIndex));
                    //yield return new Tuple<Peptide, int> new Peptide(protein, match.proteinStartPos, match.proteinEndPos, 0);//TODO add modifications
                }

                nbProteins++;
            }
            pg.Dispose();
        }

        /// <summary>
        /// Runs the NoEnzyme protein digestion on the GPU
        /// </summary>
        /// <param name="Proteins"></param>
        /// <param name="listOfQueries"></param>
        /// <returns></returns>
        public IEnumerable<Tuple<Peptide, int>> DigestProteomeOnTheFlyNoEnzyme(List<Protein> Proteins, Queries listOfQueries)
        {
            IsDone = false;
            CumulMatch = new System.Collections.Concurrent.ConcurrentQueue<Tuple<Peptide, int>>();

            Task.Factory.StartNew(() =>
            {
                RunGPUProteinDigest(options, Proteins, listOfQueries);
                IsDone = true;
            });

            Tuple<Peptide, int> item;
            while (!IsDone)
            {
                if (CumulMatch.TryDequeue(out item))
                    yield return item;
            }

            while (CumulMatch.TryDequeue(out item))
            {
                yield return item;
            }
        }

		public List<List<double>> GetProteinMasses(Protein protein, DBOptions options)
		{
			List<List<double>> masses = new List<List<double>> (protein.BaseSequence.Length);
			foreach(char aa in protein.BaseSequence)
			{
				double mass = AminoAcidMasses.GetMonoisotopicMass(aa);
				List<double> mods = new List<double> ();
				foreach (Modification mod in options.fixedModifications)
					if (mod.AminoAcid == aa)
						mods.Add (mass + mod.MonoisotopicMassShift);
				if(mods.Count == 0)
					mods.Add(mass);
				foreach (Modification mod in options.variableModifications)
					if (mod.AminoAcid == aa)
						mods.Add (mass + mod.MonoisotopicMassShift);
				masses.Add (mods);
			}

			//Protein NTerm Mods
			foreach (Modification mod in options.fixedModifications)
				if (mod.Type == ModificationType.ProteinNTerminus)
					masses [0][0] = masses [0] [0] + mod.MonoisotopicMassShift;
			foreach (Modification mod in options.variableModifications)
				if (mod.Type == ModificationType.ProteinNTerminus)
					masses [0].Add (masses [0] [0] + mod.MonoisotopicMassShift);

			//Protein CTerm Mods
			foreach (Modification mod in options.fixedModifications)
				if (mod.Type == ModificationType.ProteinCTerminus)
					masses [masses.Count - 1][0] = masses [masses.Count - 1] [0] + mod.MonoisotopicMassShift;
			foreach (Modification mod in options.variableModifications)
				if (mod.Type == ModificationType.ProteinCTerminus)
					masses [masses.Count - 1].Add (masses [masses.Count - 1] [0] + mod.MonoisotopicMassShift);
			return masses;
		}

		private IEnumerable<double[]> recur_(List<List<double>> protMasses, int startPos, int length, double[] soFar, int pos)
		{
			if(pos == length)
				yield return soFar;
			else
			{
				while(pos < length && protMasses [pos + startPos].Count == 1)
				{
					soFar [pos] = protMasses [pos + startPos][0];
					pos++;					
				}
				if (pos >= length)
					yield return soFar;
				else
				{
					foreach (double mass in protMasses[pos + startPos])
					{
						soFar[pos] = mass;
						recur_ (protMasses, startPos, length, soFar, pos + 1);
					}
				}
			}
		}

		private IEnumerable<double[]> recur_debug(List<List<double>> protMasses, int startPos, int length, List<double> soFar, int pos)
		{
			if (pos == length)
				yield return (soFar.ToArray());
			else
			{
				while(pos < length && protMasses [pos + startPos].Count == 1)
				{
					soFar [pos] = protMasses [pos + startPos][0];
					pos++;					
				}
				if (pos >= length)
					yield return (soFar.ToArray());
				else
				{
					foreach (double mass in protMasses[pos + startPos])
					{
						soFar[pos] = mass;
						recur_debug (protMasses, startPos, length, soFar, pos + 1);
					}
				}
			}
		}
		private IEnumerable<int[]> recur_index(List<List<double>> protMasses, int startPos, int length, int[] soFar, int pos)
		{
			if(pos == length)
				yield return soFar;
			else
			{
				while(pos < length && protMasses [pos + startPos].Count == 1)
				{
					soFar [pos] = 0;//protMasses [pos + startPos][0];
					pos++;					
				}
				if (pos == length)
					yield return soFar;
				else
				{
					for(int i = 0; i < protMasses[pos + startPos].Count; i++)
					{
						soFar[pos] = i;
						recur_index (protMasses, startPos, length, soFar, pos + 1);
					}
				}
			}
		}

		public Peptide GetPeptide(PeptideHit hit)
		{
			Peptide peptide = new Peptide (hit.protein, hit.startIndex, hit.startIndex + hit.masses.Length - 1, hit.missedCleavage);
			if (peptide.MonoisotopicMass != hit.mass) {
				//Find modifications
				Dictionary<int, List<Modification>> fixedMods = new  Dictionary<int, List<Modification>> ();
				Dictionary<int, Modification> varMods = new  Dictionary<int, Modification> ();

				double[] masses = peptide.GetMasses ();
				for (int i = 0; i < masses.Length; i++) {
					if (masses [i] != hit.masses [i]) {
						foreach (Modification mod in options.fixedModifications)
							if (mod.MonoisotopicMassShift + masses [i] == hit.masses[i]) {
								masses [i] += mod.MonoisotopicMassShift;
								List<Modification> list = new List<Modification> ();
								list.Add (mod);
								fixedMods.Add (i, list);
							}

						foreach (Modification mod in options.variableModifications)
							if (mod.MonoisotopicMassShift + masses [i] == hit.masses[i]) {
								masses [i] += mod.MonoisotopicMassShift;
								varMods.Add (i, mod);
							}
					}
				}
				peptide.SetFixedModifications (fixedMods);
				peptide.SetVariableModifications (varMods);

				if (peptide.MonoisotopicMass != hit.mass)
					peptide = peptide;
			}
			return peptide;
		}

		public IEnumerable<PeptideHit> DigestFast(List<Protein> proteins, Queries AllQueries, DBOptions options)
        {
			foreach(Protein protein in proteins)
			{
				List<List<double>> masses = GetProteinMasses (protein, options);

				List<int> indices = options.DigestionEnzyme.GetDigestionSiteIndices(protein);
				indices.Insert(0, 0);
				indices.Add(protein.Length - 1);

				for(int i = 0; i < indices.Count; i++)
				{
					int missedCleavage = 0;
					for (int j = i + 1; j < indices.Count; j++) { 
						int pepLength = indices [j] - indices [i];
						if (pepLength > options.MinimumPeptideLength) {
							if (pepLength > options.MaximumPeptideLength)
								break;

							foreach (double[] pepMasses in recur_debug(masses, indices[i], pepLength, new List<double>(new double[pepLength]), 0)) {
								//foreach (int[] pepIndexes in recur_(masses, indices[i], pepLength, new double[pepLength], 0)) {
								double sum = Constants.WATER_MONOISOTOPIC_MASS;
								foreach (double item in pepMasses)
									sum += item;
								int firstIndex = AllQueries.BinarySearch (MassTolerance.MzFloor (sum, options.precursorMassTolerance));
								if (firstIndex >= 0 && firstIndex < AllQueries.Count)
									yield return new PeptideHit (protein, indices[i], pepMasses, sum, firstIndex, missedCleavage);
							}
						}
						missedCleavage++;
					}
                }            
            }
        }

		public IEnumerable<Tuple<Peptide, int>> DigestProteomeOnTheFlyBKP(List<Protein> proteins, bool allowSNP, Queries AllQueries)
		{
			foreach (Peptide peptide in ProteinSearcher.ProteinDigest(options, proteins, allowSNP))
			{
				foreach (Peptide modPeptide in peptide.GetVariablyModifiedPeptides(options.variableModifications, options.maximumVariableModificationIsoforms))
				{
					modPeptide.SetFixedModifications(options.fixedModifications);
					int firstIndex = AllQueries.BinarySearch(MassTolerance.MzFloor(modPeptide.MonoisotopicMass, options.precursorMassTolerance));
					if (firstIndex >= 0 && firstIndex < AllQueries.Count)
						yield return new Tuple<Peptide, int>(modPeptide, firstIndex);
				}            
			}
		}

        public IEnumerable<Tuple<Peptide, int>> DigestProteomeOnTheFlyFast(List<Protein> proteins, bool allowSNP, Queries AllQueries)
        {
            foreach (Peptide peptide in ProteinSearcher.ProteinDigest(options, proteins, allowSNP))
            {
                //int firstIndex = AllQueries.BinarySearch(MassTolerance.MzFloor(peptide.MonoisotopicMass, options.precursorMassTolerance));
                //if (firstIndex >= 0 && firstIndex < AllQueries.Count)
                //    yield return new Tuple<Peptide, int>(peptide, firstIndex);

                //foreach (Peptide peptide in ProteinSearcher.ProteinDigestNoEnzyme(dbOptions, proteins, AllQueries))
                //if (!TargetPeptides.Contains(peptide.BaseSequence))
                //{
                foreach (Peptide modPeptide in peptide.GetVariablyModifiedPeptides(options.variableModifications, options.maximumVariableModificationIsoforms))
                {
                    modPeptide.SetFixedModifications(options.fixedModifications);
                    int firstIndex = AllQueries.BinarySearch(MassTolerance.MzFloor(modPeptide.MonoisotopicMass, options.precursorMassTolerance));
                    if (firstIndex >= 0 && firstIndex < AllQueries.Count)
                        yield return new Tuple<Peptide, int>(modPeptide, firstIndex);
                }

                //TODO check if this favors targets over decoys since proteins are sorted target->deco
                //    TargetPeptides.Add(peptide.BaseSequence);
                //}            
            }
        }
    }
}
