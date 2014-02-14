using System;
using System.Collections.Generic;

namespace PeptidAce
{
	public class PSM
	{
		public PeptideHit peptide;
		public Query query;
		public PSM (PeptideHit hit, Query pQuery)
		{
			peptide = hit;
			query = pQuery;
			Initialize ();
		}

		public int MatchingProducts;
		List<ProductMatch> AllProductMatch = new List<ProductMatch>();
		public void Initialize()
		{
			MatchingProducts = 0;
			foreach (ProductMatch match in GetProductMZs(query.options, query.spectrum.Peaks))
			{
				MatchingProducts++;
				AllProductMatch.Add (match);
			}
		}


		public IEnumerable<ProductMatch> GetProductMZs(DBOptions options, List<MsMsPeak> peaks)//, List<double> theoretical_product_mzs = null)
		{
			int num_experimental_peaks = peaks.Count;

			foreach (ProductMatch matchTheo in options.fullFragment.ComputeFragments(peptide.masses, query.precursor.Charge))
			{
				double massDiff = options.productMassTolerance.Value;
				double bestMz = -1;
				double bestInt = 0;

				foreach (int index in query.spectrum.GetIndexOfMZInRange(matchTheo.theoMz, options.productMassTolerance))
				{
					if (peaks[index].Charge <= 0 || peaks[index].Charge == matchTheo.charge)
					{

						double diff = PeptidAce.Utilities.Numerics.CalculateMassError(peaks[index].MZ, matchTheo.theoMz, options.productMassTolerance.Units);
						if (Math.Abs(diff) < options.productMassTolerance.Value)
						{
							if (Math.Abs(diff) < Math.Abs(massDiff))//TODO Priority to intensity, or precision?
							{
								massDiff = diff;
								bestMz = peaks[index].MZ;
							}
							bestInt += peaks[index].Intensity;
						}
					}
				}
				if (bestMz >= 0)
				{
					ProductMatch pMatch = new ProductMatch();
					pMatch.weight = matchTheo.weight;
					pMatch.theoMz = matchTheo.theoMz;// Utilities.MZFromMzSingleCharge(theoMass, charge);
					pMatch.obsMz = bestMz;// experimental_masses[bestIndex];
					pMatch.mass_diff = massDiff;
					pMatch.obsIntensity = bestInt;// Intensities[bestIndex];
					pMatch.charge = matchTheo.charge;// peaks[bestIndex].Charge;
					pMatch.Fragment = matchTheo.Fragment;
					pMatch.fragmentPos = matchTheo.fragmentPos;
					pMatch.normalizedIntensity = pMatch.obsIntensity / (query.spectrum.InjectionTime * query.spectrum.PrecursorIntensityPerMilliSecond);
					yield return pMatch;
				}
			}
		}
	}
}

