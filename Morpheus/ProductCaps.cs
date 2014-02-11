/*
 * Wenger CD, Coon JJ. A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra, Journal of Proteome Research, 2013; 12(3): 1377-86
 * http://www.chem.wisc.edu/~coon/software.php#morpheus
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 * Altered by Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 */
using System;
using System.IO;

namespace PeptidAce
{
    public class ProductCaps
    {
        private static readonly double[,] PRODUCT_CAP_MASSES = new double[Enum.GetNames(typeof(ProductType)).Length, Enum.GetNames(typeof(MassType)).Length];

        private static readonly ProductCaps instance = new ProductCaps();

        private ProductCaps()
        {
            using (StreamReader product_caps = new StreamReader(Path.Combine(Path.GetDirectoryName(Environment.GetCommandLineArgs()[0]), "Configuration", "product_caps.tsv")))
            {
                string header = product_caps.ReadLine();

                while(product_caps.Peek() != -1)
                {
                    string line = product_caps.ReadLine();
                    string[] fields = line.Split('\t');

                    ProductType product_type = (ProductType)Enum.Parse(typeof(ProductType), fields[0], true);
                    double cap_monoisotopic_mass = double.Parse(fields[1]);
                    PRODUCT_CAP_MASSES[(int)product_type, (int)MassType.Monoisotopic] = cap_monoisotopic_mass;
                    double cap_average_mass = double.Parse(fields[2]);
                    PRODUCT_CAP_MASSES[(int)product_type, (int)MassType.Average] = cap_average_mass;
                }
            }
        }

        public static ProductCaps Instance
        {
            get { return instance; }
        }

        public double this[ProductType productType, MassType massType]
        {
            get { return PRODUCT_CAP_MASSES[(int)productType, (int)massType]; }
        }
    }
}