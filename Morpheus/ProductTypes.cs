/*
 * Wenger CD, Coon JJ. A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra, Journal of Proteome Research, 2013; 12(3): 1377-86
 * http://www.chem.wisc.edu/~coon/software.php#morpheus
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 * Altered by Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 */
using System;
using System.Collections.Generic;
using System.IO;

namespace PeptidAce
{
    public class ProductTypes
    {
        private static readonly Dictionary<string, ProductType[]> PRODUCT_TYPES = new Dictionary<string, ProductType[]>();

        private static readonly ProductTypes instance = new ProductTypes();

        private ProductTypes()
        {
            using (StreamReader product_types = new StreamReader(Path.Combine(Path.GetDirectoryName(Environment.GetCommandLineArgs()[0]), "Configuration", "product_types.tsv")))
            {
                string header = product_types.ReadLine();
                string[] header_fields = header.Split('\t');
                Dictionary<int, ProductType> product_type_column_numbers = new Dictionary<int, ProductType>(Enum.GetNames(typeof(ProductType)).Length);
                for(int i = 1; i < header_fields.Length; i++)
                {
                    ProductType product_type = (ProductType)Enum.Parse(typeof(ProductType), header_fields[i]);
                    product_type_column_numbers.Add(i, product_type);
                }

                while(product_types.Peek() != -1)
                {
                    string line = product_types.ReadLine();
                    string[] fields = line.Split('\t');

                    string fragmentation_method = fields[0];
                    List<ProductType> frag_product_types = new List<ProductType>();
                    for(int i = 1; i < fields.Length; i++)
                    {
                        if(double.Parse(fields[i]) > 0)
                        {
                            frag_product_types.Add(product_type_column_numbers[i]);
                        }
                    }
                    PRODUCT_TYPES.Add(fragmentation_method, frag_product_types.ToArray());
                }
            }
        }

        public static ProductTypes Instance
        {
            get { return instance; }
        }

        public ProductType[] this[string fragmentationMethod]
        {
            get { return PRODUCT_TYPES[fragmentationMethod]; }
        }
    }
}