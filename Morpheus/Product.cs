/*
 * Wenger CD, Coon JJ. A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra, Journal of Proteome Research, 2013; 12(3): 1377-86
 * http://www.chem.wisc.edu/~coon/software.php#morpheus
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 * Altered by Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 */
namespace PeptidAce
{
    public class Product
    {
        public ProductType ProductType { get; private set; }

        public int Number { get; private set; }

        public double Mass { get; private set; }

        public Product(ProductType productType, int number, double mass)
        {
            ProductType = productType;
            Number = number;
            Mass = mass;
        }

        public override string ToString()
        {
            return ProductType.ToString() + Number.ToString();
        }
    }
}