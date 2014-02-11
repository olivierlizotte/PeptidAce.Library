using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace PeptidAce.Utilities.Methods
{
    public class FDRizer
    {
        /// <summary>
        /// Computes the last index of an ordered list where the FDR is below the desired value
        /// </summary>
        /// <param name="list"></param>
        /// <param name="desired_fdr"></param>
        /// <returns></returns>
        public static int Extract<_T>(List<_T> list, double desired_fdr) where _T : ITargetDecoy
        {
            int cumulDecoy = 0;
            int cumulTarget = 0;
            int bestIndex = -1;
            for (int index = 0; index < list.Count; index++)
            {
                if (list[index].Target)
                    cumulTarget++;
                else
                    cumulDecoy++;

                if (list[index].Target && cumulDecoy / (double)cumulTarget <= desired_fdr)
                    bestIndex = index;
            }
            return bestIndex;
        }//*/
    }
}
