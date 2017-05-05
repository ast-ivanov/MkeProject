using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MkeProject
{
    static class Sort
    {
        public static void SortPairs(this List<Pair> pairs, int N)
        {
            for (int i = 0; i < pairs.Count - 1; ++i)
            {
                for (int j = i + 1; j < pairs.Count; ++j)
                {
                    if (pairs[j].First * N + pairs[j].Second < pairs[i].First * N + pairs[i].Second)
                    {
                        var temp = pairs[i];
                        pairs[i] = pairs[j];
                        pairs[j] = temp;
                    }
                }
            }
        }
    }
}
