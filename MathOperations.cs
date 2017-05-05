using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MkeProject
{
    static class MathOperations
    {

        public static void Swap(ref int a, ref int b)
        {
            int temp = a;
            a = b;
            b = temp;
        }
        public static double ScalarMult(double[] a, double[] b)
        {
            double result = 0;
            if (a.Length != b.Length)
            {
                throw new ArgumentException("Arrays sizes not equal");
            }
            else
            {
                for (int i = 0; i < a.Length; i++)
                {
                    result += a[i] * b[i];
                }
            }
            return result;
        }
        public static double[] Matrix_Mult(Matrix matrix, int N, double[] vector)
        {
            double[] result = new double[vector.Length];
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    result[i] += matrix(i, j) * vector[j];
                }
            }
            return result;
        }
    }
}
