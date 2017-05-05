using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MkeProject
{
    delegate ref double Matrix(int i, int j);
    class SLAE
    {
        private int N;
        private double[] ggl, ggu, di;
        private int[] ig, jg;
        private double zero = 0;

        public double[] b { get; private set; }
        public double[] q { get; private set; }
        public SLAE(int n, double[] ggl, int[] ig, int[] jg)
        {
            this.ggl = ggl;
            this.ggu = new double[ggl.Length];
            this.di = new double[n];
            this.b = new double[n];
            this.q = new double[n];
            this.ig = ig;
            this.jg = jg;
            this.N = n;
        }
        public void Calculate_LOS()
        {
            double alpha = 0;
            double beta = 0;
            int count = 0;
            const int maxiter = 200;
            double nev = 0;
            const double eps = 1e-10;

            double[] r = MathOperations.Matrix_Mult(A, N, q);
            for (int i = 0; i < N; i++)
            {
                q[i] = 1;
                r[i] = b[i] - r[i];
            }
            double[] z = new double[N];
            r.CopyTo(z, 0);
            double[] p = MathOperations.Matrix_Mult(A, N, z);

            do
            {
                ++count;
                alpha = MathOperations.ScalarMult(p, r) / MathOperations.ScalarMult(p, p);
                for (int i = 0; i < N; i++)
                {
                    q[i] = q[i] + alpha * z[i];
                    r[i] = r[i] - alpha * p[i];
                }
                double[] Ar = MathOperations.Matrix_Mult(A, N, r);
                beta = - MathOperations.ScalarMult(p, Ar) / MathOperations.ScalarMult(p, p);
                
                for (int i = 0; i < N; i++)
                {
                    z[i] = r[i] + beta * z[i];
                    p[i] = Ar[i] + beta * p[i];
                }

                nev = MathOperations.ScalarMult(r, r);
            } while (count < maxiter && nev > eps);
            System.Console.WriteLine($"{count} iterations on los");
            System.Console.WriteLine($"{nev} - Невязка");
        }
        
        public ref double A(int i, int j)
        {
            if (i == j)
                return ref di[i];
            bool gguflag = false;
            if (i < j)
            {
                MathOperations.Swap(ref i, ref j);
                gguflag = true;
            }
            j--;

            for (int k = ig[i]; k < ig[i + 1]; k++)
            {
                if (jg[k] == j + 1)
                {
                    if (gguflag)
                        return ref ggu[k];
                    return ref ggl[k];
                }
            }

            return ref zero;
        }
    }
}