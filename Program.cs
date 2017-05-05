using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MkeProject
{
    enum Edge { Right, Top, Left, Bottom }

    struct Pair
    {
        public int First { get; set; }
        public int Second { get; set; }
        public Pair(int first, int second)
        {
            First = first;
            Second = second;
        }
    }

    class Program
    {
        static SLAE slae;
        static double lambda = 1, gamma = 1, beta = 1;
        static double[] x, y;

        static void Main(string[] args)
        {
            x = new double[] { 0, 1, 2, 3, 4, 5 };
            y = new double[] { 0, 1, 2, 3, 4, 5 };

            GeneratePortrait();

            //цикл по конечным элементам
            for (int i = 0; i < y.Length - 1; i++)
            {
                for (int j = 0; j < x.Length - 1; j++)
                {
                    int number = i * x.Length + j;
                    AddLocal(number);
                }
            }

            //нижняя граница
            for (int i = 0; i < x.Length - 1; i++)
            {
                kuslau1(i, i + 1);
            }

            //правая граница
            for (int i = x.Length - 1; i < x.Length * y.Length - 1; i = i + x.Length)
            {
                kuslau2(i, i + x.Length, Edge.Right);
            }

            //верхняя граница
            for (int i = x.Length * y.Length - x.Length; i < x.Length * y.Length - 1; i++)
            {
                kuslau3(i, i + 1, Edge.Top);
            }

            //левая граница
            for (int i = 0; i < x.Length * y.Length - x.Length; i = i + x.Length)
            {
                kuslau2(i, i + x.Length, Edge.Left);
            }

            slae.Calculate_LOS();

            for (int i = 0; i < x.Length * y.Length; i++)
            {
                System.Console.WriteLine();
                System.Console.Write($"{slae.q[i]}      {U(x[i%x.Length], y[i/x.Length])}");
            }
        }

        static void GeneratePortrait()
        {
            double[] ggl;
            int[] jg, ig;
            int number = 0;
            HashSet<Pair> pairs = new HashSet<Pair>();

            //цикл по конечным элементам
            for (int i = 0; i < y.Length - 1; i++)
            {
                for (int j = 0; j < x.Length - 1; j++)
                {
                    number = i * x.Length + j;

                    pairs.Add(new Pair(number + 1, number));
                    pairs.Add(new Pair(number + x.Length, number));
                    pairs.Add(new Pair(number + x.Length + 1, number));
                    pairs.Add(new Pair(number + x.Length + 1, number + 1));
                    pairs.Add(new Pair(number + x.Length + 1, number + x.Length));
                    pairs.Add(new Pair(number + x.Length, number + 1));
                }
            }

            List<Pair> bind = new List<Pair>(pairs);
            bind.SortPairs(x.Length * y.Length);
            ggl = new double[bind.Count];
            jg = new int[bind.Count];
            ig = new int[x.Length * y.Length + 1];

            int count = 2;
            for (int i = 0; i < bind.Count; i++)
            {
                jg[i] = bind[i].Second;
                if (i != 0)
                {
                    if (bind[i].First > bind[i - 1].First)
                    {
                        ig[count] = i;
                        count++;
                    }
                }
            }
            ig[x.Length * y.Length] = bind.Count;

            slae = new SLAE(x.Length * y.Length, ggl, ig, jg);
        }

        static void AddLocal(int number)
        {
            double x0, x1, y0, y1, hx, hy;
            x0 = x[number % x.Length];
            y0 = y[number / x.Length];
            x1 = x[number % x.Length + 1];
            y1 = y[number / x.Length + 1];
            hx = x1 - x0;
            hy = y1 - y0;

            double[,] G = Program.G(hx, hy);
            double[,] M = Program.M(hx, hy);

            slae.A(number, number) += G[0, 0] + M[0, 0];
            slae.A(number, number + 1) += G[0, 1] + M[0, 1];
            slae.A(number, number + x.Length) += G[0, 2] + M[0, 2];
            slae.A(number, number + x.Length + 1) += G[0, 3] + M[0, 3];
            slae.A(number + 1, number) += G[1, 0] + M[1, 0];
            slae.A(number + 1, number + 1) += G[1, 1] + M[1, 1];
            slae.A(number + 1, number + x.Length) += G[1, 2] + M[1, 2];
            slae.A(number + 1, number + x.Length + 1) += G[1, 3] + M[1, 3];
            slae.A(number + x.Length, number) += G[2, 0] + M[2, 0];
            slae.A(number + x.Length, number + 1) += G[2, 1] + M[2, 1];
            slae.A(number + x.Length, number + x.Length) += G[2, 2] + M[2, 2];
            slae.A(number + x.Length, number + x.Length + 1) += G[2, 3] + M[2, 3];
            slae.A(number + x.Length + 1, number) += G[3, 0] + M[3, 0];
            slae.A(number + x.Length + 1, number + 1) += G[3, 1] + M[3, 1];
            slae.A(number + x.Length + 1, number + x.Length) += G[3, 2] + M[3, 2];
            slae.A(number + x.Length + 1, number + x.Length + 1) += G[3, 3] + M[3, 3];

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    slae.b[number] += (M[0, 2 * i + j]) * f(x[number % x.Length + j], y[number / x.Length + i]);
                    slae.b[number + 1] += (M[1, 2 * i + j]) * f(x[number % x.Length + j], y[number / x.Length + i]);
                    slae.b[number + x.Length] += (M[2, 2 * i + j]) * f(x[number % x.Length + j], y[number / x.Length + i]);
                    slae.b[number + x.Length + 1] += (M[3, 2 * i + j]) * f(x[number % x.Length + j], y[number / x.Length + i]);
                }
            }
        }

        static double f(double x, double y)
        {
            return U(x, y) - 4;
        }

        static double U(double x, double y)
        {
            return x + y + x * y + x * x + y * y;
        }

        static double teta(double x, double y, Edge edge)
        {
            switch (edge)
            {
                case Edge.Right:
                    return 1 + y + 2 * x;
                case Edge.Top:
                    return 1 + x + 2 * y;
                case Edge.Left:
                    return -(1 + y + 2 * x);
                case Edge.Bottom:
                    return -(1 + x + 2 * y);
                default:
                    return 0;
            }
        }

        static double Ubeta(double x, double y, Edge edge)
        {
            return U(x, y) + teta(x, y, edge) / beta;
        }

        static void kuslau1(int number1, int number2)
        {
            double x1, x2, y1, y2;
            x1 = x[number1 % x.Length];
            y1 = y[number1 / x.Length];
            x2 = x[number2 % x.Length];
            y2 = y[number2 / x.Length];
            slae.A(number1, number1) = 1E20;
            slae.A(number2, number2) = 1E20;
            slae.b[number1] = U(x1, y1) * 1E20;
            slae.b[number2] = U(x2, y2) * 1E20;
        }

        static void kuslau2(int number1, int number2, Edge edge) 
        {
            double x1, x2, y1, y2, h;

            x1 = x[number1 % x.Length];
            y1 = y[number1 / x.Length];
            x2 = x[number2 % x.Length];
            y2 = y[number2 / x.Length];

            double teta1 = teta(x1, y1, edge);
            double teta2 = teta(x2, y2, edge);

            if (edge == Edge.Right || edge == Edge.Left)
            {
                h = y2 - y1;
            }
            else
            {
                h = x2 - x1;
            }

            slae.b[number1] += h / 6 * (2 * teta1 + teta2);
	        slae.b[number2] += h / 6 * (teta1 + 2 * teta2);
        }

        static void kuslau3(int number1, int number2, Edge edge)
        {
            double x1, y1, x2, y2, h;

            x1 = x[number1 % x.Length];
            y1 = y[number1 / x.Length];
            x2 = x[number2 % x.Length];
            y2 = y[number2 / x.Length];
            
            if (edge == Edge.Right || edge == Edge.Left)
            {
                h = y2 - y1;
            }
            else
            {
                h = x2 - x1;
            }

            double ubeta1 = Ubeta(x1, y1, edge);
            double ubeta2 = Ubeta(x2, y2, edge);

            double koef = beta * h / 6;

            slae.A(number1, number1) += 2 * koef;
            slae.A(number1, number2) += 1 * koef;
            slae.A(number2, number1) += 1 * koef;
            slae.A(number2, number2) += 2 * koef;

            slae.b[number1] += koef * (2 * ubeta1 + ubeta2);
            slae.b[number2] += koef * (ubeta1 + 2 * ubeta2);
        }

        static double[,] G(double hx, double hy)
        {
            double[,] G = new double[,]
            {
                { 2, -2, 1, -1 },
                { -2, 2, -1, 1 },
                { 1, -1, 2, -2 },
                { -1, 1, -2, 2 }
            };
            double[,] G2 = new double[,]
            {
                { 2, 1, -2, -1 },
                { 1, 2, -1, -2 },
                { -2, -1, 2, 1 },
                { -1, -2, 1, 2 }
            };
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    G[i, j] = G[i, j] * lambda * hy / 6 / hx + G2[i, j] * lambda * hx / 6 / hy;
                }
            }
            return G;
        }

        static double[,] M(double hx, double hy)
        {
            double[,] M = new double[,]
            {
                { 4, 2, 2, 1 },
                { 2, 4, 1, 2 },
                { 2, 1, 4, 2 },
                { 1, 2, 2, 4 }
            };
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    M[i, j] *= gamma * hx * hy / 36;
                }
            }
            return M;
        }
    }
}
