using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SolveSLE
{
    class Program
    {
        static void Main(string[] args)
        {
            Matrix<double> Aj = DenseMatrix.OfArray(new double[,] { { 10, -1,  2, 0 }, 
                                                                    { -1, 11, -1, 3 }, 
                                                                    {  2, -1, 10,-1 }, 
                                                                    {  0, 3 , -1, 8 }
                                                                  });
            Vector<double> Bj = DenseVector.OfArray(new double[] { 6, 25, -11, 15 });

            Matrix<double> As = DenseMatrix.OfArray(new double[,] { { 5, 2,-4 }, 
                                                                    { 2, 1,-2 }, 
                                                                    {-4,-2, 5,}
                                                                  });
            Vector<double> Bs = DenseVector.OfArray(new double[] { 1, 2, 3 });

            SleSolver sleSolver = new SleSolver();

            var A = sleSolver.generateDiagonallyDominantMatrix(4);
            Console.WriteLine(A);

            Console.WriteLine(Aj * sleSolver.JacobiMethod(Aj, Bj, 1e-16));
            Console.WriteLine(Aj * sleSolver.SeidelMethod(Aj, Bj, 1e-16));

            Console.WriteLine(As * sleSolver.JacobiMethod(As, Bs, 1e-16));
            Console.WriteLine(As * sleSolver.SeidelMethod(As, Bs, 1e-16));


        }
    }
}
