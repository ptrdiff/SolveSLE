using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace SolveSLE
{
    class SleSolver
    {
        public Vector<double> JacobiMethod(Matrix<double> A, Vector<double> B, double eps)
        {
            Vector<double> X = new DenseVector(A.RowCount);
            Vector<double> TempX = new DenseVector(A.RowCount);
            double norm = 0.0;
            int count = 0;
            do
            {
                for (int i = 0; i < A.RowCount; ++i)
                {
                    TempX[i] = B[i];
                    for (int g = 0; g < A.RowCount; ++g)
                    {
                        if (i != g)
                        {
                            TempX[i] -= A[i, g] * X[g];
                        }
                    }
                    TempX[i] /= A[i,i];
                }
                norm = Math.Abs(X[0] - TempX[0]);
                for (int h = 0; h < A.RowCount; ++h)
                {
                    if (Math.Abs(X[h] - TempX[h]) > norm)
                    {
                        norm = Math.Abs(X[h] - TempX[h]);
                    }
                    X[h] = TempX[h];
                }
                if (++count == 1000) break;
            } while (norm > eps);
            Console.WriteLine(count);
            return X;
        }
        public Vector<double> SeidelMethod(Matrix<double> A, Vector<double> B, double eps)
        {
            Vector<double> X = new DenseVector(A.RowCount);
            Vector<double> TempX = new DenseVector(A.RowCount);
            int count = 0;
            do
            {
                TempX = X.Clone();
                for (int i = 0; i < A.RowCount; ++i)
                {
                    double var = 0;
                    for (int j = 0; j < i; j++)
                    {
                        var += (A[i, j] * X[j]);
                    }
                    for (int j = i + 1; j < A.RowCount; ++j)
                    {
                        var += (A[i, j] * TempX[j]);
                    }
                    X[i] = (B[i] - var) / A[i,i];
                }
                if (++count == 1000) break;

            } while ((X - TempX).L2Norm() > eps);
            Console.WriteLine(count);
            return X;
        }
        public Matrix<double> generateDiagonallyDominantMatrix(int n)
        {
            Random random = new Random();
            Matrix<double> matrix = new DenseMatrix(n);

            for (var i = 0; i < matrix.RowCount; ++i)
            {
                matrix[i, i] = random.Next(-100, 101);
                var norm = Math.Abs((int)matrix[i, i]);
                for (var j = 0; j < matrix.ColumnCount; ++j)
                {
                    if (i != j)
                    {
                        matrix[i, j] = random.Next(-norm, norm);
                        norm = Math.Abs(norm - Math.Abs((int)matrix[i, j]));
                    }
                }
            }
            return matrix;
        }
    }
}
