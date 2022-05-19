/*
* Matrix library for Unity with projection transformation functionality.
* Copyright (c)2022 Takumi Mori, Osaka Electro-Communication University
* 
* The MIT License
* 
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
* OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/

namespace Unity.MathMatrix
{
    using System.Collections.Generic;

    public class Matrix
    {
        protected double[,] elements;

        public Matrix(int rows, int cols)
        {
            this.elements = new double[rows, cols];
        }

        public Matrix(int rows, int cols, float[] elements)
        {
            this.elements = new double[rows, cols];
            var k = 0;
            for (var i = 0; i < cols; i++)
            {
                for (var j = 0; j < rows; j++)
                {
                    this.elements[j, i] = (double)elements[k++];
                }
            }
        }

        public Matrix(int rows, int cols, double[] elements)
        {
            this.elements = new double[rows, cols];
            var k = 0;
            for (var i = 0; i < cols; i++)
            {
                for (var j = 0; j < rows; j++)
                {
                    this.elements[j, i] = elements[k++];
                }
            }
        }

        public Matrix(double[,] elements)
        {
            this.elements = elements;
        }

        public int Rows() { return elements.GetLength(0); }
        public int Cols() { return elements.GetLength(1); }

        public void Set(int x, int y, double value)
        {
            elements[x, y] = value;
        }

        public double Get(int x, int y)
        {
            return elements[x, y];
        }
        //　転置
        public Matrix T()
        {
            var _t = new Matrix(elements.GetLength(1), elements.GetLength(0));

            for (var i = 0; i < elements.GetLength(1); i++)
            {
                for (var j = 0; j < elements.GetLength(0); j++)
                {
                    _t.elements[i, j] = elements[j, i];
                }
            }
            return _t;
        }

        //　逆行列

        public Matrix Inv()
        {
            var n = elements.GetLength(0);
            var m = elements.GetLength(1);

            Matrix A = new Matrix(n, n, ToArray());
            Matrix invA = new Matrix(n, m);

            if (n == m)
            {
                int max;
                double tmp;

                for (int j = 0; j < n; j++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        invA.elements[j, i] = (i == j) ? 1 : 0;
                    }
                }

                for (int k = 0; k < n; k++)
                {
                    max = k;
                    for (int j = k + 1; j < n; j++)
                    {
                        if (System.Math.Abs(A.elements[j, k]) > System.Math.Abs(A.elements[max, k]))
                        {
                            max = j;
                        }
                    }

                    if (max != k)
                    {
                        for (int i = 0; i < n; i++)
                        {
                            // 入力行列側
                            tmp = A.elements[max, i];
                            A.elements[max, i] = A.elements[k, i];
                            A.elements[k, i] = tmp;
                            // 単位行列側
                            tmp = invA.elements[max, i];
                            invA.elements[max, i] = invA.elements[k, i];
                            invA.elements[k, i] = tmp;
                        }
                    }

                    tmp = A.elements[k, k];

                    for (int i = 0; i < n; i++)
                    {
                        A.elements[k, i] /= tmp;
                        invA.elements[k, i] /= tmp;
                    }

                    for (int j = 0; j < n; j++)
                    {
                        if (j != k)
                        {
                            tmp = A.elements[j, k] / A.elements[k, k];
                            for (int i = 0; i < n; i++)
                            {
                                A.elements[j, i] = A.elements[j, i] - A.elements[k, i] * tmp;
                                invA.elements[j, i] = invA.elements[j, i] - invA.elements[k, i] * tmp;
                            }
                        }
                    }
                }

                return invA;
            }

            return null;
        }

        // 行列の要素ごとの乗算
        public Matrix Mul(Matrix mat)
        {
            var _mul = new Matrix(elements.GetLength(0), mat.elements.GetLength(1));

            for (var i = 0; i < elements.GetLength(0); i++)
            {
                for (var j = 0; j < mat.elements.GetLength(1); j++)
                {
                    for (var k = 0; k < elements.GetLength(1); k++)
                    {
                        _mul.elements[i, j] += elements[i, k] * mat.elements[k, j];
                    }
                }
            }
            return _mul;
        }

        public float[] ToArray()
        {
            var array = new List<float>();

            for (var i = 0; i < elements.GetLength(1); i++)
            {
                for (var j = 0; j < elements.GetLength(0); j++)
                {
                    array.Add((float)elements[i, j]);
                }
            }
            return array.ToArray();
        }

        public double[] ToArrayDouble()
        {
            var array = new double[elements.Length];
            var cnt = 0;
            for (var i = 0; i < elements.GetLength(1); i++)
            {
                for (var j = 0; j < elements.GetLength(0); j++)
                {
                    array[cnt] = elements[i, j];
                    cnt++;
                }
            }
            return array;
        }

        public float[] ToArrayOpenCVforUnity(int Cvtype)
        {
            var array = new List<float>();
            array.Add(elements.GetLength(0));
            array.Add(elements.GetLength(1));
            array.Add(Cvtype);

            for (var i = 0; i < elements.GetLength(1); i++)
            {
                for (var j = 0; j < elements.GetLength(0); j++)
                {
                    array.Add((float)elements[j, i]);
                }
            }
            return array.ToArray();
        }

        // LU分解
        private Matrix LUDecomp()
        {
            if (elements.GetLength(0) == elements.GetLength(1))
            {
                var N = elements.GetLength(0);
                var A = new Matrix(N, N, ToArray());

                for (var pivot = 0; pivot < N - 1; pivot++)
                {
                    for (var row = pivot + 1; row < N; row++)
                    {
                        var s = A.elements[row, pivot] / A.elements[pivot, pivot];
                        for (var col = pivot; col < N; col++) A.elements[row, col] -= A.elements[pivot, col] * s;
                        A.elements[row, pivot] = s;
                    }
                }
                return A;
            }
            return null;
        }

        // 前進代入
        private void ForwardSub(double[,]a, ref double[]y,ref double[]b)
        {
            for (var row = 0; row < a.GetLength(0); row++)
            {
                for (var col = 0; col < row; col++) b[row] -= a[row,col] * y[col];
                y[row] = b[row];
            }
        }

        // 後退代入
        private void BackSub(double[,] a,ref double[] x,ref double[] y)
        {
            for (var row = a.GetLength(0) -1; row >= 0; row--)
            {
                for (var col = a.GetLength(0) - 1; col > row; col--)
                {
                    y[row] -= a[row, col] * x[col];
                }
                x[row] = y[row] / a[row, row];
            }
        }

        // 正規化
        private void Normarize(ref double[]x, int N)
        {
            double s = 0;
            for (var i = 0; i < N; i++) s += x[i] * x[i];
            s = System.Math.Sqrt(s);

            for (var i = 0; i < N; i++) x[i] /= s;

        }

        // 逆べき乗法 固有ベクトル・固有値の算出
        public bool InversPower(ref double[]x0, ref double lamda)
        {
            

            if (elements.GetLength(0) == elements.GetLength(1) && elements.GetLength(0) == x0.Length)
            {
                var a = LUDecomp();
                var N = elements.GetLength(0);

                // 正規化
                Normarize(ref x0, x0.Length);

                double e0 = 0;

                for (var i = 0; i < N; i++) e0 += x0[i];

                for (var k = 1; k <= 100; k++)
                {

                    // Ly = b から yを求める
                    var b = new double[N];
                    var y = new double[N];
                    for (var i = 0; i < N; i++) b[i] = x0[i];
                    ForwardSub(a.elements, ref y, ref b);

                    // Ux = y から x を求める (後退代入)
                    var x1 = new double[N];
                    BackSub(a.elements, ref x1,ref y);

                    // 内積
                    double p0 = 0;
                    double p1 = 0;

                    for (var i = 0; i < N; i++)
                    {
                        p0 += x1[i] * x1[i];
                        p1 += x1[i] * x0[i];
                    }

                    // 固有値
                    lamda = p1 / p0;

                    Normarize(ref x1, N);

                    double e1 = 0;
                    for (var i = 0; i < N; i++) e1 += x1[i];
                    if (System.Math.Abs((double)(e0 - e1)) < 0.00000000001) break;

                    for (var i = 0; i < N; i++) x0[i] = x1[i];
                    e0 = e1;
;               }
                return true;
            }
            return false;
        }

        

        //　以下Unity専用
        #if DEBUG
        public void DebugLog()
        {
            string output = "";

            for (var i = 0; i < elements.GetLength(1); i++)
            {
                for (var j = 0; j < elements.GetLength(0); j++)
                {
                    output += elements[i, j] + ",";
                }
                output += "\n";
            }

            UnityEngine.Debug.Log(output);
        }
        #endif

        public static Matrix DLT(UnityEngine.Vector2[] beforPoints, UnityEngine.Vector2[] afterPoints, bool useOpenCV)
        {

            if (beforPoints.Length == afterPoints.Length)
            {
                var N = beforPoints.Length;
                // 4点から行列導出　n = beforPoints.Length 9x2n
                var mat = new Matrix(9, 2 * N);

                var point = 0;
                for (var y = 0; y < N * 2; y += 2)
                {
                    mat.Set(0, y, 0);
                    mat.Set(1, y, 0);
                    mat.Set(2, y, 0);

                    mat.Set(3, y, beforPoints[point].x);
                    mat.Set(4, y, beforPoints[point].y);
                    mat.Set(5, y, 1);

                    mat.Set(6, y, afterPoints[point].y * beforPoints[point].x);
                    mat.Set(7, y, afterPoints[point].y * beforPoints[point].y);
                    mat.Set(8, y, afterPoints[point].y * 1);

                    mat.Set(0, y + 1, beforPoints[point].x);
                    mat.Set(1, y + 1, beforPoints[point].y);
                    mat.Set(2, y + 1, 1);

                    mat.Set(3, y + 1, 0);
                    mat.Set(4, y + 1, 0);
                    mat.Set(5, y + 1, 0);

                    mat.Set(6, y + 1, -afterPoints[point].x * beforPoints[point].x);
                    mat.Set(7, y + 1, -afterPoints[point].x * beforPoints[point].y);
                    mat.Set(8, y + 1, -afterPoints[point].x * 1);

                    point++;
                }

                // AA^T = VS^2 * V^T
                var n = mat.Mul(mat.T());

                // 固有値分解
                // 正方ベクトルでしか固有値分解できない
                var x0 = new double[n.Rows()];
                x0[0] = 1;
                double lamda = 0;
                // 逆べき乗法で固有値と固有ベクトルを算出
                n.InversPower(ref x0, ref lamda);

                // OpenCVの形式へ直す
                if (useOpenCV == true)
                {
                    var invNum = 1 / x0[x0.Length - 1];
                    for (var i = 0; i < n.Rows(); i++) x0[i] *= invNum;
                }

                var H = new Matrix((int)System.Math.Sqrt(n.Rows()), (int)System.Math.Sqrt(n.Rows()));
                var cnt = 0;

                for (var i = 0; i < H.Cols(); i++)
                {
                    for (var j = 0; j < H.Rows(); j++)
                    {
                        H.Set(i, j, x0[cnt++]);
                    }
                }
                return H;
            }
            else
            {
                UnityEngine.Debug.LogError("Different Points count! Check beforPoints count & afterPoints count");
                return null;
            }
            

        }
    }

    
}