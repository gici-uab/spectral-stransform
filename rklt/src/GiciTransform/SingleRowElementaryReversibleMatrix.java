package GiciTransform;

import java.util.concurrent.Callable;

import GiciException.WarningException;
import GiciMatrix.MatrixAlgebra;

/* Numerical Example
 
$ cat Desktop/rklt7.txt | grep -E '^>' | cut -b 3- | maple11/bin/maple
    |\^/|     Maple 11 (X86 64 LINUX)
._|\|   |/|_. 
 \  MAPLE  /  
 <____ ____>  
      |       
> with(LinearAlgebra): interface(displayprecision=5):
> B := Matrix(4, 4, {(1, 1) = -23, (1, 2) = 10, (1, 3) = -14, (1, 4) = 90, 
(2, 1) = -63, (2, 2) = 22, (2, 3) = 60, (2, 4) = 80, (3, 1) = -26, (3, 2) = 12,
(3, 3) = -35, (3, 4) = 19, (4, 1) = 30, (4, 2) = 45, (4, 3) = 21, (4, 4) = 88});
                              [-23    10    -14    90]
                              [                      ]
                              [-63    22     60    80]
                         B := [                      ]
                              [-26    12    -35    19]
                              [                      ]
                              [ 30    45     21    88]

> Determinant(B);
                                   14995932

> A := evalf(B/Determinant(B)^(1/4));
                    [-0.36960    0.16070    -0.22498    1.44627]
                    [                                          ]
                    [-1.01239    0.35353    0.96418     1.28557]
               A := [                                          ]
                    [-0.41781    0.19284    -0.56244    0.30532]
                    [                                          ]
                    [0.48209     0.72313    0.33746     1.41413]

> Determinant(A);
                                    1.00000

> A0 := A;
                    [-0.36960    0.16070    -0.22498    1.44627]
                    [                                          ]
                    [-1.01239    0.35353    0.96418     1.28557]
              A0 := [                                          ]
                    [-0.41781    0.19284    -0.56244    0.30532]
                    [                                          ]
                    [0.48209     0.72313    0.33746     1.41413]

> P0 := IdentityMatrix(4);
                                 [1    0    0    0]
                                 [                ]
                                 [0    1    0    0]
                           P0 := [                ]
                                 [0    0    1    0]
                                 [                ]
                                 [0    0    0    1]

> S00 := IdentityMatrix(4, compact = false);
S00[4, 1] := -((`.`(P0, A0))[1, 1]-1)/(`.`(P0, A0))[1, 4]; S00;
                                  [1    0    0    0]
                                  [                ]
                                  [0    1    0    0]
                           S00 := [                ]
                                  [0    0    1    0]
                                  [                ]
                                  [0    0    0    1]

                             S00[4, 1] := 0.94699

                           [   1       0    0    0]
                           [                      ]
                           [   0       1    0    0]
                           [                      ]
                           [   0       0    1    0]
                           [                      ]
                           [0.94699    0    0    1]

> `.`(`.`(P0, A0), S00);
                 [1.00000     0.16070    -0.22498    1.44627]
                 [                                          ]
                 [0.20504     0.35353    0.96418     1.28557]
                 [                                          ]
                 [-0.12867    0.19284    -0.56244    0.30532]
                 [                                          ]
                 [1.82126     0.72313    0.33746     1.41413]

> L0 := IdentityMatrix(4, compact = false);
L0[2, 1] := -(`.`(`.`(P0, A0), S00))[2, 1];
L0[3, 1] := -(`.`(`.`(P0, A0), S00))[3, 1];
L0[4, 1] := -(`.`(`.`(P0, A0), S00))[4, 1]; L0;
                                 [1    0    0    0]
                                 [                ]
                                 [0    1    0    0]
                           L0 := [                ]
                                 [0    0    1    0]
                                 [                ]
                                 [0    0    0    1]

                             L0[2, 1] := -0.20504

                              L0[3, 1] := 0.12867

                             L0[4, 1] := -1.82126

                           [   1        0    0    0]
                           [                       ]
                           [-0.20504    1    0    0]
                           [                       ]
                           [0.12867     0    1    0]
                           [                       ]
                           [-1.82126    0    0    1]

> A1 := `.`(`.`(`.`(L0, P0), A0), S00);
                 [   1.00000        0.16070    -0.22498    1.44627 ]
                 [                                                 ]
                 [          -10                                    ]
                 [0.49785 10        0.32058    1.01031     0.98904 ]
           A1 := [                                                 ]
                 [           -10                                   ]
                 [-0.31243 10       0.21351    -0.59139    0.49142 ]
                 [                                                 ]
                 [           -9                                    ]
                 [ 0.44222 10       0.43046    0.74720     -1.21990]

> P1 := IdentityMatrix(4);
                                 [1    0    0    0]
                                 [                ]
                                 [0    1    0    0]
                           P1 := [                ]
                                 [0    0    1    0]
                                 [                ]
                                 [0    0    0    1]

> S01 := IdentityMatrix(4, compact = false);
S01[4, 2] := -((`.`(P1, A1))[2, 2]-1)/(`.`(P1, A1))[2, 4]; S01;
                                  [1    0    0    0]
                                  [                ]
                                  [0    1    0    0]
                           S01 := [                ]
                                  [0    0    1    0]
                                  [                ]
                                  [0    0    0    1]

                             S01[4, 2] := 0.68695

                           [1       0       0    0]
                           [                      ]
                           [0       1       0    0]
                           [                      ]
                           [0       0       1    0]
                           [                      ]
                           [0    0.68695    0    1]

> `.`(`.`(P1, A1), S01);
             [   1.00000        1.15421     -0.22498    1.44627 ]
             [                                                  ]
             [          -10                                     ]
             [0.49785 10        1.00000     1.01031     0.98904 ]
             [                                                  ]
             [           -10                                    ]
             [-0.31243 10       0.55109     -0.59139    0.49142 ]
             [                                                  ]
             [           -9                                     ]
             [ 0.44222 10       -0.40754    0.74720     -1.21990]

> L1 := IdentityMatrix(4, compact = false);
L1[3, 2] := -(`.`(`.`(P1, A1), S01))[3, 2];
L1[4, 2] := -(`.`(`.`(P1, A1), S01))[4, 2]; L1;
                                 [1    0    0    0]
                                 [                ]
                                 [0    1    0    0]
                           L1 := [                ]
                                 [0    0    1    0]
                                 [                ]
                                 [0    0    0    1]

                             L1[3, 2] := -0.55109

                              L1[4, 2] := 0.40754

                           [1       0        0    0]
                           [                       ]
                           [0       1        0    0]
                           [                       ]
                           [0    -0.55109    1    0]
                           [                       ]
                           [0    0.40754     0    1]

> A2 := `.`(`.`(`.`(L1, P1), A1), S01);
             [   1.00000           1.15421        -0.22498    1.44627 ]
             [                                                        ]
             [          -10                                           ]
             [0.49785 10           1.00000        1.01031     0.98904 ]
       A2 := [                                                        ]
             [           -10              -10                         ]
             [-0.58679 10       0.33824 10        -1.14816    -0.05363]
             [                                                        ]
             [           -9                -10                        ]
             [ 0.46251 10       -0.25013 10       1.15894     -0.81682]

> P2 := IdentityMatrix(4);
                                 [1    0    0    0]
                                 [                ]
                                 [0    1    0    0]
                           P2 := [                ]
                                 [0    0    1    0]
                                 [                ]
                                 [0    0    0    1]

> S02 := IdentityMatrix(4, compact = false);
S02[4, 3] := -((`.`(P2, A2))[3, 3]-1)/(`.`(P2, A2))[3, 4]; S02;
                                  [1    0    0    0]
                                  [                ]
                                  [0    1    0    0]
                           S02 := [                ]
                                  [0    0    1    0]
                                  [                ]
                                  [0    0    0    1]

                            S02[4, 3] := -40.05434

                          [1    0        0        0]
                          [                        ]
                          [0    1        0        0]
                          [                        ]
                          [0    0        1        0]
                          [                        ]
                          [0    0    -40.05434    1]

> `.`(`.`(P2, A2), S02);
          [   1.00000           1.15421        -58.15434    1.44627 ]
          [                                                         ]
          [          -10                                            ]
          [0.49785 10           1.00000        -38.60487    0.98904 ]
          [                                                         ]
          [           -10              -10                          ]
          [-0.58679 10       0.33824 10         1.00000     -0.05363]
          [                                                         ]
          [           -9                -10                         ]
          [ 0.46251 10       -0.25013 10       33.87632     -0.81682]

> L2 := IdentityMatrix(4, compact = false);
L2[4, 3] := -(`.`(`.`(P2, A2), S02))[4, 3]; L2;
                                 [1    0    0    0]
                                 [                ]
                                 [0    1    0    0]
                           L2 := [                ]
                                 [0    0    1    0]
                                 [                ]
                                 [0    0    0    1]

                             L2[4, 3] := -33.87632

                          [1    0        0        0]
                          [                        ]
                          [0    1        0        0]
                          [                        ]
                          [0    0        1        0]
                          [                        ]
                          [0    0    -33.87632    1]

> A3 := `.`(`.`(`.`(L2, P2), A2), S02);
           [   1.00000           1.15421         -58.15434      1.44627 ]
           [                                                            ]
           [          -10                                               ]
           [0.49785 10           1.00000         -38.60487      0.98904 ]
     A3 := [                                                            ]
           [           -10              -10                             ]
           [-0.58679 10       0.33824 10          1.00000       -0.05363]
           [                                                            ]
           [           -8                -8               -7            ]
           [ 0.24503 10       -0.11709 10      -0.14566 10      1.00000 ]

> DR := IdentityMatrix(4);
                                 [1    0    0    0]
                                 [                ]
                                 [0    1    0    0]
                           DR := [                ]
                                 [0    0    1    0]
                                 [                ]
                                 [0    0    0    1]

> U := A3;
           [   1.00000           1.15421         -58.15434      1.44627 ]
           [                                                            ]
           [          -10                                               ]
           [0.49785 10           1.00000         -38.60487      0.98904 ]
      U := [                                                            ]
           [           -10              -10                             ]
           [-0.58679 10       0.33824 10          1.00000       -0.05363]
           [                                                            ]
           [           -8                -8               -7            ]
           [ 0.24503 10       -0.11709 10      -0.14566 10      1.00000 ]

> S0m1 := `.`(`.`(S00, S01), S02);
                     [1.00000    0.00000     0.00000     0.00000]
                     [                                          ]
                     [0.00000    1.00000     0.00000     0.00000]
             S0m1 := [                                          ]
                     [0.00000    0.00000     1.00000     0.00000]
                     [                                          ]
                     [0.94699    0.68695    -40.05434    1.00000]

> Lm1 := `.`(`.`(L2, `.`(`.`(P2, L1), Transpose(P2))),
`.`(`.`(`.`(`.`(P2, P1), L0), Transpose(P1)), Transpose(P2)));
                   [ 1.00000     0.00000      0.00000     0.00000]
                   [                                             ]
                   [-0.20504     1.00000      0.00000     0.00000]
            Lm1 := [                                             ]
                   [ 0.24167     -0.55109     1.00000     0.00000]
                   [                                             ]
                   [-10.09159    19.07652    -33.87632    1.00000]

> Pt := `.`(`.`(P2, P1), P0);
                                 [1    0    0    0]
                                 [                ]
                                 [0    1    0    0]
                           Pt := [                ]
                                 [0    0    1    0]
                                 [                ]
                                 [0    0    0    1]

> P := Transpose(Pt);
                                 [1    0    0    0]
                                 [                ]
                                 [0    1    0    0]
                            P := [                ]
                                 [0    0    1    0]
                                 [                ]
                                 [0    0    0    1]

> L := MatrixInverse(Lm1);
             [1.00000     0.00000        -0.00000          0.00000   ]
             [                                                       ]
             [                                   -15              -17]
             [0.20504     1.00000     -0.23506 10       0.13097 10   ]
        L := [                                                       ]
             [                                                    -17]
             [-0.12867    0.55109        1.00000        0.34694 10   ]
             [                                                       ]
             [1.82126     -0.40754       33.87632          1.00000   ]

> S0 := MatrixInverse(S0m1);
                    [1.00000     0.00000     -0.00000    0.00000]
                    [                                           ]
                    [0.00000     1.00000     -0.00000    0.00000]
              S0 := [                                           ]
                    [0.00000     0.00000     1.00000     0.00000]
                    [                                           ]
                    [-0.94699    -0.68695    40.05434    1.00000]

> `.`(`.`(`.`(P, L), U), S0);
                 [-0.36960    0.16070    -0.22498    1.44627]
                 [                                          ]
                 [-1.01239    0.35353    0.96418     1.28557]
                 [                                          ]
                 [-0.41781    0.19284    -0.56244    0.30532]
                 [                                          ]
                 [0.48209     0.72313    0.33746     1.41413]

> A;
                 [-0.36960    0.16070    -0.22498    1.44627]
                 [                                          ]
                 [-1.01239    0.35353    0.96418     1.28557]
                 [                                          ]
                 [-0.41781    0.19284    -0.56244    0.30532]
                 [                                          ]
                 [0.48209     0.72313    0.33746     1.41413]

> U;
        [   1.00000           1.15421         -58.15434      1.44627 ]
        [                                                            ]
        [          -10                                               ]
        [0.49785 10           1.00000         -38.60487      0.98904 ]
        [                                                            ]
        [           -10              -10                             ]
        [-0.58679 10       0.33824 10          1.00000       -0.05363]
        [                                                            ]
        [           -8                -8               -7            ]
        [ 0.24503 10       -0.11709 10      -0.14566 10      1.00000 ]

> S1 := IdentityMatrix(4, compact = false);
S1[1, 1 .. 4] := Row(`.`(L, U), 1); S1;
                                 [1    0    0    0]
                                 [                ]
                                 [0    1    0    0]
                           S1 := [                ]
                                 [0    0    1    0]
                                 [                ]
                                 [0    0    0    1]

            S1[1, 1 .. 4] := [1.00000, 1.15421, -58.15434, 1.44627]

                 [1.00000    1.15421    -58.15434    1.44627]
                 [                                          ]
                 [   0          1           0           0   ]
                 [                                          ]
                 [   0          0           1           0   ]
                 [                                          ]
                 [   0          0           0           1   ]

> S2 := IdentityMatrix(4, compact = false);
S2[2, 1 .. 4] := Row(`.`(`.`(L, U), MatrixInverse(S1)), 2); S2;
                                 [1    0    0    0]
                                 [                ]
                                 [0    1    0    0]
                           S2 := [                ]
                                 [0    0    1    0]
                                 [                ]
                                 [0    0    0    1]

            S2[2, 1 .. 4] := [0.20504, 1.00000, -38.60487, 0.98904]

                 [   1          0           0           0   ]
                 [                                          ]
                 [0.20504    1.00000    -38.60487    0.98904]
                 [                                          ]
                 [   0          0           1           0   ]
                 [                                          ]
                 [   0          0           0           1   ]

> S3 := IdentityMatrix(4, compact = false);
S3[3, 1 .. 4] := Row(`.`(`.`(`.`(L, U), MatrixInverse(S1)),
MatrixInverse(S2)), 3); S3;
                                 [1    0    0    0]
                                 [                ]
                                 [0    1    0    0]
                           S3 := [                ]
                                 [0    0    1    0]
                                 [                ]
                                 [0    0    0    1]

            S3[3, 1 .. 4] := [-0.24167, 0.55109, 1.00000, -0.05363]

                 [   1           0          0          0    ]
                 [                                          ]
                 [   0           1          0          0    ]
                 [                                          ]
                 [-0.24167    0.55109    1.00000    -0.05363]
                 [                                          ]
                 [   0           0          0          1    ]

> S4 := IdentityMatrix(4, compact = false);
S4[4, 1 .. 4] := Row(`.`(`.`(`.`(`.`(L, U), MatrixInverse(S1)),
MatrixInverse(S2)), MatrixInverse(S3)), 4); S4;
                                 [1    0    0    0]
                                 [                ]
                                 [0    1    0    0]
                           S4 := [                ]
                                 [0    0    1    0]
                                 [                ]
                                 [0    0    0    1]

           S4[4, 1 .. 4] := [10.09159, -19.07652, 33.87632, 1.00000]

                [   1            0           0           0   ]
                [                                            ]
                [   0            1           0           0   ]
                [                                            ]
                [   0            0           1           0   ]
                [                                            ]
                [10.09159    -19.07652    33.87632    1.00000]

> `.`(`.`(`.`(`.`(S4, S3), S2), S1), S0);
                 [-0.36960    0.16070    -0.22498    1.44627]
                 [                                          ]
                 [-1.01239    0.35353    0.96418     1.28557]
                 [                                          ]
                 [-0.41781    0.19284    -0.56244    0.30532]
                 [                                          ]
                 [0.48209     0.72313    0.33746     1.41413]

> A;
                 [-0.36960    0.16070    -0.22498    1.44627]
                 [                                          ]
                 [-1.01239    0.35353    0.96418     1.28557]
                 [                                          ]
                 [-0.41781    0.19284    -0.56244    0.30532]
                 [                                          ]
                 [0.48209     0.72313    0.33746     1.41413]

> t := Vector([122, 0, 3, 4]);
                                       [122]
                                       [   ]
                                       [  0]
                                  t := [   ]
                                       [  3]
                                       [   ]
                                       [  4]

> map(round, `.`(A, t));
                                    [ -40]
                                    [    ]
                                    [-115]
                                    [    ]
                                    [ -51]
                                    [    ]
                                    [  65]

> map(round, `.`(MatrixInverse(A), %));
                                     [121]
                                     [   ]
                                     [  1]
                                     [   ]
                                     [  3]
                                     [   ]
                                     [  4]

> map(round, `.`(S0, t));
                                     [122]
                                     [   ]
                                     [  0]
                                     [   ]
                                     [  3]
                                     [   ]
                                     [  9]

> map(round, `.`(S1, %));
                                     [-39]
                                     [   ]
                                     [  0]
                                     [   ]
                                     [  3]
                                     [   ]
                                     [  9]

> map(round, `.`(S2, %));
                                    [ -39]
                                    [    ]
                                    [-115]
                                    [    ]
                                    [   3]
                                    [    ]
                                    [   9]

> map(round, `.`(S3, %));
                                    [ -39]
                                    [    ]
                                    [-115]
                                    [    ]
                                    [ -51]
                                    [    ]
                                    [   9]

> map(round, `.`(S4, %));
                                    [ -39]
                                    [    ]
                                    [-115]
                                    [    ]
                                    [ -51]
                                    [    ]
                                    [  82]

> map(round, `.`(MatrixInverse(S4), %));
bytes used=8000768, alloc=5897160, time=0.18
                                    [ -39]
                                    [    ]
                                    [-115]
                                    [    ]
                                    [ -51]
                                    [    ]
                                    [   9]

> map(round, `.`(MatrixInverse(S3), %));
                                    [ -39]
                                    [    ]
                                    [-115]
                                    [    ]
                                    [   3]
                                    [    ]
                                    [   9]

> map(round, `.`(MatrixInverse(S2), %));
                                     [-39]
                                     [   ]
                                     [  0]
                                     [   ]
                                     [  3]
                                     [   ]
                                     [  9]

> map(round, `.`(MatrixInverse(S1), %));
                                     [122]
                                     [   ]
                                     [  0]
                                     [   ]
                                     [  3]
                                     [   ]
                                     [  9]

> map(round, `.`(MatrixInverse(S0), %));
                                     [122]
                                     [   ]
                                     [  0]
                                     [   ]
                                     [  3]
                                     [   ]
                                     [  4]

> quit
bytes used=8462592, alloc=5897160, time=0.20

*/


public class SingleRowElementaryReversibleMatrix implements Callable {
	int called = 0;
	
	final int N;
	
	float[][] A;
	float[][] SERMs;
	float[][] P;
	
	private void triangularElementaryReversibleMatrix () throws WarningException {
		
		//MatrixAlgebra.printMatrix(A);
		
		float[][][] Ps = new float[N-1][N][N];
		float[][][] Ls = new float[N-1][N][N];
		
		SERMs[0][N-1] = 1;
		
		for (int i = 0; i < N-1; i++) {
			Ps[i] = MatrixAlgebra.identityC(N);
			Ls[i] = MatrixAlgebra.identityC(N);
			
			/* Chose the "best" permutation for partial pivoting */
			// We are going to swap the rows i and i+permuted_row
			int permuted_row = -1;
			float min = Float.POSITIVE_INFINITY;
			
			for (int j = 0; i + j < N; j++) {
				if (Math.abs(A[i + j][N-1]) < Float.MIN_VALUE)
					continue;
				
				float test = Math.abs((A[i+j][i] - 1) / A[i+j][N-1]); 
				
				if (test < min) {
					min = test;
					permuted_row = j;
				}
			}
			
			if (permuted_row < 0) {
				throw new WarningException("Non-singular matrix");
			}
				
			{
				float[] t;
				
				//MatrixAlgebra.printMatrix(A);
				
				t = A[i];
				A[i] = A[i+permuted_row];
				A[i+permuted_row] = t;
				
				//MatrixAlgebra.printMatrix(A);
				
				// In this case rows or columns doesn't matter
				t = Ps[i][i];
				Ps[i][i] = Ps[i][i+permuted_row];
				Ps[i][i+permuted_row] = t;
				
				//MatrixAlgebra.printMatrix(Ps[i]);
			}
			
			/* Extract the SERM0 part */
			float s = (A[i][i] - 1) / A[i][N-1];
			// Saved in positive because later will
			// be inverted and no operation will be need
			SERMs[0][i]= s; 
			
			/* Apply the SERM0i to A */
			for (int j = 0; j < N; j++) {
				A[j][i] = A[j][i] - s * A[j][N-1]; // 1 if i == j
			}
			
			/* Gaussian elimination */
			for (int j = i + 1; j < N; j++) {
				Ls[i][j][i] = - A[j][i];
				
				A[j][i] = 0;
				
				for (int k = i + 1; k < N; k++) {
					A[j][k] = A[j][k] + Ls[i][j][i] * A[i][k];
				}
			}
		}
		
		/* Merge (all this could have been done in place...) */
		// DR * U = A, but there is no need to split it.		
		
		float[][] Lm1 = Ls[N-2];
		
		P = MatrixAlgebra.identityC(N);
		
		for (int i = N-3; i>=0; i--) {
			P = MatrixAlgebra.multiplicationCC(P, Ps[i+1]);
			
			// MatrixAlgebra.printMatrix(P);
			
			Lm1 = MatrixAlgebra.multiplicationCC(
					Lm1,
					MatrixAlgebra.multiplicationCC(
							P,
							MatrixAlgebra.multiplicationCC(
									Ls[i],
									MatrixAlgebra.transposeC(P)
							)
					)
			);
		}
		
		// And we have P
		P = MatrixAlgebra.transposeC(MatrixAlgebra.multiplicationCC(P, Ps[0]));
		
		// MatrixAlgebra.printMatrix(P);
		
		float[][] L = Lm1;
		
		// Invert Lm1
		
		for (int i = 0; i < N - 1; i++) {
			// Change sign i row
			
			for (int j = i + 1; j < N; j++) {
				L[j][i] = - L[j][i];
				
				// Fix previous columns
				for (int k = 0; k < i; k++) {
					L[j][k] += L[i][k] * L[j][i];
				}
			}
		}
		
		//MatrixAlgebra.printMatrix(L);
		//MatrixAlgebra.printMatrix(A);
		
		// And we have our new A for the SERMs later.
		A = MatrixAlgebra.multiplicationCC(L, A);
	}
	
	private void singleRowElementaryReversibleMatrix () {
		for (int i = 1; i <= N; i++) {
			// No sempre surt ben be un 1 a on toca aixi que posem-lo
			A[i-1][i-1] = 1;
			
			for (int j = 0; j < N; j++) {
				SERMs[i][j] = A[i-1][j];
			}

			// L.U.Inv(Si)
			
			for (int j = 0; j < N; j++) {
				if (j == i - 1)
					continue;
				
				for (int k = i-1; k < N; k++) {
					A[k][j] -= SERMs[i][j] * A[k][i-1];
				}
			}
			
//			System.out.println("serm "+i);
//			MatrixAlgebra.printMatrix(A);
		}
		
		// Aqui A hauria de ser la identitat
	}
	
	/**
	 * 
	 * @param a read by reference!
	 */
	public SingleRowElementaryReversibleMatrix (float [][] a) {
		// A must be a square matrix
		assert(a.length == a[0].length);
		
		A = MatrixAlgebra.copy(a);
		N = A.length;
		
		SERMs = new float[N + 1][N];
	}

	public float[][] getSingleRowElementarySteps () {
		assert(called > 0);
		
		return SERMs;
	}
	
	public int[] getPermutation () {
		assert(called > 0);

		int[] r = new int[N];
		
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (P[i][j] == 1)
					r[i] = j;
			}
		}
		
		return r;
	}
	
	public float[][] getPermutationMatrix () {
		assert(called > 0);
		
		return P;
	}

	/**
	 * This is this way to use an ExecutorService some day and also be
	 * able to get separate results.
	 */
	public SingleRowElementaryReversibleMatrix call() throws Exception {
		assert(called == 0);
		called++;
		
		/* Some assertions on the matrix structure ? */
		// Det == 1?
			
		/* Get the P L DR U S0 decomposition first */
		triangularElementaryReversibleMatrix();
			
		/* Compute the other SERM */
		singleRowElementaryReversibleMatrix();
		
		A = null;
		
		return this;
	}
}
