package GiciTransform;

import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.concurrent.Callable;

import GiciException.WarningException;
import GiciMatrix.MatrixAlgebra;

public class TriangularElementaryReversibleMatrix implements Callable {
	int called = 0;

	final int N;

	float[][] A;
	float[][] S;
	float[][] P;
	float[][] L;

	private class State {		
		public float[][] Ps;
		public float[][] Ls;
		public float[][] Ss;
		public float[][] A;
		public float result;
		
		public State() {
			super();
		}
		
		public State(State copyme) {
			super();
			this.Ps = MatrixAlgebra.copy(copyme.Ps);
			this.Ls = MatrixAlgebra.copy(copyme.Ls);
			this.Ss = MatrixAlgebra.copy(copyme.Ss);
			this.A = MatrixAlgebra.copy(copyme.A);
			this.result = copyme.result;
		}
	};

	class BestMove implements Comparator<BestMove> {
		public int permuted_row = -1;
		public int permuted_col = -1;
		public float min = Float.POSITIVE_INFINITY;
		
		public BestMove (){
			
		}
		
		public BestMove(int permuted_row, int permuted_col, float min) {
			super();
			this.permuted_row = permuted_row;
			this.permuted_col = permuted_col;
			this.min = min;
		}

		public int compare(BestMove o1, BestMove o2) {
			if (o1.min < o2.min)
				return -1;
			else if (o1.min > o2.min)
				return 1;
			else 
				return 0;
		}
		
		
	};
	
	private State nextState(final float[][] A, final int i, final int choice) throws WarningException {
		State r = new State();
		
		float[][] Ap = MatrixAlgebra.copy(A);
		float[][] Ps = MatrixAlgebra.identityC(N);
		float[][] Ls = MatrixAlgebra.identityC(N);
		float[][] Ss = MatrixAlgebra.identityC(N);

		/* Chose the "best" permutation for quasi-complete pivoting */

		int permuted_row = -1;
		int permuted_col = -1;

		float min = Float.POSITIVE_INFINITY;

		PriorityQueue<BestMove> pq = new PriorityQueue<BestMove>(1, new BestMove()); 
		
		for (int j = i; j < N; j++) {
			for (int k = i + 1; k < N; k++) {
				if (Math.abs(Ap[j][k]) < Float.MIN_VALUE)
					continue;

				float test = Math.abs((Ap[j][i] - 1) / Ap[j][k]); 

				if (test < min) {
					min = test;
					permuted_row = j;
					permuted_col = k;
				}
				
				pq.add(new BestMove(j, k, test));
			}
		}

//		if (permuted_row < 0 || permuted_col < 0) {
		if (pq.size() == 0) {
			MatrixAlgebra.printMatrix(Ap);
			throw new WarningException("Singular matrix (" + i + ")");
		}

		int realChoice = choice % pq.size();
		
		do {
			BestMove bm = pq.poll();
			
			assert (choice != 0 || bm.min == min);
			
			min = bm.min;
			permuted_row = bm.permuted_row;
			permuted_col = bm.permuted_col;
		} while (realChoice-- > 0);
		
		r.result = min;
		
		/* Apply the permutation (swap rows i and permuted_row)*/
		{
			float[] t;

			//MatrixAlgebra.printMatrix(A);

			t = Ap[i];
			Ap[i] = Ap[permuted_row];
			Ap[permuted_row] = t;

			//MatrixAlgebra.printMatrix(A);

			// In this case rows or columns doesn't matter
			t = Ps[i];
			Ps[i] = Ps[permuted_row];
			Ps[permuted_row] = t;

			//MatrixAlgebra.printMatrix(Ps[i]);
		}

		/* Extract the Ss[i] part */
		float s = (Ap[i][i] - 1) / Ap[i][permuted_col];
		// Change sign, and latter use substraction instead of addition.
		Ss[permuted_col][i]= s; 

		/* Apply Ss to A */
		for (int j = 0; j < N; j++) {
			Ap[j][i] = Ap[j][i] - s * Ap[j][permuted_col];
		}
//		A = MatrixAlgebra.multiplicationCC(A, Ss);

		/* Invert Ss */
		//Ss[permuted_col][i]= s;

		/* Gaussian elimination */
		for (int j = i + 1; j < N; j++) {
			Ls[j][i] = - Ap[j][i];

			Ap[j][i] = 0;

			for (int k = i + 1; k < N; k++) {
				Ap[j][k] = Ap[j][k] + Ls[j][i] * Ap[i][k];
			}
		}

		r.A = Ap;
		r.Ls = Ls;
		r.Ss = Ss;
		r.Ps = Ps;
		
		return r;
	}
	
	private void triangularElementaryReversibleMatrix () throws WarningException {
		
		//MatrixAlgebra.printMatrix(A);
		
		float[][] Ps = new float[N][N];
		float[][] Ls = new float[N][N];
		float[][] Ss = new float[N][N];
		
		P = MatrixAlgebra.identityC(N);
		S = MatrixAlgebra.identityC(N);
		float[][] Lm1 = MatrixAlgebra.identityC(N);
		
		for (int i = 0; i < N-1; i++) {
			
			int lookAhead = 1; //3
			final int lookAheadWidth = 1; //3
			float bestResult = Float.POSITIVE_INFINITY;
			State bestState = null;
			int[] lookAheadChoice = new int[lookAhead];
			int j;
			
			State[] nextState = new State[lookAhead + 1];
			for (int k = 0; k < lookAhead + 1; k++) {
				nextState[k] = new State();
			}
			nextState[0].A = A;
			nextState[0].result = 1;
			
			/* we won't look past the end */
			lookAhead = Math.min(lookAhead, N - 1 - i);
			
			for (j = 0; j < lookAhead; j++) {
				//System.out.println("j " + j + " choice " + lookAheadChoice[j]);
				
				nextState[j + 1] = nextState(nextState[j].A, i + j, lookAheadChoice[j]);
				
				/* if we are at the end test for quality */
				if (j == lookAhead - 1) {
					/* is this result better? */
					
					float newResult = 1;
					for (int k = 0; k < lookAhead + 1; k++) {
						newResult *= nextState[k].result;
					}
					
					if (bestResult > newResult) {
						bestResult = nextState[j+1].result;
						bestState = new State(nextState[1]);
					}
					
					/* increase and backtrack */
					while (j >= 0) {
						lookAheadChoice[j] = (lookAheadChoice[j] + 1) % lookAheadWidth;
						
						if(lookAheadChoice[j--] != 0)
							break;
					}
					
					if (j < 0)
						break;
				}
			}
			
			bestState = nextState(nextState[0].A, i, 0);
			
			Ps = bestState.Ps;
			Ls = bestState.Ls;
			Ss = bestState.Ss;
			A  = bestState.A;
		
			/* Online step merging.
			 * These three variables are results,
			 * and are not used in the following
			 * loop iterations. 
			 */
			
			/* S = Id . Ss[N-2] ... Ss[i] ... Ss[0] */
			S = MatrixAlgebra.multiplicationCC(Ss, S);
			
			/* P = (Id . Ps[N-2] ... Ps[i] ... Ps[0]) ^ T */
			P = MatrixAlgebra.multiplicationCC(Ps, P);
			
			/* Lm1 = Ls[N-2]P[N-2] ... Ls[0]P[0] . P */
			Lm1 = MatrixAlgebra.multiplicationCC(MatrixAlgebra.multiplicationCC(Ls, Ps), Lm1);
		}
		
		// And we have P
		P = MatrixAlgebra.transposeC(P);
		
		// DR * U = A, but there is no need to split it.		
		Lm1 = MatrixAlgebra.multiplicationCC(Lm1, P);
		
		L = Lm1;
		
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
		
		// A = U
		// And we already have P L U S
	}
	
	/**
	 * 
	 * @param a read by reference!
	 */
	public TriangularElementaryReversibleMatrix (float [][] a) {
		// A must be a square matrix
		assert(a.length == a[0].length);
		
		A = MatrixAlgebra.copy(a);
		N = A.length;
	}

	public float[][] getElementaryMatrixL () {
		assert(called > 0);
		
		return L;
	}
	
	public float[][] getElementaryMatrixU () {
		assert(called > 0);
		
		return A;
	}
	
	public float[][] getElementaryMatrixS () {
		assert(called > 0);
		
		return S;
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
	public TriangularElementaryReversibleMatrix call() throws Exception {
		assert(called == 0);
		called++;
		
		/* Some assertions on the matrix structure ? */
		// Det == 1?
			
		/* Get the P L DU S decomposition */
		if (N > 1) {
			triangularElementaryReversibleMatrix();
		} else {
			A = MatrixAlgebra.identityC(1);
			S = MatrixAlgebra.identityC(1);
			P = MatrixAlgebra.identityC(1);
			L = MatrixAlgebra.identityC(1);
		}
		
		return this;
	}
}
