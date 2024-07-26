package GiciTransform;

import java.util.concurrent.Callable;

public strictfp class TriangularElementaryReversibleTransform implements Callable {

	final int[] permutation;
	final float[][] termL;
	final float[][] termU;
	final float[][] termS;
	
	long[][][] image;
	boolean forward;
	
	public TriangularElementaryReversibleTransform(final int[] permutation,
			final float[][] termL, final float[][] termU, final float[][] termS,
			float[][][] image, final boolean forward) {

		this.permutation = permutation; 
		this.termL = termL;
		this.termU = termU;
		this.termS = termS;
		this.forward = forward;
		
		// Asserts 
		assert (termL.length == image.length
				&& termU.length == image.length
				&& termS.length == image.length
				&& termL.length == termL[0].length 
				&& termU.length == termU[0].length
				&& termS.length == termS[0].length);
		
		assert (permutation.length == image.length);
		
		// copy: no in place calculation
		this.image = new long[image.length][image[0].length][image[0][0].length];
		
		for (int z = 0; z < image.length; z++) {
			for (int y = 0; y < image[0].length; y++) {
				for (int x = 0; x < image[0][0].length; x++) {
					this.image[z][y][x] = (int)image[z][y][x];
				}
			}
		}
		// end of copy
		
	}

	public TriangularElementaryReversibleTransform(final TriangularElementaryReversibleMatrix term,
			float[][][] image, final boolean forward) {

		permutation = term.getPermutation();
		termL = term.getElementaryMatrixL();
		termU = term.getElementaryMatrixU();
		termS = term.getElementaryMatrixS();
		this.forward = forward;
		
		// Asserts 
		assert (termL.length == image.length
				&& termU.length == image.length
				&& termS.length == image.length
				&& termL.length == termL[0].length 
				&& termU.length == termU[0].length
				&& termS.length == termS[0].length);
		
		assert (permutation.length == image.length);
		
		// copy: no in place calculation
		this.image = new long[image.length][image[0].length][image[0][0].length];
		
		for (int z = 0; z < image.length; z++) {
			for (int y = 0; y < image[0].length; y++) {
				for (int x = 0; x < image[0][0].length; x++) {
					this.image[z][y][x] = (int)image[z][y][x];
				}
			}
		}
		// end of copy
		
	}

	private long[] applyLowerTerm (float[][] term, long[] data) {
		int N = data.length;		
		assert (term.length == N && term[0].length == N);
		
		// Bottom-up application
		for (int i = N-1; i >= 0; i--) {
			double acc = 0;
			
			for (int j = 0; j < i; j++) {
				acc += (double) term[i][j] * (double) data[j];
			}
			
			data[i] += Math.rint(acc);
		}
		
		return data;
	}
	
	private long[] applyUpperTerm (float[][] term, long[] data) {
		int N = data.length;		
		assert (term.length == N && term[0].length == N);
		
		// Top-down application
		for (int i = 0; i < N; i++) {
			double acc = 0;
			
			for (int j = i + 1; j < N; j++) {
				acc += (double) term[i][j] * (double) data[j];
			}
			
			data[i] += Math.rint(acc);
		}
		
		// Careful with the possible sign change in the last component!
		if (term[N-1][N-1] < 0) {
			data[N-1] = -data[N-1];
		}
		
		return data;	
	}
	
	private long[] applyPermutation (int[] permutation, long[] data) {
		// Apply the permutation
		long[] t = new long[data.length];
		
		for (int z = 0; z < data.length; z++) {
			t[z] = data[permutation[z]];
		}
		
		return t;
	}
	
	private long[] removeLowerTerm (float[][] term, long[] data) {
		int N = data.length;		
		assert (term.length == N && term[0].length == N);
				
		for (int i = 0; i < N; i++) {
			double acc = 0;
			
			for (int j = 0; j < i; j++) {
				acc += (double) term[i][j] * (double) data[j];
			}

			data[i] -= Math.rint(acc);
		}
		
		return data;
	}
	
	private long[] removeUpperTerm (float[][] term, long[] data) {
		int N = data.length;		
		assert (term.length == N && term[0].length == N);
		
		// Careful with the possible sign change in the last component!
		if (term[N-1][N-1] < 0) {
			data[N-1] = -data[N-1];
		}
		
		for (int i = N - 1; i >= 0; i--) {
			double acc = 0;
			
			for (int j = i + 1; j < N; j++) {
				acc += (double) term[i][j] * (double) data[j];
			}
			
			data[i] -= Math.rint(acc);
		}
		
		return data;	
	}
	
	private int[] inversePermutation (int[] permutation) {
		int[] ipermutation = new int[permutation.length];
		
		for (int i = 0; i < permutation.length; i++) {
			ipermutation[permutation[i]] = i;
		}
		
		return ipermutation;
	}
	
	public float[][][] call() {

		int N = image.length;
		
		//MatrixAlgebra.printMatrix(serms);
		
		if (forward) {
			for (int y = 0; y < image[0].length; y++) {
				for (int x = 0; x < image[0][0].length; x++) {

					long[] tl = new long[N];
					
					for (int z = 0; z < N; z++) {
						tl[z] = image[z][y][x];
					}

					tl = applyLowerTerm(termS, tl);
					tl = applyUpperTerm(termU, tl);
					tl = applyLowerTerm(termL, tl);
					tl = applyPermutation(permutation, tl);
					
					for (int z = 0; z < N; z++) {
						image[z][y][x] = tl[z];
					}
				}
			}
		} else {
			// Reverse

			int[] ipermutation = inversePermutation(permutation);
			
			for (int y = 0; y < image[0].length; y++) {
				for (int x = 0; x < image[0][0].length; x++) {

					long[] tl = new long[N];
					
					for (int z = 0; z < N; z++) {
						tl[z] = image[z][y][x];
					}

					tl = applyPermutation(ipermutation, tl);
					tl = removeLowerTerm(termL, tl);
					tl = removeUpperTerm(termU, tl);
					tl = removeLowerTerm(termS, tl);
					
					for (int z = 0; z < N; z++) {
						image[z][y][x] = tl[z];
					}
				}
			}
		}
		
		// copy: no in place calculation
		float[][][] t = new float[image.length][image[0].length][image[0][0].length];
		
		for (int z = 0; z < image.length; z++) {
			for (int y = 0; y < image[0].length; y++) {
				for (int x = 0; x < image[0][0].length; x++) {
					
					// Anything outside this boundary is not lossless
					assert(image[z][y][x] < (1 << 23) &&
							image[z][y][x] > -(1 << 23));
					
					t[z][y][x] = image[z][y][x];
				}
			}
		}
		
		return t;
	}


}
