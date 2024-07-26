package rklt;

import GiciAnalysis.ImageCovariance;
import GiciMatrix.MatrixAlgebra;
import GiciTransform.KarhunenLoeveTransform;
import GiciTransform.TriangularElementaryReversibleMatrix;
import GiciTransform.TriangularElementaryReversibleTransform;
import GiciTransform.ZeroMean;

public class ApplyRKLT implements ApplyTransform {
	
	final int size;
	float [][][] image;
	
	float[][] termL;
	float[][] termU;
	float[][] termS;
	float[] means;
	int[] permutation;

	float [][] klt;
	float [] eigenvalues;
	
	final int dimension;
	final float subsampling;
	
	public ApplyRKLT(final int _size, final int _dimension) {
		this(_size, _dimension, 1);
	}
	
	public ApplyRKLT(final int _size, final int _dimension, final float _subsampling) {
		size = _size;
		dimension = _dimension;
		
		termL = null;
		termU = null;
		termS = null;
		means = null;
		permutation = null;
		
		subsampling = _subsampling;
		
		assert (subsampling > 0 && subsampling <= 1);
	}

//	private void zeroRunLength (float[] a, float zero) {
//		int zc = 0;
//		float max = Float.NEGATIVE_INFINITY;
//		
//		System.out.println("Testing: ");
//		for (int i = 0; i < a.length; i++) {
//			System.out.print(a[i] + ", ");
//		}
//		System.out.println("");
//		
//		for (int i = 0; i < a.length; i++) {
//			if (Math.abs(a[i]) < zero) {
//				zc++;
//				max = Math.max(max, Math.abs(a[i]));
//			} else {
//				if (zc != 0) {
//					System.out.println("zero run of length: " + zc + " (max: " + max + ")");
//				}
//				zc = 0;
//				max = Float.NEGATIVE_INFINITY;
//			}
//		}
//		
//		if (zc != 0) {
//			System.out.println("zero run of length: " + zc + " (max: " + max + ")");
//		}
//	}
	
	/* (non-Javadoc)
	 * @see rklt.ApplyTransform#apply()
	 */
	public void apply() throws Exception {
		// Correct to zero mean data
		ZeroMean zm = new ZeroMean(image, dimension, true);
		zm.applyMeansInPlace();
		means = zm.getMeans();
		zm = null;
		
		//System.out.println("covar");
		//Compute the covariance matrix
		final float [][] covariance;
		
		if (subsampling == 1) 
			covariance = ImageCovariance.generateCovarianceMatrix(image, dimension);
		else
			covariance = ImageCovariance.generateCovarianceMatrix(image, dimension, subsampling, subsampling);

		//System.out.println("klt");
		// KLT Transform
		float[][][] kltResult = KarhunenLoeveTransform.generateTransform(covariance); 
		klt = kltResult[0];
		klt = MatrixAlgebra.transposeC(klt);
		eigenvalues = MatrixAlgebra.matrixDiagonal(kltResult[1]);
		
		TriangularElementaryReversibleMatrix term = (new TriangularElementaryReversibleMatrix(klt)).call();

		//System.out.println("apply");
		assert(dimension == 0);
		image = (new TriangularElementaryReversibleTransform(term, image, true)).call();
		//System.out.println("done apply");
		
		permutation = term.getPermutation();
		termL = term.getElementaryMatrixL();
		termU = term.getElementaryMatrixU();
		termS = term.getElementaryMatrixS();
		
//		/* count (+/-)zero runs */
//		for (int i = 0; i < termL.length; i++) {
//			zeroRunLength(termL[i], 0.00001f);
//		}
//		for (int i = 0; i < termL.length; i++) {
//			zeroRunLength(termU[i], 0.00001f);
//		}
//		for (int i = 0; i < termL.length; i++) {
//			zeroRunLength(termS[i], 0.00001f);
//		}
		
	}
	
	/* (non-Javadoc)
	 * @see rklt.ApplyTransform#remove()
	 */
	public void remove() {
		assert(dimension == 0);
		image = (new TriangularElementaryReversibleTransform(permutation, termL, termU, termS, image, false)).call();

		ZeroMean zm = new ZeroMean(image, dimension, means);
		zm.applyMeansInPlace();
	}
	
	/* (non-Javadoc)
	 * @see rklt.ApplyTransform#getPackedData()
	 */
	public float[][] getPackedData () {
		float[][] r = new float[3*size+2][size];
		
		for (int i = 0; i < size; i++) {
			System.arraycopy(termL[i], 0, r[i + size * 0], 0, size);
			System.arraycopy(termU[i], 0, r[i + size * 1], 0, size);
			System.arraycopy(termS[i], 0, r[i + size * 2], 0, size);

			r[size * 3][i] = means[i];
			r[size * 3 + 1][i] = (float)permutation[i];
		}
	
		return r;
	}
	
	/* (non-Javadoc)
	 * @see rklt.ApplyTransform#setPackedData(float[][])
	 */
	public void setPackedData (final float[][] r) {
		assert(r.length == 3*size+2 && r[0].length == size);
		
		termL = new float[size][size];
		termU = new float[size][size];
		termS = new float[size][size];
		means = new float[size];
		permutation = new int[size];
		
		for (int i = 0; i < size; i++) {
			System.arraycopy(r[i + size * 0], 0, termL[i], 0, size);
			System.arraycopy(r[i + size * 1], 0, termU[i], 0, size);
			System.arraycopy(r[i + size * 2], 0, termS[i], 0, size);

			means[i] = r[size * 3][i];
			permutation[i] = (int)r[size * 3 + 1][i];
		}
	}
	
	/* (non-Javadoc)
	 * @see rklt.ApplyTransform#getImage()
	 */
	public float[][][] getImage() {
		return image;
	}

	/* (non-Javadoc)
	 * @see rklt.ApplyTransform#setImage(float[][][])
	 */
	public void setImage(float[][][] image) {
		this.image = image;
	}
	
	/* some helper functions to help deal with packed data */
//	static public float[][] head (final float[][] a, final int size) {
//		float r[][] = new float[size * 3 + 2][];
//		
//		System.arraycopy(a, 0, r, 0, r.length);
//		
//		return r;
//	}
//	
//	static public float[][] tail (final float[][] a, final int size) {
//		float r[][] = new float[a.length - (size * 3 + 2)][];
//		
//		System.arraycopy(a, size * 3 + 2, r, 0, r.length);
//
//		return r;
//	}
//	
//	static public float[][] cat (final float[][] a, final float[][] b) {
//		float r[][] = new float[a.length + b.length][];
//		
//		System.arraycopy(a, 0, r, 0, a.length);
//		System.arraycopy(b, 0, r, a.length, b.length);
//		
//		return r;
//	}
	
	/* debug */
	
	/* (non-Javadoc)
	 * @see rklt.ApplyTransform#dumpInternalState()
	 */
	public void dumpInternalState () {
		MatrixAlgebra.printMatrix(termL);
		MatrixAlgebra.printMatrix(termU);
		MatrixAlgebra.printMatrix(termS);
		MatrixAlgebra.printVector(means);
		
		float[] r = new float[permutation.length];
		for (int i = 0; i < permutation.length; i++) {
			r[i] = (float) permutation[i];
		}
		
		MatrixAlgebra.printVector(r);
	}

	/* (non-Javadoc)
	 * @see rklt.ApplyTransform#getKltMatrix()
	 */
	public float[][] getKltMatrix() {
		return klt;
	}
	
	/* (non-Javadoc)
	 * @see rklt.ApplyTransform#getEigenvalues()
	 */
	public float[] getEigenvalues() {
		return eigenvalues;
	}
}
