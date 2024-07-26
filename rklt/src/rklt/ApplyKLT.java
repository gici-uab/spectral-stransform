package rklt;

import GiciAnalysis.ImageCovariance;
import GiciMatrix.MatrixAlgebra;
import GiciTransform.KarhunenLoeveTransform;
import GiciTransform.LinearTransform;
import GiciTransform.ZeroMean;

public class ApplyKLT implements ApplyTransform {
	
	final int size;
	float [][][] image;
	
	float[] means;

	float [][] klt;
	float [] eigenvalues;
	
	final int dimension;
	final float subsampling;
	
	public ApplyKLT(final int _size, final int _dimension) {
		this(_size, _dimension, 1);
	}
	
	public ApplyKLT(final int _size, final int _dimension, final float _subsampling) {
		size = _size;
		dimension = _dimension;
		
		means = null;
		
		subsampling = _subsampling;
		
		assert (subsampling > 0 && subsampling <= 1);
	}
	
	/* (non-Javadoc)
	 * @see rklt.ApplyTransform#apply()
	 */
	public void apply() throws Exception {
		// Correct to zero mean data
		ZeroMean zm = new ZeroMean(image, dimension, false);
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

		//System.out.println("apply");
		assert(dimension == 0);
		
		LinearTransform lt = new LinearTransform(klt, dimension);
		lt.transformInPlace(image);
	}
	
	/* (non-Javadoc)
	 * @see rklt.ApplyTransform#remove()
	 */
	public void remove() {
		assert(dimension == 0);
		
		LinearTransform lt = new LinearTransform(MatrixAlgebra.transposeC(klt), dimension);
		lt.transformInPlace(image);

		ZeroMean zm = new ZeroMean(image, dimension, means);
		zm.applyMeansInPlace();
	}
	
	/* (non-Javadoc)
	 * @see rklt.ApplyTransform#getPackedData()
	 */
	public float[][] getPackedData () {
		float[][] r = new float[size+1][size];
		
		for (int i = 0; i < size; i++) {
			System.arraycopy(klt[i], 0, r[i], 0, size);
			r[size][i] = means[i];
		}
	
		return r;
	}
	
	/* (non-Javadoc)
	 * @see rklt.ApplyTransform#setPackedData(float[][])
	 */
	public void setPackedData (final float[][] r) {
		assert(r.length == size + 1 && r[0].length == size);
	
		klt = new float[size][size];
		means = new float[size];
		
		for (int i = 0; i < size; i++) {
			System.arraycopy(r[i], 0, klt[i], 0, size);
			means[i] = r[size][i];
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
	
	/* debug */
	
	/* (non-Javadoc)
	 * @see rklt.ApplyTransform#dumpInternalState()
	 */
	public void dumpInternalState () {
		MatrixAlgebra.printMatrix(klt);
		MatrixAlgebra.printVector(means);
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
