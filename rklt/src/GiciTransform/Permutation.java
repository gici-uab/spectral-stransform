package GiciTransform;

import java.util.Arrays;

import GiciMatrix.MatrixAlgebra;

/**
 * A permutation defined as: "2, 0, 1", means that if applied, the resulting vector will have first the 3rd (2) component, and so on. 
 */
public class Permutation {
	int [] permutation;

	/**
	 * Creates an identity permutation.
	 * @param size
	 */
	public Permutation(final int size) {
		permutation = new int[size];
		this.setToIdentity();
	}
	
	/**
	 * Initializes this object with an arbritrary permutation.
	 * @param p
	 */
	public Permutation(final int[] p) {
		permutation = new int[p.length];
		this.setPermutation(p);
	}
	
	public void setToIdentity () {
		for (int i = 0; i < permutation.length; i++) {
			permutation[i] = i;
		}
	}
	
	public boolean isIdentity () {
		for (int i = 0; i < permutation.length; i++) {
			if (permutation[i] != i)
				return false;
		}
		
		return true;
	}
	
	/**
	 * Appends a permutation to the current one that does this.
	 * @param from
	 * @param to
	 * @param size
	 */
	public void permutationBringTo(final int from, final int to, final int size) {
		assert (permutation.length >= to + size);
		assert (permutation.length >= from + size);
		assert (from >= 0 && to >= 0 && size >= 0);
	
		if (to == from)
			return;
		
		// Calculate the zone to be vacated and its destination
		final int vacateFrom, vacateTo, vacateSize;

		if (to + size <= from || from + size <= to) {
			// No overlap (just swap the zones)
			vacateFrom = to;
			vacateTo = from;
			vacateSize = size;
		} else if (to < from) {
			vacateFrom = to;
			vacateTo = to + size;
			vacateSize = from - to;
		} else { //if (to > from) {
			vacateFrom = from + size;
			vacateTo = from;
			vacateSize = to - from;
		}
	
		assert(vacateFrom + vacateSize <= permutation.length);
		assert(vacateSize >= 0);
		
		int[] vacate = Arrays.copyOfRange(permutation, vacateFrom, vacateFrom + vacateSize);
		
		System.arraycopy(permutation, from, permutation, to, size);
		System.arraycopy(vacate, 0, permutation, vacateTo, vacateSize);
	}
	
	public void appendPermutation(Permutation p) {
		Permutation t = new Permutation(p.getPermutation());
		permutation = t.applyToVector(permutation);
	}
	
	public void setPermutation(int[] permutation) {
		assert(this.permutation.length == permutation.length);
		this.permutation = permutation;
		
		/* Check permutation for correctness */
		int[] permTest = new int[permutation.length];
		System.arraycopy(permutation, 0, permTest, 0, permutation.length);
		Arrays.sort(permTest);
		for (int i = 0; i < permTest.length; i++) {
			assert(permTest[i] == i);
		}
	}
	
	public int[] getPermutation() {
		return permutation;
	}
	
	public float[][] getPermutationAsMatrix() {
		return applyToMatrixRows(MatrixAlgebra.identityC(permutation.length));
	}
	
	/**
	 * This permutation is modified in a way that applying it before and
	 * after modification would result in no permutation performed.
	 */
	public void reversePermutation() {
		int[] newPermutation = new int[permutation.length];
		
		for (int i = 0; i < permutation.length; i++) {
			newPermutation[permutation[i]] = i;
		}
		
		permutation = newPermutation;
	}

	public void permuteInPlace(final float[][][] image, final int dimension){
		int zSize;
		int ySize;
		int xSize;
		
		final float newImage[][][];
			
		zSize = image.length;
		ySize = image[0].length;
		xSize = image[0][0].length;
		
		assert (dimension >= 0 && dimension <= 2);
		assert((dimension == 0 && image.length == permutation.length)
				|| (dimension == 1 && image[0].length == permutation.length)
				|| (dimension == 2 && image[0][0].length == permutation.length));		
		
		switch(dimension) {
		case 0:
			assert(permutation.length == zSize);
			
			newImage = new float[zSize][][];
			
			// Careful with trying to do it in place here
			// as overlaps may yield unexpected results
			for (int z = 0; z < zSize; z++) {
				newImage[z] = image[permutation[z]];
			}
			
			break;
		case 1:
			assert(permutation.length == ySize);
			
			newImage = new float[zSize][ySize][];
			
			for (int z = 0; z < zSize; z++) {
				for (int y = 0; y < ySize; y++) {
					newImage[z][y] = image[z][permutation[y]];
				}
			}
			
			break;
		case 2:
			assert(permutation.length == xSize);
			
			newImage = new float[zSize][ySize][xSize];
			
			for (int z = 0; z < zSize; z++) {
				for (int y = 0; y < ySize; y++) {
					for (int x = 0; y < ySize; y++) {
						newImage[z][y][x] = image[z][y][permutation[x]];
					}
				}
			}
			
			break;
			
			default:
				throw new IndexOutOfBoundsException();
		}
		
		//image = newImage;
		System.arraycopy(newImage, 0, image, 0, newImage.length);
	}
	
	public float[][] applyToMatrixRows(final float[][] a) {
		float[][][] b = {MatrixAlgebra.copy(a)};
		
		assert(permutation.length == a.length);
		
		Permutation p = new Permutation(permutation);
		p.permuteInPlace(b, 1);
		
		return b[0];
	}
	
	public float[][] applyToMatrixCols(final float[][] a) {
		float[][][] b = {MatrixAlgebra.copy(a)};
		
		assert(permutation.length == a[0].length);
		
		Permutation p = new Permutation(permutation);
		p.permuteInPlace(b, 2);
		
		return b[0];
	}
	
	public int[] applyToVector(int[] v) {
		int[] r = new int[v.length];
		
		for (int z = 0; z < r.length; z++) {
			r[z] = v[permutation[z]];
		}
		
		return r;
	}
	
	/* debug */
	public void scream() {
		System.out.print("Permutation " + this + " ");
		MatrixAlgebra.printVector(permutation);
		
		(new Exception()).printStackTrace(System.out);
	}
}
