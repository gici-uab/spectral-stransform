package rklt;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Date;
import java.util.LinkedList;
import GiciMatrix.MatrixAlgebra;

class TransformInformation implements Serializable {
	private static final long serialVersionUID = 3L;
	
	final LinkedList<float[][]> operations;
	final LinkedList<int[]> permutations; 
	final LinkedList<Integer> transformTypes; 
	final String creationDate = (new Date()).toString();
	
	public TransformInformation() {
		operations = new LinkedList<float[][]>();
		permutations = new LinkedList<int[]>();
		transformTypes = new LinkedList<Integer>();
	}
	
	// Deep copy constructor
	public TransformInformation(final TransformInformation copyFrom) {
		this();
		
		for (float[][] i : copyFrom.operations) {
			float[][] j = new float[i.length][];
			
			for (int k = 0; k < i.length; k++) {
				j[k] = Arrays.copyOf(i[k], i[k].length); 
			}
			
			this.operations.add(j);
		}
		
		for (int[] i : copyFrom.permutations) {
			int [] j = Arrays.copyOf(i, i.length); 
			
			this.permutations.add(j);
		}
		
		for (int i: copyFrom.transformTypes) {
			this.transformTypes.add(i);
		}
	}
	
	void pushTransformType(int t) {
		transformTypes.add(t);
	}
	
	int popTransformType() {
		return transformTypes.removeLast();
	}
	
	float[][] popDataPacket() {
		return operations.removeLast();
	}
	
	void pushDataPacket(final float[][] a) {
		operations.add(MatrixAlgebra.copy(a));
	}
	
	int[] popPermutation() {
		return permutations.removeLast();
	}
	
	void pushPermutation(final int[] a) {
		permutations.add(Arrays.copyOfRange(a, 0, a.length));
	}
	
	int getSteps() {
		assert (operations.size() == permutations.size());
		assert (operations.size() == transformTypes.size());
		return operations.size();
	}
	
	int getTransformWidth() {
		assert (permutations.size() > 0);
		int[] t = permutations.get(0);
		return t.length;
	}
}