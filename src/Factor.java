/**
 * Created by marc_leef on 5/20/15.
 * Factor construction and calculations.
 */
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;


public class Factor {
	private double[] phi;
	private int[] scope;
	private int curPosition = 0;
	private int size;
	private int[] cardinalities;
	private Factor facMessageOrigin;
	private HashMap<Integer, Integer> strides = new HashMap();
	private int varMessageOrigin;
	
    /**
     * Factor constructor.
     * @param v Variables in the factor.
     * @param cards Cardinalities of variables in the network.
     * @return New Factor.
     */
	public Factor(int[] v, int[] cards) {
		scope = v;
		
		// Reverse variables for stride construction
		for(int i = 0; i < scope.length / 2; i++)
		{
		    int temp = scope[i];
		    scope[i] = scope[scope.length - i - 1];
		    scope[scope.length - i - 1] = temp;
		}
		cardinalities = cards;
	}
	
    /**
     * Used as a copy constructor in loopy bp.
     * @param v Variables in the factor.
     * @param cards Cardinalities of variables in the network.
     * @param str Strides of variables in the factor.
     * @param p Values in the factor.
     * @param s Size of the factor.
     * @return New Markov Network.
     */
	public Factor(int[] v, int[] cards, HashMap<Integer, Integer> str, double[] p, int s) {
		scope = v;
		cardinalities = cards;
		strides = str;
		phi = p;
		size = s;
	}
	
    /**
     * Initializes the factor so it can store values.
     * @param length Length of the value array.
     */
	public void intializePhi(int length) {
		phi = new double[length];
		size = length;
		
		
		// Calculate stride lengths of reversed list.
		for(int i = 0; i < scope.length; i++) {
			if(i == 0) {
				strides.put(scope[i], 1);
			}
			else{
				strides.put(scope[i], strides.get(scope[i-1]) * cardinalities[scope[i-1]]);
			}
		}
	}
	
    /**
     * Adds values to factor.
     * @param vals Values to add to factor (strings).
     */
	public void addToPhi(String[] vals) {
		for(int i = 0; i < vals.length; i++) {
			phi[curPosition] = Double.parseDouble(vals[i]);
			curPosition++;
		}

	}
	
    /**
     * Adds values to factor.
     * @param vals Values to add to factor (doubles).
     */
	public void addToPhi(double[] vals) {
		for(int i = 0; i < vals.length; i++) {
			phi[curPosition] = vals[i];
			curPosition++;
		}
	}
	
    /**
     * Gets the stride of input variable, 0 if variable is not in scope of factor.
     * @param var Variable to find stride of.
     * @return Stride of variable.
     */
	public int getStride(int var) {
		if(strides.get(var) != null) {
			return strides.get(var);
		}
		return 0;
	}
	
    /**
     * Gets value from factor.
     * @param index Index of wanted value.
     * @return Value at specified index within the factor.
     */
	public double getValue(int index) {
		return phi[index];
	}
	
    /**
     * Gets size of factor.
     * @return Size of factor.
     */
	public int getSize() {
		return size;
	}
	
    /**
     * Returns all variables in scope of factor.
     * @return Scope of factor.
     */
	public int[] getScope() {
		return scope;
	}
	
    /**
     * Helper function to print factor values.
     */
	public void printFactor() {
		System.out.print("Variables: ");
		for(int i = 0; i < scope.length; i++) {
			System.out.print(scope[i] + " ");
		}
		System.out.println();
		System.out.print("Values: ");
		for(int i = 0; i < phi.length; i++) {
			System.out.print(phi[i] + " ");
		}
		System.out.println();
	}
	
    /**
     * Factor marginalization algorithm.
     * @param var Variable to sum out of factor.
     */
	public void sumOut(int var) {
		boolean inScope = false;
		for(int i = 0; i < scope.length; i++) {
			if(var == scope[i]) {
				inScope = true;
			}
		}
		if(!inScope) {
			return;
		}
		

		int tempSize = this.getSize()/cardinalities[var];
		int[] tempScope = new int[scope.length - 1];
		double[] tempPhi = new double[tempSize];
		HashMap<Integer, Integer> tempStrides = new HashMap();
		ArrayList<double[]> subLists = new ArrayList<>();
			
		double accum = 0;
		int curIndex = 0;

		
		HashMap<Integer, Integer> assignment = new HashMap<>();

		// Remove variable from scope.
		int j = 0;
		for(int i = 0; i < scope.length; i++) {
			if(scope[i] != var) {
				tempScope[j] = scope[i];
				assignment.put(tempScope[j], 0);
				j++;
			}
		}

		
		
		
		// Calculate new strides.
		for(int i = 0; i < tempScope.length; i++) {
			if(i == 0) {
				tempStrides.put(tempScope[i], 1);
			}
			else{
				tempStrides.put(tempScope[i], tempStrides.get(tempScope[i-1]) * cardinalities[tempScope[i-1]]);
			}
		}
		
		tempStrides.put(var, 0);
		assignment.put(var, 0);
		
		// Generate assignments and sum all values 
		for(int i = 0; i < size; i++) {			
			tempPhi[curIndex] += phi[i];
			for(int l : scope) {
				assignment.replace(l, assignment.get(l) + 1);
				if(assignment.get(l) == cardinalities[l]) {
					assignment.replace(l, 0);
					curIndex = curIndex - (cardinalities[l] - 1) * tempStrides.get(l);
				}
				else {
					curIndex = curIndex + tempStrides.get(l);
					break;
				}
			}
			
		}
		
		// Update lists with new values.
		phi = tempPhi;
		scope = tempScope;
		size = tempSize;
		strides = tempStrides;
	}
	
    /**
     * Summs up all values in factor.
     * @return Sum of values in factor.
     */
	public double sumUp() {
		double result = 0;
		for(int i = 0; i < phi.length; i++) {
			result += phi[i];
		}
		return result;
	}
	
    /**
     * Checks if input variable is in the scope of the factor.
     * @return True if variable is in scope, false otherwise.
     */
	public boolean inScope(int var) {
		for(int i = 0; i < scope.length; i++) {
			if(scope[i] == var) {
				return true;
			}
		}
		return false;
	}
	
    /**
     * Creates a copy of factor.
     * @param Copied factor.
     */
	public Factor copy() {
		return new Factor(scope, cardinalities, strides, phi, size);
	}
	
    /**
     * Setter function for message origin for use in Loopy BP.
     * @param f Factor that message is coming from.
     */
	public void setFactorSender(Factor f) {
		facMessageOrigin = f;
	}
	
    /**
     * Getter function for message origin for use in Loopy BP.
     * @return Factor message origin.
     */
	public Factor getFactorSender() {
		return facMessageOrigin;
	}
	
    /**
     * Setter function for message origin for use in Loopy BP.
     * @param v Variable that message is coming from.
     */
	public void setVariableSender(int v) {
		varMessageOrigin = v;
	}
	
    /**
     * Getter function for message origin for use in Loopy BP.
     * @return Variable message origin.
     */
	public int getVariableSender() {
		return varMessageOrigin;
	}
	
    /**
     * Normalizes values in factor for use in Loopy BP.
     */
	public void normalize() {
		double max = Double.MIN_VALUE;
		for(int i = 0; i < phi.length; i++) {
			if(phi[i] > max) {
				max = phi[i];
			}
		}
		for(int i = 0; i < phi.length; i++) {
			phi[i] = phi[i]/max;
		}
	}
	
	
	

}
