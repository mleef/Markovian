/**
 * Created by marc_leef on 5/20/15.
 * Markov Network construction, variable elimination, and loopy belief propogation.
 */

import java.util.ArrayList;
import java.util.HashMap;


public class MarkovNetwork {
	private int[] cardinalities;
	private int[] variables;
	private ArrayList<Factor> factors;
	private HashMap<Integer, int[]> neighbors;
	private static final double CONVERSION_THRESHOLD = .000001;
	private static final int MAX_ITERS = 100;

    /**
     * Markov Network constructor.
     * @param v Variables in the network.
     * @param c Cardinalities of variables in the network.
     * @param name f List of factors that make up the network.
     * @return New Markov Network.
     */
	public MarkovNetwork(int[] v, int[] c, ArrayList<Factor> f) {
		variables = v;
		cardinalities = c;
		factors = f;
		neighbors = new HashMap<>();
		
		// Calculate neighbors for min neighbors heuristic.
		for(Factor fac : factors) {
			int scope[] = fac.getScope();
			for(int i = 0; i < scope.length; i++) {
				if(neighbors.containsKey(scope[i])) {
					neighbors.put(scope[i], getUnion(neighbors.get(scope[i]), scope));
				}
				else {
					neighbors.put(scope[i], scope);
				}
			}
		}
		
	}
	

    /**
     * Multiplies all factors in the network together. Can be used as a partition function.
     * @return Total value of all factors. Brute forced partition value if no variables eliminated.
     */
	public double multiplyAllFactors() {
		double result = 0.0;
		if(factors.size() == 1) {
			return factors.get(0).sumUp();
		}
		else if(factors.size() == 0) {
			return -1;
		}
		
		Factor f = multiplyFactors(factors.get(0), factors.get(1));
		for(int i = 2; i < factors.size(); i++) {
			f = multiplyFactors(f, factors.get(i));
		}
		
		for(int i = 0; i < f.getSize(); i++){ 
			result += f.getValue(i);
		}
		return result;
	}
	
    /**
     * Constructs full joint distribution of network.
     * @return Factor over all variables in the network.
     */
	public Factor getFullJointDistribution() {
		Factor f = multiplyFactors(factors.get(0), factors.get(1));
		for(int i = 2; i < factors.size(); i++) {
			f = multiplyFactors(f, factors.get(i));
		}
		return f;
	}
	
    /**
     * Markov Network constructor. Adapted from pseudocde found on page 359 of Diane Koller's "Probablistic Graphical Models"
     * @param factor1 First factor to multiply.
     * @param factor2 Second factor to multiply.
     * @return New factor that is product of both inputs.
     */
	public Factor multiplyFactors(Factor factor1, Factor factor2) {
		int j = 0;
		int k = 0;
		int size = 1;
		HashMap<Integer, Integer> assignment = new HashMap<>();
		
		
		// Get union variables
		int[] scope = getUnion(factor1.getScope(), factor2.getScope());
		
		// Calculate size of resulting  factor and initialize assignments.
		for(int i = 0; i < scope.length; i++) {
			assignment.put(scope[i], 0);
			size *= cardinalities[scope[i]];
		}
		
		double[] psi = new double[size];
		Factor result = new Factor(scope, cardinalities);
		result.intializePhi(size);

		for(int i = 0; i < size; i++) {
			// Multiply values using pre-computed indices into both factors and produce new value.
			psi[i] = factor1.getValue(j) * factor2.getValue(k);
			for(int l : scope) {
				assignment.replace(l, assignment.get(l) + 1);
				if(assignment.get(l) == cardinalities[l]) {
					assignment.replace(l, 0);
					// Cardinality reached, move back indices. 
					j = j - (cardinalities[l] - 1) * factor1.getStride(l);
					k = k - (cardinalities[l] - 1) * factor2.getStride(l);
				}
				else {
					// Move indices forward.
					j = j + factor1.getStride(l);
					k = k + factor2.getStride(l);
					break;
				}
			}
			
		}
		
		result.addToPhi(psi);
		return result;
	}
	
    /**
     * Helper function for multiplying list of factors together.
     * @param factorList List of factors to multiply together.
     * @return Product result of list of factors.
     */
	public Factor multiplyFactors(ArrayList<Factor> factorList) {
		Factor result = null;
		if(factorList.size() == 1) {
			return factorList.get(0);
		}
		else {
			result = multiplyFactors(factorList.get(0), factorList.get(1));
			for(int i = 2; i < factorList.size(); i++) {
				result = multiplyFactors(result, factorList.get(i));
			}
		}
		
		return result;
	}
	
    /**
     * Helper function to obtain the union of two factors' variable scopes.
     * @param a1 Scope of first factor.
     * @param a2 Scope of second factor.
     * @return Union of two input scopes.
     */
	public int[] getUnion(int[] a1, int[] a2) {
		int[] tmp = new int[a1.length + a2.length];
		System.arraycopy(a1, 0, tmp, 0, a1.length);
		System.arraycopy(a2, 0, tmp, a1.length, a2.length);
		HashMap<Integer, Boolean> map = new HashMap<>();
		
		for(int i = 0; i < tmp.length; i++) {
			if(map.get(tmp[i]) == null) {
				map.put(tmp[i], true);
			}
			else {
				tmp[i] = -1;
			}
		}
		
		int[] result = new int[map.keySet().size()];
		int j = 0;
		for(int i = 0; i < tmp.length; i++) {
			if(tmp[i] != -1) {
				result[j] = tmp[i];
				j++;
			}
		}
		
		return result;
	}
	
    /**
     * Eliminates list of variables from the network by way of factor marginilization.
     * @param vars List of variables to eliminate.
     */
	public void eliminateVariables(int...vars) {
		for(int var : vars) {
			ArrayList<Factor> targetFactors = new ArrayList<>();
			
			// Aggregate all factors that contain variable for multiplication
			for(int i = 0; i < factors.size(); i++) {
				if(factors.get(i).inScope(var)) {
					targetFactors.add(factors.get(i));
					factors.remove(i);
					i--;
				}
			}
			Factor f;
			
			// Multiply factors together then sum out variable from the resulting factor.
			if(targetFactors.size() > 1) {
				f = multiplyFactors(targetFactors.get(0), targetFactors.get(1));
				for(int i = 2; i < targetFactors.size(); i++) {
					f = multiplyFactors(f, targetFactors.get(i));
				}
				f.sumOut(var);
				factors.add(f);
				continue;
			}
			else if(targetFactors.size() == 1) {
				f = targetFactors.get(0);
				f.sumOut(var);
				factors.add(f);
				continue;
			}
			else {
				continue;
			}
			
		}
		
	}
	
    /**
     * Eliminates variable with the fewest neighbors from the network.
     */
	public void eliminateMinNeighborVariable() {
		// Get next minimum neighbor.
		int var = getMinNeighbor();
		ArrayList<Factor> targetFactors = new ArrayList<>();
		
		// Aggregate all factors that contain variable for multiplication
		for(int i = 0; i < factors.size(); i++) {
			if(factors.get(i).inScope(var)) {
				targetFactors.add(factors.get(i));
				factors.remove(i);
				i--;
			}
		}
		Factor f;
		
		// Multiply factors together then sum out variable from the resulting factor.
		if(targetFactors.size() > 1) {
			f = multiplyFactors(targetFactors.get(0), targetFactors.get(1));
			for(int i = 2; i < targetFactors.size(); i++) {
				f = multiplyFactors(f, targetFactors.get(i));
			}
			f.sumOut(var);
			factors.add(f);
			return;
		}
		else if(targetFactors.size() == 1) {
			f = targetFactors.get(0);
			f.sumOut(var);
			factors.add(f);
			return;
		}
		else {
			return;
		}
	}
	
    /**
     * Prints out all factors in the network.
     */
	public void printFactors() {
		for(Factor f : factors) {
			f.printFactor();
			System.out.println();
		}
	}
	
	
    /**
     * Finds the variable in the network that has the fewest neighbors.
     * @return Variable with the fewest neighbors.
     */
	public int getMinNeighbor() {
		if(neighbors.size() == 0) {
			return -1;
		}
		
		int low = Integer.MAX_VALUE;
		int lowIndex = -1;
		
		for(int key : neighbors.keySet()) {
			if(neighbors.get(key).length < low) {
				low = neighbors.get(key).length;
				lowIndex = key;
			}
		}
		
		neighbors.remove(lowIndex);
		
		return lowIndex;
	}
	
    /**
     * Loopy belief propagation algorithm implementation. Passes messages back and forth
     * between factors and variables until messages converge. By default, the algorithm is
     * capped at 100 iterations.
     */
	public void loopyBP() {
		HashMap<Integer, ArrayList<Factor>> varMessages = new HashMap<>();
		HashMap<Factor, ArrayList<Factor>> factorMessages = new HashMap<>();

		HashMap<Integer, ArrayList<Factor>> previousVarMessages = new HashMap<>();

		// Initialize message lists for both variables and factors.
		int curVariable = 0;
		for(int i = 0; i < variables.length; i++) {
			varMessages.put(variables[i], new ArrayList<Factor>());
		}
		for(int i = 0; i < factors.size(); i++) {
			factorMessages.put(factors.get(i), new ArrayList<Factor>());
		}
		
		int numIterations = 0;
		boolean converged = false;
		
		while(!converged) {
			// Iterate through all variables and generate messages for factors.
			for(int var : variables) {
				for(Factor f : factors) {
					if(f.inScope(var)) {
						
						// Initialize message to with all states equally probable.
						Factor newMessage = new Factor(new int[] {var}, cardinalities);
						newMessage.intializePhi(cardinalities[var]);
						double[] vals = new double[newMessage.getSize()];
						for(int i = 0; i < vals.length; i++) {
							vals[i] = .5;
						}
						newMessage.addToPhi(vals);
						
						// Multiply incoming messages together, excluding those from destination factor.
						for(Factor prevMessage : varMessages.get(var)) {
							if(prevMessage.getFactorSender() != f) {
								newMessage = multiplyFactors(newMessage, prevMessage);
							}
						}
						
						newMessage.setVariableSender(var);
						factorMessages.get(f).add(newMessage);

					}
				}
			}
			
			// Delete previous messages to variables.
			for(int i = 0; i < variables.length; i++) {
				varMessages.replace(variables[i], new ArrayList<Factor>());
			}
						
			// Iterate through all factors and generate new messages for variables.
			for(Factor bigFac : factors) {
				for(int var : bigFac.getScope()) {
					// Start with own factor table
					Factor newMessage = bigFac.copy();

					// If there is one incoming message to the factor, multiply it by the factor itself.
					if (factorMessages.get(bigFac).size() == 1) {
						// Check that message is not from variable we are sending message to.
						if(factorMessages.get(bigFac).get(0).getVariableSender() != var) {
							newMessage = multiplyFactors(factorMessages.get(bigFac).get(0), newMessage);
						}
					}
					// If there are multiple messages, multiply them all together with the original factor.
					else {
						boolean first = true;
						for(Factor prevMessage : factorMessages.get(bigFac)) {
							// Check that message is not from variable we are sending message to.
							if(prevMessage.getVariableSender() != var)  {
								newMessage = multiplyFactors(newMessage, prevMessage);
							}
						}
					}
					
					
					for(int newMessageVar : newMessage.getScope()) {
						// Sum out all variables except of destination node.
						if(newMessageVar != var) {
							newMessage.sumOut(newMessageVar);
						}
					}
					// Add message origin and then send.
					newMessage.setFactorSender(bigFac);
					varMessages.get(var).add(newMessage);
				}
			}
			
			
			
			// Normalize messages
			for(ArrayList<Factor> messages : varMessages.values()) {
				for(Factor f : messages) {
					f.normalize();
				}
			}
			
			// Store previous messages.
			if(numIterations == 0) {
				for(int var : variables) {
					previousVarMessages.put(var, varMessages.get(var));
				}
			}
			// Test for convergence.
			else {
				boolean sigDiff = false;
				for(int var : variables) {
					ArrayList<Factor> vars = varMessages.get(var);
					ArrayList<Factor> prevVars = previousVarMessages.get(var);
					for(int i = 0; i < vars.size(); i++) {
						if(messageDifference(vars.get(i), prevVars.get(i)) > CONVERSION_THRESHOLD) {
							sigDiff = true;
						}
					}
					previousVarMessages.replace(var, varMessages.get(var));
				}
				
				// If no significant difference in messages, algorithm has converged.
				if(!sigDiff) {
					converged = true;
				}
			}
			
			
			// Delete previous messages to factors.
			for(int i = 0; i < factors.size(); i++) {
				factorMessages.replace(factors.get(i), new ArrayList<Factor>());
			}
			
			numIterations++;
			
			// Cap iterations at 100.
			if(numIterations > 100) {
				System.out.println("Exceeded " + MAX_ITERS + " iterations... breaking and calculating marginals.");
				converged = true;
			}
		}
		
		// Compute marginals.
		for(int var : variables) {
			Factor result = multiplyFactors(varMessages.get(var));
			double total = result.sumUp();
			for(int i = 0; i < result.getSize(); i++) {
				System.out.print(result.getValue(i)/total + " ");
			}
			System.out.println();
		}
		
		
	}
	
    /**
     * Calculates the variance between messages for use in loopy BP.
     * @param a First message to analyze.
     * @param b Second message to analyze.
     * @return Averaged differences in values between the two messages.
     */
	public double messageDifference(Factor a, Factor b) {
		double result = 0.0;
		for(int i = 0; i < a.getSize(); i++) {
			result += Math.abs(a.getValue(i) - b.getValue(i));
		}
		
		return result/a.getSize();
		
	}
	
}
