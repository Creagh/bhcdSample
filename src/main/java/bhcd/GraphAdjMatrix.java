/**
 * 
 */
package bhcd;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.Arrays;

import com.google.common.primitives.Ints;
import com.rits.cloning.Immutable;

import blang.distributions.Generators;
import blang.inits.ConstructorArg;
import blang.inits.DesignatedConstructor;
import briefj.BriefIO;
import briefj.collections.UnorderedPair;


/**
 * A simple, undirected, unweighted graph, represented as an adjacency matrix.
 * 
 * @author Creagh Briercliffe (creagh.briercliffe@gmail.com)
 *
 */
@Immutable
public class GraphAdjMatrix {
	
	public int[] vertices; // an array containing all vertices
  
	protected int e; // number of edges
	
	public int[][] adj; // the adjacency matrix
	
	private static int[] vertices(int size) {
	  int [] vertices = new int[size];
	  for (int i = 0; i< size ;i++) {
		  vertices[i] = i+1; // follow weird 1-index convention for vertices
	  }
	  return vertices;
	}
	
	private GraphAdjMatrix(int[] vertices, int e, int[][] adj, boolean modifiable) {
		if (vertices.length != adj.length) throw new RuntimeException();
		this.vertices = vertices;
		this.e = e;
		this.adj = adj;
	}
	
    /**
     * Construct an empty graph from vertices
     * @param int[] vertices
     */
    public GraphAdjMatrix(int[] vertices) {
        this.vertices = vertices;
        this.e = 0;
        int v = vertices.length;
        this.adj = new int[v][v];
    }
    
    /**
     * Get the number of vertices and edges
     * @return
     */
    public int getNumVertices() { return vertices.length; }
    public int getNumEdges() { return e; }
    public int getMaxVertexNum() { return Arrays.stream(vertices).max().getAsInt();}
    
    /**
     * Get the array of vertices
     * @return
     */
    public int[] getVertices() {
		return vertices;
	}

    /**
     * Check if the graph contains the edge u-v
     * @param int u, where u is a vertex ID (and not the index of the matrix)
     * @param int v
     * @return boolean edgeExists
     */
    public boolean containsEdge(int u, int v) {
		boolean edgeExists = false;
    		if(adj[u-1][v-1] == 1) {
    			edgeExists = true;
    		}
    		return(edgeExists);
    }

    /**
     * Get the list of neighbours of v
     * @param int v - a vertex
     * @return List<Integer> neighbours
     */
    public List<Integer> getNeighbours(int v) {
    	List<Integer> neighbours = new ArrayList<Integer>();
    		
        for(int i=0; i < vertices.length; i++) {
        	if(adj[v-1][i] == 1) {
        		neighbours.add(i+1); // index position in the adjacency matrix
        		// corresponds to the vertex number minus one, since Java arrays
        		// start from position 0
        	}
        }
		return neighbours;
    }
    
    /**
     * Count the number of edges between two lists of vertices (left & right)
     * @param List<Integer> left
     * @param List<Integer> right
     * @return int count
     * 
     * Note: if 'left' and 'right' are the sets of leaf nodes descended from an internal node r,
     * then 'count' is E_r
     */
    public int countEdgesBetween(List<Integer> left, List<Integer> right) {
    		int count = 0;
    		
    		// loop through rows given by left, 
    		// counting the number of 1's at the columns specified by right
    		for(int i : left) {
    			for(int j : right) {
    				if(adj[i-1][j-1] == 1) {
    					count++;
    				}
    			}
    		}
		return count; 	
    }
    
    /**
	 * Read in a graph from a CSV file of edges and store as an GraphAdjMatrix
	 * 
	 * @param file - CSV file with two columns; each row defines an undirected edge between two vertices,
	 * represented by their integer ID values. E.g. "3, 5" denotes an edge between vertex 3 and vertex 5.
	 * 
	 * Conventions: vertex numbering can start from 0 or 1. Suppose m is the maximum integer included in the file.
	 * Then any integer k, such that 0 < k < m, which does not appear in the file, is assumed to a singleton vertex.
	 * All singleton vertices will be included as leaves of the dendrogram and the adjacency matrix will have
	 * rows 1, 2,...,m or 1, 2,...,m+1 (depending on whether numbering started at 0 or 1).
	 */
	@DesignatedConstructor
	public static GraphAdjMatrix parseToGraphAdjMatrix(@ConstructorArg(value = "file") File file) {
		
		// First, iterate through the input file and create a set of edges and record the max and min
		// vertex numbers. Then convert to an adjacency matrix.
		
		Set<UnorderedPair<Integer, Integer>> edges = new LinkedHashSet<>();
		int e = 0; // number of edges
		int max = 0; // the maximum observed vertex number
		int min = 1; // the minimum observed vertex number
		
		// iterate through the lines in the CSV file
		for (List<String> line : BriefIO.readLines(file).splitCSV()) {
			int firstVertex = Integer.parseInt(line.get(0));
			int secondVertex = Integer.parseInt(line.get(1));
			// add an edge
			edges.add(UnorderedPair.of(firstVertex, secondVertex));
			e++;

			// update the current maximum
			if(firstVertex > max) {
				max = firstVertex; 
			}
			if(secondVertex > max) {
				max = secondVertex;
			}
			
			// update the current minimum
			if(firstVertex < min) {
				min = firstVertex;
			}
			if(secondVertex < min) {
				min = secondVertex;
			}
		}
		
		// create an array of all vertices -- including singletons 
	    int[] vertices = new int[max - min + 1];
	    int v = vertices.length;
	    
		// Check to see if the vertex numbers start from 1 or 0.
		boolean startAtZero = false; // true if vertices start at 0, false if they start at 1
		if(min == 0) {
			startAtZero = true; 
			max++; // max vertex number increased by 1
		}
		for (int i = 0; i < v; i++) {
		    vertices[i] = i+1; // start vertex numbering from 1, regardless of whether min = 0 or 1
		}
		
		// initialize the adjacency matrix
		// note: adj matrix is initialized to the size of the maximum vertex number
		// appearing in the data.
		int[][]adj = new int[max][max];
    
		// Diagnostic check:
		System.out.println("The length of the vertex array is: " + v);
		System.out.println("The max vertex number is: " + max);
    
		// iterate through the set of edges and add each to the adj matrix
    
		if(startAtZero) { 
			for(UnorderedPair<Integer, Integer> edge : edges) {
				int first = edge.getFirst();
				int second = edge.getSecond();
		
				// create binary entry in the adjacency matrix, but note that Java
				// counts array positions starting from 0.
				adj[first][second] = 1;
				adj[second][first] = 1; // adj matrix is symmetric for undirected graph
			}
			
		} else { // vertex numbers start at 1
    			for(UnorderedPair<Integer, Integer> edge : edges) {
    				int first = edge.getFirst();
    				int second = edge.getSecond();
    		
    				// create binary entry in the adjacency matrix, but note that Java
    				// counts array positions starting from 0.
    				adj[first - 1][second - 1] = 1;
    				adj[second - 1][first - 1] = 1; // adj matrix is symmetric for undirected graph

    			}
		}
        
		return new GraphAdjMatrix(vertices, e, adj, false);
	}
	
	/**
	 * Print an adjacency matrix
	 * 
	 * @param 2d integer matrix
	 */
	public static void printMatrix(int[][] matrix) {
		
		for (int i = 0; i < matrix.length; i++) {
		    for (int j = 0; j < matrix[i].length; j++) {
		        System.out.print(matrix[i][j] + " ");
		    }
		    System.out.println();
		}
	}
	
	/**
	 * Convert the adjacency matrix to a 2d array of edges. Each row contains an edge.
	 * 
	 */
	public int[][] convertToEdgeList() {
		int v = this.getMaxVertexNum();
		int[][] edges = new int[e][2];
		
		int eCounter = 0;
		
		outerloop:
		for(int r = 0; r < v; r++) {
			for(int c = 0; c < r; c++) { // need to check only lower-triangle of adj matrix, due to symmetry
				if(eCounter >= e) {
					break outerloop; // leave loop early when we've encountered all edges
				}
				if(adj[r][c] == 1) {
					edges[eCounter][0] = r+1; // add one so vertexID's start from 1
					edges[eCounter][1] = c+1;
					eCounter++;
				}
			}
		}
		return edges;
	}
	
	 public static class MutableGraphAdjMatrix extends GraphAdjMatrix {
	    /**
	     * Create a latent GraphAdjMatrix (only use for exact invariance tests)
	     * @param size
	     */
	    public MutableGraphAdjMatrix(int size) {
	    	this(vertices(size));
	    }
	    
	    public MutableGraphAdjMatrix(int [] vertices) {
	    	super(vertices, 0, 
	    			new int[Arrays.stream(vertices).max().getAsInt()][Arrays.stream(vertices).max().getAsInt()], true);
	    	// Use the max vertex number in vertices array to create the adj matrix.
	    	// This prevents issues where certain numbers are excluded from vertices, like in the case of 
	    	// graphs with singletons.
	    }
	    
	    public void clearEdges() {
	      e = 0;
	      for (int i = 0; i < vertices.length; i++)
	        for (int j = 0; j < vertices.length; j++)
	          adj[i][j] = 0;
	    }
	    
	    /**
	     * Add an undirected edge u-v
	     * @param int u
	     * @param int v
	     */
	    protected void addEdge(int u, int v) {
	        if (adj[u-1][v-1] == 0) {
	            e++;
	            adj[u-1][v-1] = 1;
	            adj[v-1][u-1] = 1;
	        }
	    }
	    
	    /**
	     * Sample a graph from the likelihood by instantiating internal node probabilities for all internal nodes r. 
	     * Theta represents the mixture parameter. rho represents the prob. parameter of the truncated Geometric likelihood.
	     * 
	     * Cases: theta = 0 => standard model
	     * 		  theta = 1 and rho = 0 => Zero-Inflated mixture model (point mass & Binomial)
	     * 		  theta = 1 and rho in (0, 1] => (truncated) Geometric & Binomial mixture model
	     * 
	     * @param Random rand
	     * @param double alpha, beta
	     * @param double theta -- 0 < theta <=1 for mixture model. theta = 0 for standard model
	     * @param double rho -- 0 < rho <= 1 for Geometric & Binomial mixture model. rho = 0 for one of the other models.
	     * @param BasicDendrogram dend
	     */
	    public void sampleGraph(Random rand, double alpha, double beta, double theta, double rho, BasicDendrogram dend) {
	      
	      // extract the leaves of dend and convert to an integer array of vertices
	      List<Integer> leaves = dend.getLeaves(dend.getRoot());
	      int[] vertices = Ints.toArray(leaves);
	      
	      clearEdges();
	      
	      // instantiate internal node probabilities
	      int n = vertices.length; // number of vertices
	      double[] p = new double[n-1]; // array containing internal node probabilities
	      
	      
	      // Get all internal tree nodes in dend
	      // Note: there are n-1 internal nodes for n vertices
	      Set<BasicTreeNode> internalNodes = dend.getInternalNodes();
	      
	      for(BasicTreeNode r: internalNodes) {
	    	  
	    	// extract the vertex IDs of the leaf nodes from the left and right subtrees rooted at r
		    List<Integer> left = dend.getLRSubtree(r, true);
		    List<Integer> right = dend.getLRSubtree(r, false);
		    // by default, r must be the LCA between the leaf nodes in the sets 'left' and 'right'
		    
		    // extract the vertex IDs of the leaf nodes from the left and right subtrees rooted at r
	        List<Integer> left = dend.getLRSubtree(r, true);
	        List<Integer> right = dend.getLRSubtree(r, false);
		    
		    if(rho == 0) { // use the standard or Zero-inflated model
		    	
		    	if(theta == 0) { //use the standard model
		    		
		    		generateBinomialEdges(Random rand, double alpha, double beta, List<Integer> left, List<Integer> right)
		    		
		    	} else { // use the Zero-inflated model
		    		
		    		// with probability 1-theta, p_r ~ Beta(a, b), and we generate edges from Binomial component
		    		if(Generators.bernoulli(rand, 1-theta)) { 
		    			generateBinomialEdges(Random rand, double alpha, double beta, List<Integer> left, List<Integer> right)
		    		}
		    		// otherwise (with prob. theta), p_r = 0, and hence we generate no edges
		    		
		    	}
		    	
		    	
		    } else { // use the Geom & Binom mixture model
		    	
		    	// with probability 1-theta, p_r ~ Beta(a, b), and we generate edges from Binomial component
	    		if(Generators.bernoulli(rand, 1-theta)) { 
	    			
	   
	    			generateBinomialEdges(Random rand, double alpha, double beta, List<Integer> left, List<Integer> right)
	    		} else { // with 
		    	
	    		}
	    		
		    }
	    	  
	      }
	      
	      
	      
	      
	      
	      
	      
	      
	      
	      
	      
	      // Create an ArrayList to store which internal nodes will follow Geom. component.
	      // This will remain empty if one of the other models is used.
	      ArrayList<Integer> geom_nodes = new ArrayList<Integer>();
	      
	      for(int r = 0; r < n-1; r++) { 
	    	  p[r] = Generators.beta(rand, alpha, beta);
	    	  
	    	  if(theta > 0) { // use one of the mixture models
	    		  
	    		  if(rho == 0) { // use the Zero-inflated model
	    			  // with prob. theta, p_r is 0
	    			  if(Generators.bernoulli(rand, theta)) { p[r] = 0; } 
	    		  } else { // use the Geom & Binom model
	    			  // with prob. theta, p_r is rho
	    			  if(Generators.bernoulli(rand, theta)) { 
	    				  p[r] = rho; 
	    				  geom_nodes.add(r); // store the indices of the internal nodes following the Geom. component
	    			  } 
	    		  }
	    		  
	    	  }
	      }

	      // get all internal tree nodes in dend
	      // Note: there are n-1 internal nodes for n vertices
	      Set<BasicTreeNode> internalNodes = dend.getInternalNodes();
	      
	      int idx = 0; // the current index for the internal node prob. p_r from the array
	      
	      for(BasicTreeNode r: internalNodes) {
	        
	        // extract the vertex IDs of the leaf nodes from the left and right subtrees rooted at r
	        List<Integer> left = dend.getLRSubtree(r, true);
	        List<Integer> right = dend.getLRSubtree(r, false);
	        
	        // by default, r must be the LCA between the leaf nodes in the sets 'left' and 'right'
	        
	        if(rho == 0) { // use the standard or zero-inflated model
	        	
	        } else { // use the geom & binomial mixture model
	        	if(r in geom_nodes) { // this internal node uses Geom component
	        		
	        	} else { // use Binomial component
	        		
	        	}
	        }
	        // BINOMIAL COMPONENT:
	        // for every pair of leaves from the left and right sets, flip a weighted-coin with prob. p_r,
	        // to determine if they will be joined by an edge
	        // Note: this code will also generate the Zero-Inflated component, since p[r] will be 0, so Bernoulli(0) will return 0 for that
	        // edge.
	        for(int u: left) {
	          for(int v: right) {
	              Boolean edge = Generators.bernoulli(rand, p[idx]);
	              if(edge) { addEdge(u, v); }
	          }
	        }
	        idx++; // move to next internal node prob.  
	      }
	      
	    }
	  }
}
