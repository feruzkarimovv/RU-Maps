package rumaps;

import java.util.*;

/**
 * This class represents the information that can be attained from the Rutgers University Map.
 * 
 * The RUMaps class is responsible for initializing the network, streets, blocks, and intersections in the map.
 * 
 * You will complete methods to initialize blocks and intersections, calculate block lengths, find reachable intersections,
 * minimize intersections between two points, find the fastest path between two points, and calculate a path's information.
 * 
 * Provided is a Network object that contains all the streets and intersections in the map
 * 
 * @author Vian Miranda
 * @author Anna Lu
 */
public class RUMaps {
    
    private Network rutgers;

    /**
     * **DO NOT MODIFY THIS METHOD**
     * 
     * Constructor for the RUMaps class. Initializes the streets and intersections in the map.
     * For each block in every street, sets the block's length, traffic factor, and traffic value.
     * 
     * @param mapPanel The map panel to display the map
     * @param filename The name of the file containing the street information
     */
    public RUMaps(MapPanel mapPanel, String filename) {
        StdIn.setFile(filename);
        int numIntersections = StdIn.readInt();
        int numStreets = StdIn.readInt();
        StdIn.readLine();
        rutgers = new Network(numIntersections, mapPanel);
        ArrayList<Block> blocks = initializeBlocks(numStreets);
        initializeIntersections(blocks);

        for (Block block: rutgers.getAdjacencyList()) {
            Block ptr = block;
            while (ptr != null) {
                ptr.setLength(blockLength(ptr));
                ptr.setTrafficFactor(blockTrafficFactor(ptr));
                ptr.setTraffic(blockTraffic(ptr));
                ptr = ptr.getNext();
            }
        }
    }

    /**
     * **DO NOT MODIFY THIS METHOD**
     * 
     * Overloaded constructor for testing.
     * 
     * @param filename The name of the file containing the street information
     */
    public RUMaps(String filename) {
        this(null, filename);
    }

    /**
     * **DO NOT MODIFY THIS METHOD**
     * 
     * Overloaded constructor for testing.
     */
    public RUMaps() { 
        
    }

    /**
     * Initializes all blocks, given a number of streets.
     * the file was opened by the constructor - use StdIn to continue reading the file
     * @param numStreets the number of streets
     * @return an ArrayList of blocks
     */
    public ArrayList<Block> initializeBlocks(int numStreets) {
        // WRITE YOUR CODE HERE
        ArrayList<Block> list = new ArrayList<Block>();

        for(int i = 0; i < numStreets; i++){
        String StreetName = StdIn.readLine();
        int NumOfblocks = StdIn.readInt(); 
        StdIn.readLine();

        for(int j = 0; j < NumOfblocks; j++){
            int blockNum = StdIn.readInt();
            StdIn.readLine();
            int points = StdIn.readInt(); 
            StdIn.readLine();
            double roadSize = StdIn.readDouble();
            StdIn.readLine();
            Block block = new Block(roadSize, StreetName, blockNum);
        for(int z = 0; z < points; z++){
            int x = StdIn.readInt();
            int y = StdIn.readInt();
            StdIn.readLine();
            Coordinate coordinate = new Coordinate(x, y);
                if(z == 0){
                    block.startPoint(coordinate);
                } 
                else{
                    block.nextPoint(coordinate);
                }
            }
        list.add(block);
        }
    }
        return list; // Replace this line, it is provided so the code compiles
    }

    /**
     * This method traverses through each block and finds
     * the block's start and end points to create intersections. 
     * 
     * It then adds intersections as vertices to the "rutgers" graph if
     * they are not already present, and adds UNDIRECTED edges to the adjacency
     * list.
     * 
     * Note that .addEdge(__) ONLY adds edges in one direction (a -> b). 
     */
    public void initializeIntersections(ArrayList<Block> blocks) {
        // WRITE YOUR CODE HERE
        for (int i = 0; i < blocks.size(); i++){
            Block block = blocks.get(i);
            ArrayList<Coordinate> coordinates = block.getCoordinatePoints();
            Coordinate startP = coordinates.get(0);
            Coordinate endP = coordinates.get(coordinates.size() - 1);
            int startPIndex = rutgers.findIntersection(startP);
            int endPIndex = rutgers.findIntersection(endP);
            if(startPIndex == -1){
                Intersection intersection = new Intersection(startP);
                block.setFirstEndpoint(intersection);
                rutgers.addIntersection(intersection);
            }
            if(endPIndex == -1){
                Intersection intersection = new Intersection(endP);
                block.setLastEndpoint(intersection);
                rutgers.addIntersection(intersection);
            }
            if(startPIndex != -1){
                Intersection intersection = rutgers.getIntersections()[startPIndex];
                block.setFirstEndpoint(intersection);
            }
            if(endPIndex != -1){
                Intersection intersection = rutgers.getIntersections()[endPIndex];
                block.setLastEndpoint(intersection);
            }
            startPIndex = rutgers.findIntersection(startP);
            endPIndex = rutgers.findIntersection(endP);

            Block blockA = block.copy();
            rutgers.addEdge(startPIndex, blockA);

            Block blockB = block.copy();
            blockB.setFirstEndpoint(block.getLastEndpoint());
            blockB.setLastEndpoint(block.getFirstEndpoint());
            rutgers.addEdge(endPIndex, blockB);
        }
     }

    /**
     * Calculates the length of a block by summing the distances between consecutive points for all points in the block.
     * 
     * @param block The block whose length is being calculated
     * @return The total length of the block
     */
    public double blockLength(Block block) {
        double totalLength = 0.0;
        ArrayList<Coordinate> points = block.getCoordinatePoints();
        
        for (int i = 0; i < points.size() - 1; i++) {
            totalLength += coordinateDistance(points.get(i), points.get(i + 1));
        }
        
        return totalLength;
    }

    /**
     * Use a DFS to traverse through blocks, and find the order of intersections
     * traversed starting from a given intersection (as source).
     * 
     * Implement this method recursively, using a helper method.
     */
    public ArrayList<Intersection> reachableIntersections(Intersection source) {
        ArrayList<Intersection> visited = new ArrayList<>();
        boolean[] marked = new boolean[rutgers.getIntersections().length];
        dfs(source, visited, marked);
        return visited;
    }

    private void dfs(Intersection current, ArrayList<Intersection> visited, boolean[] marked) {
        int currentIndex = rutgers.findIntersection(current.getCoordinate());
        if (currentIndex == -1 || marked[currentIndex]) {
            return;
        }
        
        marked[currentIndex] = true;
        visited.add(current);
        
        Block ptr = rutgers.adj(currentIndex);
        while (ptr != null) {
            Intersection next = ptr.other(current);
            dfs(next, visited, marked);
            ptr = ptr.getNext();
        }
    }

    /**
     * Finds and returns the path with the least number of intersections (nodes) from the start to the end intersection.
     * 
     * - If no path exists, return an empty ArrayList.
     * - This graph is large. Find a way to eliminate searching through intersections that have already been visited.
     * 
     * @param start The starting intersection
     * @param end The destination intersection
     * @return The path with the least number of turns, or an empty ArrayList if no path exists
     */
    public ArrayList<Intersection> minimizeIntersections(Intersection start, Intersection end) {
        ArrayList<Intersection> path = new ArrayList<>();
        if (start.equals(end)) {
            path.add(start);
            return path;
        }

        int startIndex = rutgers.findIntersection(start.getCoordinate());
        int endIndex = rutgers.findIntersection(end.getCoordinate());
        
        if (startIndex == -1 || endIndex == -1) {
            return path;
        }

        boolean[] marked = new boolean[rutgers.getIntersections().length];
        Intersection[] edgeTo = new Intersection[rutgers.getIntersections().length];
        ArrayList<Intersection> queue = new ArrayList<>();
        
        queue.add(start);
        marked[startIndex] = true;
        
        while (!queue.isEmpty()) {
            Intersection current = queue.remove(0);
            int currentIndex = rutgers.findIntersection(current.getCoordinate());
            
            Block ptr = rutgers.adj(currentIndex);
            while (ptr != null) {
                Intersection next = ptr.other(current);
                int nextIndex = rutgers.findIntersection(next.getCoordinate());
                
                if (!marked[nextIndex]) {
                    marked[nextIndex] = true;
                    edgeTo[nextIndex] = current;
                    queue.add(next);
                    
                    if (next.equals(end)) {
                        // Reconstruct path
                        Intersection v = end;
                        while (v != null) {
                            path.add(v);
                            int vIndex = rutgers.findIntersection(v.getCoordinate());
                            v = edgeTo[vIndex];
                        }
                        Collections.reverse(path);
                        return path;
                    }
                }
                ptr = ptr.getNext();
            }
        }
        
        return path;
    }

    /**
     * Finds the path with the least traffic from the start to the end intersection using a variant of Dijkstra's algorithm.
     * The traffic is calculated as the sum of traffic of the blocks along the path.
     * 
     * What is this variant of Dijkstra?
     * - We are using traffic as a cost - we extract the lowest cost intersection from the fringe.
     * - Once we add the target to the done set, we're done. 
     * 
     * @param start The starting intersection
     * @param end The destination intersection
     * @return The path with the least traffic, or an empty ArrayList if no path exists
     */
    public ArrayList<Intersection> fastestPath(Intersection start, Intersection end) {
        ArrayList<Intersection> path = new ArrayList<>();
        if (start.equals(end)) {
            path.add(start);
            return path;
        }

        int startIndex = rutgers.findIntersection(start.getCoordinate());
        int endIndex = rutgers.findIntersection(end.getCoordinate());
        
        if (startIndex == -1 || endIndex == -1) {
            return path;
        }

        double[] distTo = new double[rutgers.getIntersections().length];
        Intersection[] pred = new Intersection[rutgers.getIntersections().length];
        boolean[] done = new boolean[rutgers.getIntersections().length];
        
        for (int i = 0; i < distTo.length; i++) {
            distTo[i] = Double.POSITIVE_INFINITY;
        }
        distTo[startIndex] = 0.0;
        
        ArrayList<Intersection> fringe = new ArrayList<>();
        fringe.add(start);
        
        while (!fringe.isEmpty()) {
            // Find intersection with minimum distance
            Intersection current = fringe.get(0);
            int currentIndex = rutgers.findIntersection(current.getCoordinate());
            double minDist = distTo[currentIndex];
            
            for (int i = 1; i < fringe.size(); i++) {
                Intersection v = fringe.get(i);
                int vIndex = rutgers.findIntersection(v.getCoordinate());
                if (distTo[vIndex] < minDist) {
                    current = v;
                    currentIndex = vIndex;
                    minDist = distTo[vIndex];
                }
            }
            
            fringe.remove(current);
            done[currentIndex] = true;
            
            if (current.equals(end)) {
                // Reconstruct path
                Intersection v = end;
                while (v != null) {
                    path.add(v);
                    int vIndex = rutgers.findIntersection(v.getCoordinate());
                    v = pred[vIndex];
                }
                Collections.reverse(path);
                return path;
            }
            
            Block ptr = rutgers.adj(currentIndex);
            while (ptr != null) {
                Intersection next = ptr.other(current);
                int nextIndex = rutgers.findIntersection(next.getCoordinate());
                
                if (!done[nextIndex]) {
                    double newDist = distTo[currentIndex] + ptr.getTraffic();
                    
                    if (newDist < distTo[nextIndex]) {
                        distTo[nextIndex] = newDist;
                        pred[nextIndex] = current;
                        
                        if (!fringe.contains(next)) {
                            fringe.add(next);
                        }
                    }
                }
                ptr = ptr.getNext();
            }
        }
        
        return path;
    }

    /**
     * Calculates the total length, average experienced traffic factor, and total traffic for a given path of blocks.
     * 
     * You're given a list of intersections (vertices); you'll need to find the edge in between each pair.
     * 
     * Compute the average experienced traffic factor by dividing total traffic by total length.
     *  
     * @param path The list of intersections representing the path
     * @return A double array containing the total length, average experienced traffic factor, and total traffic of the path (in that order)
     */
    public double[] pathInformation(ArrayList<Intersection> path) {
        double totalLength = 0.0;
        double totalTraffic = 0.0;
        
        if (path.size() < 2) {
            return new double[] {0.0, 0.0, 0.0};
        }
        
        for (int i = 0; i < path.size() - 1; i++) {
            Intersection current = path.get(i);
            Intersection next = path.get(i + 1);
            int currentIndex = rutgers.findIntersection(current.getCoordinate());
            
            Block ptr = rutgers.adj(currentIndex);
            while (ptr != null) {
                if (ptr.other(current).equals(next)) {
                    totalLength += ptr.getLength();
                    totalTraffic += ptr.getTraffic();
                    break;
                }
                ptr = ptr.getNext();
            }
        }
        
        double avgTrafficFactor = totalLength == 0.0 ? 0.0 : totalTraffic / totalLength;
        return new double[] {totalLength, avgTrafficFactor, totalTraffic};
    }

    /**
     * Calculates the Euclidean distance between two coordinates.
     * PROVIDED - do not modify
     * 
     * @param a The first coordinate
     * @param b The second coordinate
     * @return The Euclidean distance between the two coordinates
     */
    private double coordinateDistance(Coordinate a, Coordinate b) {
        // PROVIDED METHOD

        double dx = a.getX() - b.getX();
        double dy = a.getY() - b.getY();
        return Math.sqrt(dx * dx + dy * dy);
    }

    /**
     * **DO NOT MODIFY THIS METHOD**
     * 
     * Calculates and returns a randomized traffic factor for the block based on a Gaussian distribution.
     * 
     * This method generates a random traffic factor to simulate varying traffic conditions for each block:
     * - < 1 for good (faster) conditions
     * - = 1 for normal conditions
     * - > 1 for bad (slower) conditions
     * 
     * The traffic factor is generated with a Gaussian distribution centered at 1, with a standard deviation of 0.2.
     * 
     * Constraints:
     * - The traffic factor is capped between a minimum of 0.5 and a maximum of 1.5 to avoid extreme values.
     * 
     * @param block The block for which the traffic factor is calculated
     * @return A randomized traffic factor for the block
     */
    public double blockTrafficFactor(Block block) {
        double rand = StdRandom.gaussian(1, 0.2);
        rand = Math.max(rand, 0.5);
        rand = Math.min(rand, 1.5);
        return rand;
    }

    /**
     * Calculates the traffic on a block by the product of its length and its traffic factor.
     * 
     * @param block The block for which traffic is being calculated
     * @return The calculated traffic value on the block
     */
    public double blockTraffic(Block block) {
        // PROVIDED METHOD
        
        return block.getTrafficFactor() * block.getLength();
    }

    public Network getRutgers() {
        return rutgers;
    }
}
