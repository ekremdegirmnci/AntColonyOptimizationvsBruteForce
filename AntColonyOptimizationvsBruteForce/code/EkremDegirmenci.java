import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;
/**
 * @author ekremdegirmenci
 * @version 1.0
 * @since Date: 13.05.2024
 * Migros Delivery: In this assignment we are trying to
 * find the shortest path by using 2 different methods (ACO and
 * Brute Force Methods).
 */

public class EkremDegirmenci {
    // Store the shortest distance found so far
    private static double shortestDistance = Double.MAX_VALUE;
    // Store the path corresponding to the shortest distance
    private static List<Integer> shortestPath = new ArrayList<>();
    // Record start time of the program to calculate total runtime
    private static final long startTime = System.currentTimeMillis();
    // Matrix to hold pheromone levels between houses
    private static double[][] pheromoneLevels;
    // Parameters for the Ant Colony Optimization algorithm
    private static double alpha = 0.8;  // Influence of pheromone on direction
    private static double beta = 2.5;   // Influence of heuristic value (distance)
    private static double evaporationRate = 0.9; // Rate at which pheromones evaporate
    private static double Q = 0.0001;   // Constant used to calculate pheromone increase
    private static int numberOfAnts = 50; // Number of ants used in the simulation
    private static int iterations = 500;  // Number of iterations for the algorithm
    private static final double initialPheromoneIntensity = 0.1; // Initial pheromone level on all paths
    static int chosenMethod = 1; // Determines the method to use: 1 for brute-force, 2 for ACO
    private static int displayMode = 1; // 1 to visualize shortest path, 2 to visualize pheromone levels



    public static void main(String[] args) {
        String filename = "input01.txt"; // Filename for input coordinates
        List<double[]> coordinates = readCoordinates(filename); // Load coordinates from file
        int n = coordinates.size(); // Number of nodes in the problem

        // Decide method based on user choice and execute
        if (chosenMethod == 1) {
            List<Integer> path = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                path.add(i);
            }
            calculatePermutations(coordinates, path, 1);
        } else if (chosenMethod == 2) {
            antColonyOptimization(coordinates, n);
        }

        long endTime = System.currentTimeMillis();
        double timeTaken = (endTime - startTime) / 1000.0;

        // Output the results of the computation
        System.out.println("Method: " + (chosenMethod == 1 ? "Brute-Force Method" : "Ant Colony Optimization"));
        System.out.println("Shortest Distance: " + String.format("%.5f", shortestDistance));
        System.out.println("Shortest Path: " + formatPath(shortestPath));
        System.out.println("Time it takes to find the shortest path: " + timeTaken + " seconds.");
        displayRoute(coordinates, shortestPath);
        StdDraw.pause(6500); // Pause before closing to view results
        System.exit(0);
    }

    /**
     * Reads coordinates from a file and stores them in a list of double arrays.
     * @param filename The name of the file containing the coordinates.
     * @return A list of coordinates where each coordinate is represented as a double array [x, y].
     */
    private static List<double[]> readCoordinates(String filename) {
        List<double[]> coordinates = new ArrayList<>();
        File file = new File(filename);
        try (Scanner scanner = new Scanner(file)) {
            while (scanner.hasNextLine()) {
                String fileLine = scanner.nextLine();
                String[] elems = fileLine.split(",");
                double x = Double.parseDouble(elems[0]);
                double y = Double.parseDouble(elems[1]);
                coordinates.add(new double[]{x, y});
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return coordinates;
    }

    /**
     * Generates all permutations of a given path to find the shortest possible route.
     * @param coordinates The list of coordinates.
     * @param path The current permutation of path indices being considered.
     * @param k The current position in the path to permute.
     */
    private static void calculatePermutations(List<double[]> coordinates, List<Integer> path, int k) {
        if (k == path.size()) {
            double distance = calculateTotalDistance(coordinates, path);
            if (distance < shortestDistance) {
                shortestDistance = distance;
                shortestPath = new ArrayList<>(path);
            }
        } else {
            for (int i = k; i < path.size(); i++) {
                Collections.swap(path, k, i);
                calculatePermutations(coordinates, path, k + 1);
                Collections.swap(path, k, i);
            }
        }
    }

    /**
     * Calculates the total distance for a given path.
     * @param coordinates List of coordinates.
     * @param path List of indices representing the path.
     * @return Total distance of the path.
     */
    private static double calculateTotalDistance(List<double[]> coordinates, List<Integer> path) {
        double totalDistance = 0.0;
        for (int i = 0; i < path.size() - 1; i++) {
            int from = path.get(i);
            int to = path.get(i + 1);
            totalDistance += distance(coordinates.get(from), coordinates.get(to));
        }
        // Ensure returning to the starting point is considered
        totalDistance += distance(coordinates.get(path.get(path.size() - 1)), coordinates.get(0));
        return totalDistance;
    }

    /**
     * Computes the Euclidean distance between two points.
     * @param a The first point [x1, y1].
     * @param b The second point [x2, y2].
     * @return The distance between point a and b.
     */
    private static double distance(double[] a, double[] b) {
        return Math.sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]));
    }

    /**
     * Formats the path for output by converting index to 1-based and creating a readable string.
     * @param path The path of indices.
     * @return A string representation of the path in 1-based index format.
     */
    private static String formatPath(List<Integer> path) {
        StringBuilder stringB = new StringBuilder("[");
        for (int i = 0; i < path.size(); i++) {
            stringB.append(path.get(i) + 1);
            if (i < path.size() - 1) {
                stringB.append(", ");
            } else {
                stringB.append(", 1]");
            }
        }
        return stringB.toString();
    }

    /**
     * Implements the Ant Colony Optimization algorithm to find the shortest path between points.
     * @param coordinates The coordinates of the locations.
     * @param n The number of locations.
     */
    private static void antColonyOptimization(List<double[]> coordinates, int n) {
        firstPheromones(n);  // Initialize pheromones
        Random rand = new Random();
        for (int t = 0; t < iterations; t++) {
            List<List<Integer>> allPath = new ArrayList<>();
            double[][] localPheromoneChanged = new double[n][n];
            for (int ant = 0; ant < numberOfAnts; ant++) {
                List<Integer> path = createPaths(coordinates, rand, n);
                allPath.add(path);
                double pathDistance = calculateTotalDistance(coordinates, path);
                double deltaPheromone = Q / pathDistance;
                for (int i = 0; i < path.size() - 1; i++) {
                    localPheromoneChanged[path.get(i)][path.get(i + 1)] += deltaPheromone;
                }
            }
            changePheromones(localPheromoneChanged, n);
            evaluatePaths(coordinates, allPath);
        }
    }

    /**
     * Initializes the pheromone levels for all paths.
     * @param n The number of locations.
     */
    private static void firstPheromones(int n) {
        pheromoneLevels = new double[n][n];
        for (double[] row : pheromoneLevels) {
            Arrays.fill(row, initialPheromoneIntensity); // Set initial pheromone level
        }
    }

    /**
     * Generates a path for an ant based on pheromone levels and heuristic information.
     * @param coordinates The coordinates of locations.
     * @param random A random number generator for probabilistic decisions.
     * @param n The number of locations.
     * @return A complete path starting and ending at Migros.
     */
    private static List<Integer> createPaths(List<double[]> coordinates, Random random, int n) {
        List<Integer> path = new ArrayList<>();
        path.add(0); // Start at Migros.

        List<Integer> emptyNodes = new ArrayList<>();
        for (int i = 1; i < n; i++) {
            emptyNodes.add(i);
        }

        while (!emptyNodes.isEmpty()) {
            int lastNode = path.get(path.size() - 1);
            int nextNode = selectNextNode(lastNode, emptyNodes, coordinates, random);
            path.add(nextNode);
            emptyNodes.remove(Integer.valueOf(nextNode));
        }
        return path;
    }

    /**
     * Selects the next point for an ant to visit based on the pheromone level and distance.
     * @param lastNode The last point visited.
     * @param emptyNodes The list of points not yet visited.
     * @param coordinates The coordinates of all points.
     * @param random A random number generator to help with probabilistic selection.
     * @return The index of the next point to visit.
     */
    private static int selectNextNode(int lastNode, List<Integer> emptyNodes, List<double[]> coordinates, Random random) {
        double[] probabilities = new double[emptyNodes.size()];
        double sumOfProbabilities = 0;
        for (int i = 0; i < emptyNodes.size(); i++) {
            int node = emptyNodes.get(i);
            double edgeValue = Math.pow(pheromoneLevels[lastNode][node], alpha) * Math.pow(1.0 / distance(coordinates.get(lastNode), coordinates.get(node)), beta);
            probabilities[i] = edgeValue;
            sumOfProbabilities += edgeValue;
        }
        for (int i = 0; i < probabilities.length; i++) {
            probabilities[i] /= sumOfProbabilities;
        }
        return emptyNodes.get(selectIndexBasedOnProbability(probabilities, random));
    }

    /**
     * Selects an index from a list of probabilities.
     * @param probabilities An array of probabilities corresponding to each possible choice.
     * @param random A random number generator for selection.
     * @return The selected index based on the probabilities provided.
     */
    private static int selectIndexBasedOnProbability(double[] probabilities, Random random) {
        double randomValue = random.nextDouble();
        double cumulativeProbability = 0.0;
        for (int i = 0; i < probabilities.length; i++) {
            cumulativeProbability += probabilities[i];
            if (cumulativeProbability > randomValue) {
                return i;
            }
        }
        return probabilities.length - 1; // Fallback to the last index if no other is chosen
    }

    /**
     * Updates the pheromone levels across all paths.
     * @param localPheromoneUpdate The changes to pheromone levels determined by the paths ants took.
     * @param n The number of locations.
     */
    private static void changePheromones(double[][] localPheromoneUpdate, int n) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                pheromoneLevels[i][j] = pheromoneLevels[i][j] * evaporationRate + localPheromoneUpdate[i][j];
            }
        }
    }

    /**
     * Evaluates all generated paths to find the shortest one.
     * @param coords The coordinates of all points.
     * @param allPaths A list of all paths generated in the current iteration.
     */
    private static void evaluatePaths(List<double[]> coords, List<List<Integer>> allPaths) {
        for (List<Integer> path : allPaths) {
            double distance = calculateTotalDistance(coords, path);
            if (distance < shortestDistance) {
                shortestDistance = distance;
                shortestPath = new ArrayList<>(path);
            }
        }
    }

    /**
     * Visualizes the route on a graphical interface using StdDraw.
     * @param coordinates The coordinates of all points.
     * @param optimalPath The optimal path found.
     */
    public static void displayRoute(List<double[]> coordinates, List<Integer> optimalPath) {
        StdDraw.enableDoubleBuffering();
        int canvasWidth = 600;
        int canvasHeight = 600;
        StdDraw.setCanvasSize(canvasWidth, canvasHeight);
        StdDraw.setXscale(0, 1);
        StdDraw.setYscale(0, 1);
        StdDraw.clear(StdDraw.WHITE);

        // Depending on the display mode, draw either the path or pheromone levels
        if (displayMode == 1) {
            StdDraw.setPenColor(StdDraw.BLACK);
            StdDraw.setPenRadius(0.005);
            for (int i = 0; i < optimalPath.size() - 1; i++) {
                int from = optimalPath.get(i);
                int to = optimalPath.get(i + 1);
                double[] fromCoord = coordinates.get(from);
                double[] toCoord = coordinates.get(to);
                StdDraw.line(fromCoord[0], fromCoord[1], toCoord[0], toCoord[1]);
            }
            double[] startCoord = coordinates.get(optimalPath.get(0));
            double[] endCoord = coordinates.get(optimalPath.get(optimalPath.size() - 1));
            StdDraw.line(endCoord[0], endCoord[1], startCoord[0], startCoord[1]);
        } else if (displayMode == 2) {
            double maxPheromone = maxPheromoneLevel();
            for (int i = 0; i < pheromoneLevels.length; i++) {
                for (int j = 0; j < pheromoneLevels[i].length; j++) {
                    double intensity = pheromoneLevels[i][j] / maxPheromone;
                    StdDraw.setPenRadius(0.003 + 0.01 * intensity);
                    StdDraw.setPenColor(new Color(0, 0, 0, (int) (255 * intensity)));
                    StdDraw.line(coordinates.get(i)[0], coordinates.get(i)[1], coordinates.get(j)[0], coordinates.get(j)[1]);
                }
            }
            if (!shortestPath.isEmpty()) {
                int lastNodeIndex = shortestPath.get(shortestPath.size() - 1);
                int firstNodeIndex = shortestPath.get(0);
                double[] startingCoordinate = coordinates.get(firstNodeIndex);
                double[] endingCoordinate = coordinates.get(lastNodeIndex);
                StdDraw.setPenColor(StdDraw.BLACK); 
                StdDraw.setPenRadius(0.01); // Slightly thicker line for the return path
                StdDraw.line(endingCoordinate[0], endingCoordinate[1], startingCoordinate[0], startingCoordinate[1]);
            }
        }

        // Draw all locations, change Migros color based on the chosen method
        for (int i = 0; i < coordinates.size(); i++) {
            double[] coordinate = coordinates.get(i);
            if (i == 0 && chosenMethod == 2) {
                StdDraw.setPenColor(StdDraw.LIGHT_GRAY);
            } else {
                StdDraw.setPenColor(i == 0 ? StdDraw.PRINCETON_ORANGE : StdDraw.LIGHT_GRAY);
            }
            StdDraw.filledCircle(coordinate[0], coordinate[1], 0.02);
            StdDraw.setPenColor(StdDraw.BLACK);
            StdDraw.text(coordinate[0], coordinate[1], String.valueOf(i + 1));
        }

        StdDraw.show();
    }
    
    private static double maxPheromoneLevel() {
        double maxLevel = 0;
        for (double[] row : pheromoneLevels) {
            for (double level : row) {
                if (level > maxLevel) {
                    maxLevel = level;
                }
            }
        }
        return maxLevel;
    }
}
