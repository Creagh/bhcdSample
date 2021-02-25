package bhcd;

public class Pair {
    private int[][] array1;
    private double[][] array2;
    
    public Pair(int[][] array1, double[][] array2)
    {
        this.array1 = array1;
        this.array2 = array2;

    }
    public int[][] getArray1() { return array1; }
    public double[][] getArray2() { return array2; }
}