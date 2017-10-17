package ru.avtsoy;

/**
 * Hello world!
 *
 */
public class App {
    private static final double TAU = 1;
    private static final double A = 1;
    private static final double B = 1;
    private static final double H_X = 0.1;
    private static final double H_Y = 0.1;
    private static final double MU = 1/(2*Math.PI);
    private static final double T = Math.PI/2;
    private static final double ALPHA = 0.0001;
    private static final int N_X = 10;
    private static final int N_Y = 10;
    private static double[][] phi  = new double[N_X][N_Y];
    public static void main( String[] args ) {
        System.out.println( "Hello World!" );

    }
}
