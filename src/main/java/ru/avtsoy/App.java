package ru.avtsoy;

import java.util.Arrays;

/**
 * Hello world!
 */
public class App {
    private static final double TAU = 1;
    private static final double A = 1;
    private static final double B = 1;
    private static final int N_X = 11;
    private static final int N_Y = 11;
    private static final double H_X = A / (N_X - 1);
    private static final double H_Y = B / (N_Y - 1);
    private static final double MU = 1 / (2 * Math.PI);
    private static final double T = Math.PI / 2;
    private static final double ALPHA = 0.0001;
    private static final int N_TIME_STEPS = 11;
    private static final double TIME_STEP = T / (N_TIME_STEPS - 1);
    private static final double[][][] f = new double[N_TIME_STEPS][N_X][N_Y];

    public static void main(String[] args) {
        System.out.println("Hello World!");
        double[][][] phi = new double[N_TIME_STEPS][N_X][N_Y];
        phi = eqForPhi(phi);
        System.out.println(Arrays.deepToString(phi[N_TIME_STEPS - 1]));
        double[][][] q = new double[N_TIME_STEPS][N_X][N_Y];
        q = eqForQ(q,phi);
        System.out.println(Arrays.deepToString(q[0]));
    }

    private static double[][][] eqForPhi(double[][][] phi) {
        for (int j = 1; j < N_X - 1; j++) {
            for (int k = 1; k < N_Y - 1; k++) {
                phi[0][j][k] = u0(j, k);
            }
        }
        phi[0] = corner(phi[0]);
        for (int j = 1; j < N_X - 1; j++) {
            for (int k = 1; k < N_Y - 1; k++) {
                phi[1][j][k] = u1(j, k) * TIME_STEP + phi[0][j][k];
            }
        }
        phi[1] = corner(phi[1]);
        for (int i = 2; i < N_TIME_STEPS; i++) {
            for (int j = 1; j < N_X - 1; j++) {
                for (int k = 1; k < N_Y - 1; k++) {
                    phi[i][j][k] = Math.pow(TIME_STEP, 2) * (f[i - 1][j][k] +
                            MU * (differenceForX(j, k, phi[i - 1]) +
                                    differenceForY(j, k, phi[i - 1]))) +
                            2 * phi[i - 1][j][k] - phi[i - 2][j][k];
                }
            }
            phi[i] = corner(phi[i]);
        }
        return phi;
    }

    private static double[][][] eqForQ(double[][][] q, double[][][] phi) {
        for (int j = 1; j < N_X - 1; j++) {
            for (int k = 1; k < N_Y - 1; k++) {
                q[N_TIME_STEPS - 1][j][k] = (phi[N_TIME_STEPS - 1][j][k] - phi[N_TIME_STEPS - 2][j][k]) / TIME_STEP - varphi1(j, k);
            }
        }
        q[N_TIME_STEPS - 1] = corner(q[N_TIME_STEPS - 1]);
        for (int j = 1; j < N_X - 1; j++) {
            for (int k = 1; k < N_Y - 1; k++) {
                q[N_TIME_STEPS - 2][j][k] = q[N_TIME_STEPS - 1][j][k] - TIME_STEP * (varphi0(j, k) - phi[N_TIME_STEPS - 1][j][k]);
            }
        }
        q[N_TIME_STEPS - 2] = corner(q[N_TIME_STEPS - 2]);
        for (int i = N_TIME_STEPS - 3; i >= 0; i--) {
            for (int j = 1; j < N_X - 1; j++) {
                for (int k = 1; k < N_Y - 1; k++) {
                    q[i][j][k] = Math.pow(TIME_STEP, 2) * (f[i + 1][j][k] +
                            MU * (differenceForX(j, k, q[i + 1]) +
                                    differenceForY(j, k, q[i + 1]))) +
                            2 * q[i + 1][j][k] - q[i + 2][j][k];
                }
            }
            q[i] = corner(q[i]);
        }
        return q;
    }

    private static double u0(int j, int k) {
        return 0;
    }

    private static double u1(int j, int k) {
        return Math.sin(j * H_X) * Math.sin(k * H_Y);
    }

    private static double varphi0(int j, int k) {
        return Math.sin(j * H_X) * Math.sin(k * H_Y);
    }

    private static double varphi1(int j, int k) {
        return 0;
    }

    private static double differenceForX(int j, int k, double[][] arr) {
        return (arr[j - 1][k] - 2 * arr[j][k] + arr[j + 1][k]) / Math.pow(H_X, 2);
    }

    private static double differenceForY(int j, int k, double[][] arr) {
        return (arr[j][k - 1] - 2 * arr[j][k] + arr[j][k + 1]) / Math.pow(H_Y, 2);
    }

    private static double[][] corner(double[][] arr) {
        for (int k = 0; k < N_Y; k++) {
            arr[0][k] = 0;
            arr[N_X - 1][k] = 0;
        }
        for (int j = 0; j < N_X; j++) {
            arr[j][0] = 0;
            arr[j][N_Y - 1] = 0;
        }
        return arr;
    }
}
