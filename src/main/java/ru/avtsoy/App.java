package ru.avtsoy;

import org.math.plot.Plot2DPanel;
import org.math.plot.Plot3DPanel;

import javax.swing.*;

import static org.math.array.DoubleArray.increment;

public class App {
    private static final double TAU = 1;
    private static final double A = 1;
    private static final double B = 1;
    private static final int N_X = 51;
    private static final int N_Y = 51;
    private static final double H_X = A / (N_X - 1);
    private static final double H_Y = B / (N_Y - 1);
    private static final double MU = 1 / (2 * Math.pow(Math.PI, 2));
    private static final double T = Math.PI / 2;
    private static final int N_TIME_STEPS = 110;
    private static final double TIME_STEP = T / (N_TIME_STEPS - 1);
    private static final double[][][] f = new double[N_TIME_STEPS][N_X][N_Y];

    private static class Runner {
        double[][] u0 = generateU0();
        double[][] u1 = generateU1();
        double[][][] phi = new double[N_TIME_STEPS][N_X][N_Y];
        double[][][] q = new double[N_TIME_STEPS][N_X][N_Y];

        Runner(int M, double alpha) {
            for (int l = 0; l < M; l++) {
                phi = eqForPhi(phi, u0, u1);
                q = eqForQ(q, phi);
                u0 = newU0(u0, q, alpha);
                u1 = newU1(u1, q, alpha);
            }
        }

    }


    public static void main(String[] args) {
        double[][] z = new double[N_X][N_Y];
        double[][] z0 = new double[N_X][N_Y];
        for (int j = 1; j < N_X - 1; j++) {
            for (int k = 1; k < N_Y - 1; k++) {
                z[j][k] = Math.sin(Math.PI * j * H_X) * Math.sin(Math.PI * k * H_Y);
            }
        }
        double[] errors = new double[8];
        double[] alphas = {0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001};
//        for (int i = 0; i < 8; i++) {
//            double[][] phiT_t = new double[N_X][N_Y];
//            Runner runner = new Runner(3,alphas[i]);
//            for (int j = 1; j < N_X - 1; j++) {
//                for (int k = 1; k < N_Y - 1; k++) {
//                    phiT_t[j][k] = (runner.phi[N_TIME_STEPS - 1][j][k] - runner.phi[N_TIME_STEPS - 2][j][k]) / TIME_STEP;
//                }
//            }
//            double varerror = 0;
//            double weight = 0;
//            for (int j = 1; j < N_X - 1; j++) {
//                for (int k = 1; k < N_Y - 1; k++) {
//                    varerror += Math.pow(runner.phi[N_TIME_STEPS - 1][j][k] - z[j][k], 2) * H_X * H_Y;
//                    varerror += Math.pow(phiT_t[j][k] - 0, 2) * H_X * H_Y;
//                    weight += Math.pow(z[j][k], 2) * H_X * H_Y;
//                }
//            }
//            varerror = Math.sqrt(varerror / weight);
//            errors[7-i] = varerror;
//            System.out.println(varerror + ", " + alphas[i]);
//        }
//        for (int i = 0; i < 8; i++) {
//            Runner runner = new Runner(3,alphas[i]);
//            double varerror = 0;
//            double weight = 0;
//            for (int i1 = 0; i1 < N_TIME_STEPS; i1++) {
//                for (int j = 1; j < N_X - 1; j++) {
//                    for (int k = 1; k < N_Y - 1; k++) {
//                        varerror += Math.pow(runner.phi[i1][j][k] -
//                                Math.sin(Math.PI * j * H_X) * Math.sin(Math.PI * k * H_Y) *
//                                        Math.sin((i1) * TIME_STEP), 2) * H_X * H_Y;
//                        weight += Math.pow(Math.sin(Math.PI * j * H_X) * Math.sin(Math.PI * k * H_Y) *
//                                Math.sin((i1) * TIME_STEP), 2) * H_X * H_Y;
//                    }
//                }
//            }
//            varerror = Math.sqrt(varerror / weight);
//            errors[7 - i] = varerror;
//            System.out.println(varerror + ", " + alphas[i]);
//        }

//        for (int i = 0; i < 8; i++) {
//            Runner runner = new Runner(3,alphas[i]);
//            double varerror = 0;
//            double weight = 0;
//                for (int j = 1; j < N_X - 1; j++) {
//                    for (int k = 1; k < N_Y - 1; k++) {
//                        varerror += Math.pow(runner.u0[j][k], 2) * H_X * H_Y;
//                        varerror += Math.pow(runner.u1[j][k] -
//                                Math.sin(Math.PI * j * H_X) * Math.sin(Math.PI * k * H_Y), 2) * H_X * H_Y;
//                        weight += Math.pow(Math.sin(Math.PI * j * H_X) * Math.sin(Math.PI * k * H_Y), 2) * H_X * H_Y;
//                    }
//                }
//
//            varerror = Math.sqrt(varerror / weight);
//            errors[7 - i] = varerror;
//            System.out.println(varerror + ", " + alphas[i]);
//        }
        Runner runner = new Runner(3,0.001);
        double[] x = increment(0.0, H_X, 1.0); // x = 0.0:0.1:1.0
        double[] y = increment(0.0, H_Y, 1.0);// y = 0.0:0.05:1.0
        double[] alpha_x = {-8, -7, -6, -5, -4, -3, -2, -1};
//        Plot2DPanel plot = new Plot2DPanel("SOUTH");
        Plot3DPanel plot = new Plot3DPanel("SOUTH");
        plot.addGridPlot("u0", x, y, runner.u0);
        plot.addGridPlot("u0orig", x, y, z0);
//        plot.addGridPlot("phiT", x, y, runner.phi[N_TIME_STEPS-1]);
//        plot.addGridPlot("phiT_t", x, y, phiT_t);
//        plot.addGridPlot("phi0", x, y, phi[0]);
//        plot.addGridPlot("u1", x, y, runner.u1);
//        plot.addGridPlot("u1orig", x, y, z);
//        plot.addGridPlot("phiTorig", x, y, z);
//        plot.addLinePlot("u/u0", alpha_x, errors);
//        plot.setAxisLabel(0, "alpha");
//        plot.setAxisLabel(1, "");
        JFrame frame = new JFrame("a plot panel");
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        frame.setSize(600, 600);
        frame.setContentPane(plot);
        frame.setVisible(true);
    }

    private static double[][][] eqForPhi(double[][][] phi, double[][] u0, double[][] u1) {
        for (int j = 1; j < N_X - 1; j++) {
            System.arraycopy(u0[j], 1, phi[0][j], 1, N_Y - 1 - 1);
        }
        phi[0] = corner(phi[0]);
        for (int j = 1; j < N_X - 1; j++) {
            for (int k = 1; k < N_Y - 1; k++) {
                phi[1][j][k] = u1[j][k] * TIME_STEP + phi[0][j][k];
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

    private static double[][] generateU0() {
        double[][] u0 = new double[N_X][N_Y];
        for (int j = 1; j < N_X - 1; j++) {
            for (int k = 1; k < N_Y - 1; k++) {
                u0[j][k] = 0;
            }
        }
        return u0;
    }

    private static double[][] generateU1() {
        double[][] u1 = new double[N_X][N_Y];
        for (int j = 1; j < N_X - 1; j++) {
            for (int k = 1; k < N_Y - 1; k++) {
                u1[j][k] = 0;
            }
        }
        return u1;
    }

    private static double[][] newU0(double[][] u0, double[][][] q, double alpha) {
        for (int j = 1; j < N_X - 1; j++) {
            for (int k = 1; k < N_Y - 1; k++) {
                u0[j][k] = u0[j][k] - TAU * (alpha * u0[j][k] - (q[1][j][k] - q[0][j][k]) / TIME_STEP);
            }
        }
        return u0;
    }

    private static double[][] newU1(double[][] u1, double[][][] q, double alpha) {
        for (int j = 1; j < N_X - 1; j++) {
            for (int k = 1; k < N_Y - 1; k++) {
                u1[j][k] = u1[j][k] - TAU * (alpha * u1[j][k] + q[0][j][k]);
            }
        }
        return u1;
    }

    private static double varphi0(int j, int k) {
        return Math.sin(Math.PI * j * H_X) * Math.sin(Math.PI * k * H_Y);
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
