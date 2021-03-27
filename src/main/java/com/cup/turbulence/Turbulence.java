/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.turbulence;

import com.cup.field.Field;
import com.cup.mesh.Mesh;
import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

/**
 *
 * @author winsw
 */
public interface Turbulence {

    final public double CAPPA = 0.41;

    final public double ELOG = 8.342;

    /**
     * 计算应变率张量，采用内推法、没有进行开根号
     *
     * @param W
     * @param mesh
     * @param S shear strain tensor
     */
    static public void cacS(Field W, Mesh mesh, double[][][] S) {
        double leta;
        double a, b;
        for (int z = 1; z <= mesh.getNZ(); ++z) {
            for (int y = 1; y <= mesh.getNY(); ++y) {
                for (int x = 1; x <= mesh.getNX(); ++x) {
                    leta = sqrt(mesh.alpha(mesh.gettPx()[x], mesh.gettPy()[y]));
                    if (x < 0.5 * W.getNewField().length) {
                        a = (W.getOldField()[x][y][z] - W.getOldField()[x - 1][y][z]) / (mesh.getDXP()[x - 1]) / leta;
                    } else {
                        a = (W.getOldField()[x + 1][y][z] - W.getOldField()[x][y][z]) / (mesh.getDXP()[x]) / leta;
                    }
                    if (y < 0.5 * W.getNewField()[0].length) {
                        b = (W.getOldField()[x][y][z] - W.getOldField()[x][y - 1][z]) / (mesh.getDYP()[y - 1]) / leta;
                    } else {
                        b = (W.getOldField()[x][y + 1][z] - W.getOldField()[x][y][z]) / (mesh.getDYP()[y]) / leta;
                    }
                    S[x][y][z] = a * a + b * b;
                }
            }
        }
    }

    /**
     * center difference
     *
     * @param W
     * @param mesh
     * @param S shear strain tensor
     */
    static public void cacScenter(Field W, Mesh mesh, double[][][] S) {
        double alpha;
        double aw, ae, bw, be;
        double a, b;
        for (int z = 1; z <= mesh.getNZ(); ++z) {
            for (int y = 1; y <= mesh.getNY(); ++y) {
                for (int x = 1; x <= mesh.getNX(); ++x) {
                    alpha = mesh.alpha(mesh.gettPx()[x], mesh.gettPy()[y]);
                    aw = (W.getOldField()[x][y][z] - W.getOldField()[x - 1][y][z])
                            / (mesh.getDXP()[x - 1]);
                    ae = (W.getOldField()[x + 1][y][z] - W.getOldField()[x][y][z])
                            / (mesh.getDXP()[x]);
                    bw = (W.getOldField()[x][y][z] - W.getOldField()[x][y - 1][z])
                            / (mesh.getDYP()[y - 1]);
                    be = (W.getOldField()[x][y + 1][z] - W.getOldField()[x][y][z])
                            / (mesh.getDYP()[y]);
                    a = 0.5 * (aw + ae);
                    b = 0.5 * (bw + be);
                    S[x][y][z] = (a * a + b * b) / alpha;
                }
            }
        }
    }

    /**
     * center difference
     *
     * @param W
     * @param mesh
     * @param S shear strain tensor
     */
    static public void cacScenter1(Field W, Mesh mesh, double[][][] S) {
        double alpha;
        double a, b;
        for (int Z = 1; Z <= mesh.getNZ(); ++Z) {
            for (int Y = 1; Y <= mesh.getNY(); ++Y) {
                for (int X = 1; X <= mesh.getNX(); ++X) {
                    alpha = mesh.alpha(mesh.gettPx()[X], mesh.gettPy()[Y]);
                    a = (W.getOldField()[X + 1][Y][Z] - W.getOldField()[X - 1][Y][Z])
                            / (mesh.getDXP()[X - 1] + mesh.getDXP()[X]);
                    b = (W.getOldField()[X][Y + 1][Z] - W.getOldField()[X][Y - 1][Z])
                            / (mesh.getDYP()[Y - 1] + mesh.getDYP()[Y]);
                    S[X][Y][Z] = (a * a + b * b) / alpha;
                }
            }
        }
    }

    static public double cacTauw(double dpdz, double radia) {
        return 0.25 * 2.0 * radia * dpdz;
    }

    static public double cacYplus(double dw, double rho, double tauw, double mu) {
        return dw * sqrt(rho * tauw) / mu;
    }

    /**
     * get the distance from wall to point.
     *
     * @param eta
     * @param xi
     * @param mesh
     * @param radia
     * @param theta1
     * @return
     */
    static public double cacDw(double eta, double xi,
            Mesh mesh, double radia, double theta1) {
        double x0, y0, hL;
        x0 = mesh.realX(eta, xi);
        y0 = mesh.realY(eta, xi);
        hL = radia - radia * cos(theta1);
        double dw
                = radia
                - sqrt(pow(x0, 2) + pow((y0 - radia + hL), 2));
        return dw;

    }

    static public double turbulenceLengthscale(double dh) {
        return 0.07 * dh;
    }

    static public double turbulenceIntensity(double Re) {
        return 0.0550 * pow(Re, -0.0407);
    }

    static public double getRe(double rho, double V, double D, double mu) {
        double re = rho * V * D / mu;
        return re;
    }

    static public double setK(double W) {
        return 1e-4 * pow(W, 2.0);
    }

    static public double setEpsilon(double K, double R) {
        double temp
                = pow(0.09, 0.75) * pow(K, 3.0 / 2.0)
                / Turbulence.turbulenceLengthscale(R * 2);
        return temp;
    }

    static public double setOmega(double K, double R) {
        double temp;
        temp
                = pow(0.09, 0.75) * pow(K, 3.0 / 2.0)
                / Turbulence.turbulenceLengthscale(R * 2);
        temp = temp / K;
        return temp;
    }

    static public double getTAU(double vis, double W1, double W, double DN) {
        return vis * (W1 - W) / DN;
    }

    /**
     *
     * @param rho
     * @param mu
     * @param CK
     * @param DN
     * @return rho * CK * DN / mu;
     */
    static public double getYPL(double rho, double mu, double CK, double DN) {
        return rho * CK * DN / mu;
    }

    /**
     *
     * @param YPL
     * @param mu
     * @return YPL * mu * CAPPA / Math.log(ELOG * YPL);
     */
    static public double getVISCW(double YPL, double mu) {
        return YPL * mu * CAPPA / Math.log(ELOG * YPL);
    }

    static public double getGEN(double TAU, double K, double DN) {
        double CMU = 0.09;
        double CMU25 = sqrt(sqrt(CMU));
        return abs(TAU) * CMU25 * sqrt(max(0, K)) / (DN * Turbulence.CAPPA);
    }

    /**
     *
     * @param K
     * @return CMU25 * sqrt(max(0, K));
     */
    static public double getCK(double K) {
        double CMU = 0.09;
        double CMU25 = sqrt(sqrt(CMU));
        return CMU25 * sqrt(max(0, K));
    }

    static public double getF1(double KSI) {
        return Math.tanh(pow(KSI, 4.0));
    }

    static public double getF2(double K, double Omega, double mu, double rho, double DN) {
        double A, B;
        double betastar = 0.09;
        A = 2.0 * Math.sqrt(K) / (betastar * DN * Omega);
        B = (500.0 * mu / rho) / (DN * DN * Omega);
        return Math.tanh(pow(Math.max(A, B), 2.0));
    }

    static public double getCDKW(double rho, double Omega, double DKW) {
        double temp;
        double sigmaw2 = 0.856;
        temp = Math.max(2.0 * rho * sigmaw2 / Omega * DKW, 1.E-10);
        return temp;
    }

    static public double getKSI(double K, double Omega, double mu, double rho,
            double DN, double CDKW) {
        double temp;
        double A, B, C;
        double betastar = 0.09;
        double sigmaw2 = 0.856;
        A = Math.sqrt(K) / (betastar * Omega * DN);
        B = (500.0 * mu) / (rho * Omega * DN * DN);
        C = 4.0 * rho * sigmaw2 * K / (CDKW * DN * DN);
        temp = Math.min(Math.max(A, B), C);
        return temp;
    }

}
