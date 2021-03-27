/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.io;

import com.cup.system.SystemControl;

/**
 *
 * @author winsway
 */
public class DealFields {

    static public double[][][] colocatedMap(SystemControl sys, double[][][] B, String var) {
        double[][][] A;
        int k = 1;
        A = new double[sys.cell.numPx()][sys.cell.numPy()][sys.cell.numPz()];
        for (int i = 0; i < sys.cell.numPx(); i++) {
            for (int j = 0; j < sys.cell.numPy(); j++) {
                if ("U".equals(var)) {
                    if (i == 0) {
                        A[i][j][k] = B[i][j][k];
                    } else if (i == (sys.cell.numPx() - 1)) {
                        A[i][j][k] = B[i - 1][j][k];
                    } else {
                        A[i][j][k] = (B[i - 1][j][k] + B[i][j][k]) / 2;
                    }
                } else if ("V".equals(var)) {
                    if (j == 0) {
                        A[i][j][k] = B[i][j][k];
                    } else if (j == (sys.cell.numPy() - 1)) {
                        A[i][j][k] = B[i][j - 1][k];
                    } else {
                        A[i][j][k] = (B[i][j - 1][k] + B[i][j][k]) / 2;
                    }
                }
            }
        }
        return A;
    }

    static public double[][] colocatedMap(SystemControl sys, double[][] B, String var) {
        double[][] A;
        A = new double[sys.cell.numPx()][sys.cell.numPy()];
        for (int i = 0; i < sys.cell.numPx(); i++) {
            for (int j = 0; j < sys.cell.numPy(); j++) {
                if ("U".equals(var)) {
//                u
                    if (i == 0) {
                        A[i][j] = B[i][j];
                    } else if (i == (sys.cell.numPx() - 1)) {
                        A[i][j] = B[i - 1][j];
                    } else {
                        A[i][j] = (B[i - 1][j] + B[i][j]) / 2;
                    }
                } else if ("V".equals(var)) {
//              v
                    if (j == 0) {
                        A[i][j] = B[i][j];
                    } else if (j == (sys.cell.numPy() - 1)) {
                        A[i][j] = B[i][j - 1];
                    } else {
                        A[i][j] = (B[i][j - 1] + B[i][j]) / 2;
                    }
                }
            }
        }
        return A;
    }
}
