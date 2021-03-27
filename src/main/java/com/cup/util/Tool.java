/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.util;

import com.cup.field.Field;
import com.cup.mesh.Mesh;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

/**
 *
 * @author winsw
 */
public final class Tool {

    public static double clip(double value, double minClip, double maxClip) {
        return value < minClip ? minClip
                : value > maxClip ? maxClip
                        : value;
    }

    static public void setVaule(double[][][] temp, double value) {
        for (int i = 0; i < temp.length; i++) {
            for (int j = 0; j < temp[0].length; j++) {
                for (int k = 0; j < temp[0][0].length; j++) {
                    temp[i][j][k] = value;
                }
            }
        }
    }

    static public void setVaule(double[][][] temp, double[][][] value) {
        for (int i = 0; i < temp.length; i++) {
            for (int j = 0; j < temp[0].length; j++) {
                for (int k = 0; j < temp[0][0].length; j++) {
                    temp[i][j][k] = value[i][j][k];
                }
            }
        }
    }

    static public void copyAtoB(double[][][] A, double[][][] B) {
        for (int k = 0; k < A[0][0].length; ++k) {
            for (int j = 0; j < A[0].length; ++j) {
                for (int i = 0; i < A.length; ++i) {
                    B[i][j][k] = A[i][j][k];
                }
            }
        }
    }

    public static double getTotalVol(Field W, Mesh mesh) {
        double temp = 0;
        for (int z = 1; z < W.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < W.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < W.getNewField().length - 1; ++x) {
                    temp
                            = temp
                            + W.getNewField()[x][y][z]
                            * mesh.getVol(x, y, z);
                }
            }
        }
        return temp;
    }

    public static double getTotalVol(Field W, Mesh mesh, double a) {
        double temp = 0;
        for (int z = 1; z < W.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < W.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < W.getNewField().length - 1; ++x) {
                    temp
                            = temp
                            + a * mesh.getVol(x, y, z);
                }
            }
        }
        return temp;
    }

    public static double distance(double x, double y, double x1, double y1) {
        double dx, dy;
        dx = x1 - x;
        dy = y1 - y;
        return sqrt(pow(dx, 2.0) + pow(dy, 2.0));
    }

}
