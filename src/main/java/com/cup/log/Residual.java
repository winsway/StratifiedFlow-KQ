/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.log;

import static java.lang.Math.max;

/**
 *
 * @author winsw
 */
public interface Residual {

    /**
     * 无穷范数
     *
     * @param a
     * @param b
     * @return
     */
    static public double getRes(double[][][] a, double[][][] b) {
        double temp = 0;
        double M = 0;
        for (int z = 1; z < a[0][0].length - 1; ++z) {
            for (int y = 1; y < a[0].length - 1; ++y) {
                for (int x = 1; x < a.length - 1; ++x) {
                    temp = temp + Math.abs(a[x][y][z] - b[x][y][z]);
                }
                M = max(temp, M);
                temp = 0;
            }
        }
        return M;
    }
}
