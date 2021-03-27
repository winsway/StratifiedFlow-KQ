/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.util;

import static java.lang.Math.log;
import static java.lang.Math.sqrt;

/**
 *
 * @author winsway
 */
public final class Math {

    public static double arcosh(double x) {
        try {
            if (x >= 1) {
                return log(x + sqrt(x * x - 1.0));
            }
        } catch (Exception e) {
            System.out.println("x must >= 1");
        }
        return 0;
    }

    /**
     * summer = sum(Math.pow(temp[i][j], times))
     *
     * @param temp 二维数组
     * @param times 幂次，n
     * @return 返回求和的值
     */
    static public double sum(double[][][] temp, int times) {
        double summer = 0.0;
        for (int i = 0; i < temp.length; i++) {
            for (int j = 0; j < temp[0].length; j++) {
                for (int k = 0; j < temp[0][0].length; j++) {
                    summer = summer + java.lang.Math.pow(temp[i][j][k], times);
                }
            }
        }
        return summer;
    }
}
