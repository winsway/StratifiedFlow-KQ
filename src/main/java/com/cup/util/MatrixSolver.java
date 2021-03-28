/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.util;

import com.cup.system.SystemControl;
import com.cup.field.CoeffMatrix;
import com.cup.log.ResidualLog;

/**
 * solve Matrix
 *
 * @author winsway
 */
public class MatrixSolver {

    /**
     *
     * @param phi
     * @param coe
     * @param sys
     * @return
     */
    static public ResidualLog SIPSOL(double[][][] phi, CoeffMatrix coe, SystemControl sys) {
        ResidualLog resLog = new ResidualLog();
        double alpha = 0.92;
        double URF = 0.7;
        double URFFI = 1. / URF;
        double temp;
// need to complete
        for (int i = 1; i < phi.length - 1; i++) {
            for (int j = 1; j < phi[0].length - 1; j++) {
                for (int k = 1; k < phi[0][0].length - 1; k++) {

                }
            }
        }
        return resLog;
    }

    static public ResidualLog gaussSeidel(double[][][] phi, CoeffMatrix coe, SystemControl sys) {
        ResidualLog resLog = new ResidualLog();
        double URF = 0.9;
        double URFFI = 1. / URF;
        double temp;
        do {
            resLog.intialErro();
            for (int z = 1; z < phi[0][0].length - 1; ++z) {
                for (int y = 1; y < phi[0].length - 1; ++y) {
                    for (int x = 1; x < phi.length - 1; ++x) {
                        temp = (coe.getAw(x, y, z) * phi[x - 1][y][z]
                                + coe.getAe(x, y, z) * phi[x + 1][y][z]
                                + coe.getAs(x, y, z) * phi[x][y - 1][z]
                                + coe.getAn(x, y, z) * phi[x][y + 1][z]
                                + coe.getAb(x, y, z) * phi[x][y][z - 1]
                                + coe.getAt(x, y, z) * phi[x][y][z + 1]
                                + coe.getB(x, y, z)
                                + (URFFI - 1.0) * coe.getAp(x, y, z) * phi[x][y][z])
                                / (URFFI * coe.getAp(x, y, z));
                        resLog.residualDeal(phi[x][y][z], temp);
                        phi[x][y][z] = temp;
                    }
                }
            }
            resLog.No_Iteration++;
            resLog.getResidual();
//            System.out.println("resLog.getResidual();=" + resLog.getMaxResidual());
        } while (resLog.getMaxResidual() > sys.fvSolution.getRtol()
                && resLog.No_Iteration <= sys.fvSolution.getMaxIter());
        return resLog;
    }

    static public ResidualLog gaussSeidel(double[][][] phi, CoeffMatrix coe, SystemControl sys, int num) {
        ResidualLog resLog = new ResidualLog();
        double URF = 1.0;
        double URFFI = 1. / URF;
        double temp;
        do {
            resLog.intialErro();
            for (int z = 1; z < phi[0][0].length - 1; ++z) {
                for (int y = 1; y < phi[0].length - 1; ++y) {
                    for (int x = 1; x < phi.length - 1; ++x) {
                        temp = (coe.getAw(x, y, z) * phi[x - 1][y][z]
                                + coe.getAe(x, y, z) * phi[x + 1][y][z]
                                + coe.getAs(x, y, z) * phi[x][y - 1][z]
                                + coe.getAn(x, y, z) * phi[x][y + 1][z]
                                + coe.getAb(x, y, z) * phi[x][y][z - 1]
                                + coe.getAt(x, y, z) * phi[x][y][z + 1]
                                + coe.getB(x, y, z)
                                + (URFFI - 1.0) * coe.getAp(x, y, z) * phi[x][y][z])
                                / (URFFI * coe.getAp(x, y, z));
                        resLog.residualDeal(phi[x][y][z], temp);
                        phi[x][y][z] = temp;
                    }
                }
            }
            resLog.No_Iteration++;
            resLog.getResidual();
//            System.out.println("resLog.getResidual();=" + resLog.getMaxResidual());
        } while (resLog.getMaxResidual() > sys.fvSolution.getRtol()
                && resLog.No_Iteration < num);

        return resLog;
    }

}
