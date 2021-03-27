/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.matrix.solve;

import com.winswe.matrix.Matrix;
import com.winswe.mesh.Structed2D;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月8日 下午4:25:47
 */
public class GaussSeidel extends MatrixSolver {

    public GaussSeidel(
            Matrix coe,
            double[] FI,
            Structed2D grid,
            SolverPerformance solverPerformance
    ) {
        super(coe, FI, grid, solverPerformance);
    }

    @Override
    public void solve() {

        double[] RES;

        RES = new double[FI.length];

        for (int loopi = 0; loopi < solverPerformance.getSubLoopNumber(); loopi++) {

            for (int J = 1; J <= grid.getNY(); J++) {
                for (int I = 1; I <= grid.getNX(); I++) {
                    int IJ = grid.getCellIndex(I, J);
                    int IJN = grid.getCellIndex(I, J + 1);
                    int IJS = grid.getCellIndex(I, J - 1);
                    int IJE = grid.getCellIndex(I + 1, J);
                    int IJW = grid.getCellIndex(I - 1, J);
                    RES[IJ] = coe.getS()[IJ] - coe.getAP()[IJ] * FI[IJ]
                            - coe.getAN()[IJ] * FI[IJN] - coe.getAS()[IJ] * FI[IJS]
                            - coe.getAE()[IJ] * FI[IJE] - coe.getAW()[IJ] * FI[IJW];
                    FI[IJ] = (RES[IJ] + coe.getAP()[IJ] * FI[IJ]) / coe.getAP()[IJ];
                }
            }

            //check convergence
            double RESN = 0;
            for (int resi = 0; resi < RES.length; resi++) {
                RESN = RESN + Math.abs(RES[resi]);
            }

            if (solverPerformance.getLoopNumber() == 0) {
                solverPerformance.setFirstResidual(RESN);
            }

            solverPerformance.setRelativeResidual(RESN / solverPerformance.getFirstResidual());

            if (solverPerformance.convergence()) {
                break;
            }

        }

        solverPerformance.plus();
    }

}
