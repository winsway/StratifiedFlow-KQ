/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.matrix.solve;

import com.winswe.matrix.Matrix;
import com.winswe.mesh.Structed2D;
import static com.winswe.util.PARAMETER.SMALL;
import static java.lang.Math.abs;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月1日 上午5:20:50
 */
public class SIPSOL extends MatrixSolver {

    private final double ALFA = 0.92;

    public SIPSOL(
            Matrix coe,
            double[] FI,
            Structed2D grid,
            SolverPerformance solverPerformance
    ) {
        super(coe, FI, grid, solverPerformance);
    }

    @Override
    public void solve() {
        double[] LW, LS, LPR;
        double[] UE, UN, RES;

        LW = new double[FI.length];
        LS = new double[FI.length];
        LPR = new double[FI.length];
        UE = new double[FI.length];
        UN = new double[FI.length];
        RES = new double[FI.length];

        double P1, P2;

        for (int J = 1; J <= grid.getNY(); J++) {
            for (int I = 1; I <= grid.getNX(); I++) {
                int IJ = grid.getCellIndex(I, J);
                int IJW = grid.getCellIndex(I - 1, J);
                int IJS = grid.getCellIndex(I, J - 1);
                LW[IJ] = coe.getAW()[IJ] / (1.0 + ALFA * UN[IJW] + SMALL);
                LS[IJ] = coe.getAS()[IJ] / (1.0 + ALFA * UE[IJS] + SMALL);
                P1 = ALFA * LW[IJ] * UN[IJW];
                P2 = ALFA * LS[IJ] * UE[IJS];
                LPR[IJ]
                        = 1.0 / (coe.getAP()[IJ] + P1 + P2
                        - LW[IJ] * UE[IJW] - LS[IJ] * UN[IJS] + SMALL);
                UN[IJ] = (coe.getAN()[IJ] - P1) * LPR[IJ];
                UE[IJ] = (coe.getAE()[IJ] - P2) * LPR[IJ];
            }
        }

        for (int loopi = 0; loopi < solverPerformance.getSubLoopNumber(); loopi++) {
            double RESN = 0.;
            for (int I = 1; I <= grid.getNX(); I++) {
                for (int J = 1; J <= grid.getNY(); J++) {
                    int IJ = grid.getCellIndex(I, J);
                    int IJN = grid.getCellIndex(I, J + 1);
                    int IJS = grid.getCellIndex(I, J - 1);
                    int IJE = grid.getCellIndex(I + 1, J);
                    int IJW = grid.getCellIndex(I - 1, J);
                    RES[IJ]
                            = coe.getS()[IJ] - coe.getAP()[IJ] * FI[IJ]
                            - coe.getAN()[IJ] * FI[IJN] - coe.getAS()[IJ] * FI[IJS]
                            - coe.getAE()[IJ] * FI[IJE] - coe.getAW()[IJ] * FI[IJW];
                }
            }

            //FORWARD SUBSTITUTION
            for (int I = 1; I <= grid.getNX(); I++) {
                for (int J = 1; J <= grid.getNY(); J++) {
                    int IJ = grid.getCellIndex(I, J);
                    int IJS = grid.getCellIndex(I, J - 1);
                    int IJW = grid.getCellIndex(I - 1, J);

                    RESN = RESN + abs(RES[IJ]);

                    RES[IJ]
                            = (RES[IJ]
                            - LS[IJ] * RES[IJS]
                            - LW[IJ] * RES[IJW])
                            * LPR[IJ];
                }
            }

            //STORE INITIAL RESIDUAL SUM FOR CHECKING CONV. OF OUTER ITER.
            if (solverPerformance.getLoopNumber() == 0) {
                solverPerformance.setFirstResidual(RESN);
            }
            solverPerformance.setRelativeResidual(RESN / solverPerformance.getFirstResidual());

            //BACK SUBSTITUTION AND CORRECTION
            for (int I = grid.getNX(); I >= 1; I--) {
                for (int J = grid.getNY(); J >= 1; J--) {
                    int IJ = grid.getCellIndex(I, J);
                    int IJN = grid.getCellIndex(I, J + 1);
                    int IJE = grid.getCellIndex(I + 1, J);
                    RES[IJ]
                            = RES[IJ]
                            - UN[IJ] * RES[IJN]
                            - UE[IJ] * RES[IJE];

                    FI[IJ] = FI[IJ] + RES[IJ];

                }
            }
            //CHECK CONVERGENCE OF INNER ITERATIONS
            if (solverPerformance.convergence()) {
                break;
            }

        }

        solverPerformance.plus();
    }

}
