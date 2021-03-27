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
 * @date 2021年3月24日 下午4:16:36
 */
public abstract class MatrixSolver {

    protected final Matrix coe;
    protected final Structed2D grid;
    protected final double[] FI;
    protected final SolverPerformance solverPerformance;

    public MatrixSolver(
            Matrix coe,
            double[] FI,
            Structed2D grid,
            SolverPerformance solverPerformance) {
        this.coe = coe;
        this.grid = grid;
        this.FI = FI;
        this.solverPerformance = solverPerformance;
    }

    public abstract void solve();

    /**
     *
     * @param name solver name
     * @param coe Matrix Coefficient
     * @param FI field
     * @param grid mesh
     * @param solverPerformance solver performance
     * @return matrix solver
     */
    public static MatrixSolver factory(
            String name,
            Matrix coe,
            double[] FI,
            Structed2D grid,
            SolverPerformance solverPerformance) {
        if ("GaussSeidel".equals(name)) {
            return new GaussSeidel(coe, FI, grid, solverPerformance);
        } else if ("SIPSOL".equals(name)) {
            return new SIPSOL(coe, FI, grid, solverPerformance);
        } else {
            throw new ArithmeticException("without the solver name");
        }

    }

}
