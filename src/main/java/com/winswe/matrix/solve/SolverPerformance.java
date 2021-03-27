/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.matrix.solve;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月24日 下午4:27:18
 */
public class SolverPerformance {

    private final String name;

    /**
     * relative residual
     */
    private double relativeResidual;

    /**
     * sub loop number
     */
    private final int subLoopNumber;

    /**
     * convergence criterion
     */
    private final double convergenceCriterion;

    public SolverPerformance(
            String name,
            int subLoopNumber,
            double convergenceCriterion
    ) {
        this.name = name;
        this.subLoopNumber = subLoopNumber;
        this.convergenceCriterion = convergenceCriterion;
    }

    private double firstResidual;

    private int loopNumber;

    public double getRelativeResidual() {
        return relativeResidual;
    }

    public void setRelativeResidual(double relativeResidual) {
        this.relativeResidual = relativeResidual;
    }

    public int getSubLoopNumber() {
        return subLoopNumber;
    }

    public double getFirstResidual() {
        return firstResidual;
    }

    public void setFirstResidual(double firstResidual) {
        this.firstResidual = firstResidual;
    }

    public void plus() {
        this.loopNumber++;
    }

    public int getLoopNumber() {
        return loopNumber;
    }

    public boolean convergence() {
        return relativeResidual < convergenceCriterion;
    }

    @Override
    public String toString() {
        return "SolverPerformance{"
                + "name=" + name
                + ", loopNumber = " + loopNumber
                + ", firstResidual = " + firstResidual
                + ", relativeResidual=" + relativeResidual
                + ", convergenceCriterion = " + convergenceCriterion
                + ", subLoopNumber=" + subLoopNumber
                + '}';
    }

    public String getName() {
        return name;
    }

    public double getConvergenceCriterion() {
        return convergenceCriterion;
    }

}
