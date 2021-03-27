/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.boundary;

/**
 * It is a mixing boundary condition
 *
 * @author winsway
 */
public class RobinBC {

    /**
     * the coefficient of partial phi over partial n
     */
    private double a;
    /**
     * the coefficient of phi
     */
    private double b;
    /**
     * right side constant "c"
     */
    private double c;
    /**
     * funciton
     */
    private double fx;
    /**
     * M=C/fx
     */
    private double M;

    /**
     *
     * @param a the coefficient of partial phi over partial n
     * @param b the coefficient of phi
     * @param c right side constant "c"
     * @param fx funciton
     */
    public RobinBC(double a, double b, double c, double fx) {
        this.a = a;
        this.b = b;
        this.c = c;
        this.fx = fx;
        this.M = c / fx;
    }

    /**
     *
     * @param robinBc
     */
    public RobinBC(RobinBC robinBc) {
        if (robinBc != null) {
            this.a = robinBc.a;
            this.b = robinBc.b;
            this.c = robinBc.c;
            this.fx = robinBc.fx;
            this.M = robinBc.M;
        }
    }

    public double newValue(double dx, double phiP) {
        double temp = dx * M / (a + b * dx) + a / (a + b * dx) * phiP;
        return temp;
    }

    @Override
    public String toString() {
        return (" a = " + a + " b = " + b + " c =" + c + " fx = " + fx);
    }

    @Override
    public RobinBC clone() throws CloneNotSupportedException {

        return new RobinBC(this);
    }

    public double getA() {
        return a;
    }

    public double getB() {
        return b;
    }

    public double getC() {
        return c;
    }

    public double getFx() {
        return fx;
    }

    public double getM() {
        return M;
    }

}
