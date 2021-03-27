/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.boundary.factory;

/**
 *
 * @author winsway
 */
public class RobinBC extends AbstractBC {

    public double a, b, c, fx;
    public double M;

    /**
     *
     * @param a
     * @param b
     * @param c
     * @param fx
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
        System.out.println(" a = " + a + " b = " + b + " c =" + c + " fx = " + fx);
        return null;
    }

    @Override
    public RobinBC clone() throws CloneNotSupportedException {
        RobinBC clone = new RobinBC(this);
        return clone;
    }

}
