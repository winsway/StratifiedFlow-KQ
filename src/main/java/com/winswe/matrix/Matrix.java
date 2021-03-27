/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.matrix;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月1日 上午4:58:20
 */
public class Matrix {

    private double[] AE, AW, AN, AS, AT, AB;
    private double[] AR, AL;
    private double[] S;
    private double[] AP, APR;

    private final int maxI, maxJ;
    private final int size;

    public Matrix(int I, int J) {
        this.maxI = I;
        this.maxJ = J;
        size = (I + 2) * (J + 2);
    }

    public void intialMatrix() {
        AE = new double[size];
        AW = new double[size];
        AN = new double[size];
        AS = new double[size];
        AP = new double[size];
        S = new double[size];
    }

    public double[] getAE() {
        return AE;
    }

    public double[] getAW() {
        return AW;
    }

    public double[] getAN() {
        return AN;
    }

    public double[] getAS() {
        return AS;
    }

    public double[] getAT() {
        return AT;
    }

    public double[] getAB() {
        return AB;
    }

    public double[] getAR() {
        return AR;
    }

    public double[] getAL() {
        return AL;
    }

    public double[] getS() {
        return S;
    }

    public double[] getAP() {
        return AP;
    }

    public double[] getAPR() {
        return APR;
    }

    public int getSize() {
        return size;
    }

    public int getMaxI() {
        return maxI;
    }

    public int getMaxJ() {
        return maxJ;
    }

}
