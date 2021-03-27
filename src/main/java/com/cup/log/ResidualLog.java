/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.log;

/**
 * <li>残差的一系列方法</li>
 * <li>1.矩阵的无穷范数、2范数、1范数</li>
 * <li>2.相对残差</li>
 * <li>3.绝对残差</li>
 *
 * @author winsway
 */
public class ResidualLog {

    public double Initial_residual, Final_residual;
    public double Relative_residual;
    public double errInite, err2, err1;
    double irInifit, ir2, ir1;
    double frInifit, fr2, fr1;

    public int No_Iteration;

    public void residualDeal(double vNew, double vOld) {
        this.inifiteNorm(vNew, vOld);
        this.norm2(vNew, vOld);
        this.norm1(vNew, vOld);
    }

    public double getMaxResidual() {
        double max = (frInifit > fr2) ? frInifit : fr2;
        max = (max > fr1) ? max : fr1;
        Final_residual = max;
        return max;
    }

    public void getResidual() {
        if (No_Iteration == 1) {
            irInifit = errInite;
            ir2 = err2;
            ir1 = err1;
        }
//        frInifit = errInite / irInifit;
//        fr2 = err2 / ir2;
//        fr1 = err1 / ir1;

        frInifit = errInite;
        fr2 = err2;
        fr1 = err1;

    }

    public void intialErro() {
        errInite = 0;
        err2 = 0;
        err1 = 0;
    }

    public void inifiteNorm(double vNew, double vOld) {
        if (Math.abs(vNew - vOld) > errInite) {
            errInite = Math.abs(vNew - vOld);
        }
    }

    public void norm2(double vNew, double vOld) {
        err2 = err2 + Math.pow((vNew - vOld), 2);
    }

    public void norm1(double vNew, double vOld) {
        err1 = err1 + Math.abs(vNew - vOld);
    }

    public void outPut() {
        System.out.printf("inifiteNorm=%e   Norm2=%e   Norm1=%e\n", frInifit, fr2, fr1);
    }

    public void setZero() {
        this.Initial_residual = 0;
        this.Final_residual = 0;
        this.No_Iteration = 0;
    }

}
