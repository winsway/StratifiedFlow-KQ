/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.tur.model.sstkomega.tf;

import com.cup.boundary.factory.RobinBC;
import com.winswe.io.Writer;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import static java.lang.Math.abs;

/**
 * 功能： 1.整体控制；2.牛顿迭代；3.动网格；
 *
 * @author winsway
 */
public class Control {

    int old = 0, now = 1;
    /**
     * 项目名称
     */
    String Title;
    /**
     * 文件的位置
     */
    String position;
    /**
     * 物性参数
     */
    double[] testO, testW;
    /**
     * 初始压降和液位高度
     */
    double dpdz, hl;
    /**
     * 圆管半径,M
     */
    double Radia;
    /**
     * 设定入口流量，单位是m³/s
     */
    double Qoil, Qwater;
    /**
     * 设定入口的温度，单位都是K
     */
    double Toil, Twater;
    /**
     *
     */
    SSTKOmegaTF[][] X = new SSTKOmegaTF[10][2];
    RobinBC[] bc = new RobinBC[4];
    int Nx, Ny;

    void initial() {
        position = position + this.toString();
        int W = 0, K = 1, Omega = 2, T = 3;
        for (int i = 0; i < X.length; i++) {
            X[i][old] = new SSTKOmegaTF();
            X[i][old].position = position;
            X[i][old].Title = "Xold" + i + "\\";
            X[i][old].setFluid(testO, testW);
            X[i][old].setBoundary(bc);
            X[i][old].setMesh(Nx, Ny);
            //
            X[i][now] = new SSTKOmegaTF();
            X[i][now].position = position;
            X[i][now].Title = "Xnow" + i + "\\";
            X[i][now].setFluid(testO, testW);
            X[i][now].setBoundary(bc);
            X[i][now].setMesh(Nx, Ny);
        }
    }

    public Control(String position, String Title) {
        this.Title = Title;
        this.position = position;
    }

    public void setFluid(double[] testO, double[] testW) {
        this.testO = testO;
        this.testW = testW;
    }

    public void setBoundary(RobinBC[] bc) {
        this.bc = bc;
    }

    public void setMesh(int Nx, int Ny) {
        this.Nx = Nx;
        this.Ny = Ny;
    }

    public void setQ(double Qoil, double Qwater) {
        this.Qoil = Qoil;
        this.Qwater = Qwater;
    }

    public void setT(double Toil, double Twater) {
        this.Toil = Toil;
        this.Twater = Twater;
    }

    public void setRadia(double Radia) {
        this.Radia = Radia;
    }
    int index;

    public void application(int index) throws FileNotFoundException, CloneNotSupportedException {
        initial();
        if (index == 0) {
            scheme0();
        }
        if (index == 1) {
            scheme1();
        }
        if (index == 2) {
            scheme2();
        }
    }

    /**
     * 单独测试速度场，引入牛顿迭代
     *
     * @throws FileNotFoundException
     * @throws CloneNotSupportedException
     */
    void scheme0()
            throws FileNotFoundException,
            CloneNotSupportedException {
        dpdz = 3;
        hl = 0.5;
        index = 1;
        this.cacNewtonIteration(X[0][now], dpdz, hl, index);

    }

    /**
     * 单独测试温度场
     */
    void scheme1() {

    }

    /**
     * 测试瞬态温度场+速度场
     */
    void scheme2() {

    }

    void cacNewtonIteration(SSTKOmegaTF F1, double dpdz, double hl, int index)
            throws FileNotFoundException {
        double Gwater = Qwater * testW[0], Goil = Qoil * testO[0];
        PrintWriter pw;
        String dirName = F1.position + F1.toString();
        Writer.createDir(dirName);
        String title = "pressure and liquid high";
        String fileName = dirName + title + "." + "txt";
        pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(fileName)));
//        
        pw.printf("Gwater = %e\t Goil = %e\t\n", Gwater, Goil);
        pw.flush();
//        
        SSTKOmegaTF Fx1, Fy1;
        double xold = dpdz, yold = hl;
        double xnew, ynew;
        double F, G, Fy, Fx, Gy, Gx;
        double errF, errG;
        do {
//
            Fx1 = new SSTKOmegaTF();
            Fx1.setFluid(testO, testW);
            Fx1.setBoundary(bc);
            Fx1.setMesh(Nx, Ny);
            Fx1.application(xold * 1.01, yold, Radia, index);
            pw.printf("Fx1 xold =%e\t yold =%e\t\n", xold * 1.01, yold);
            pw.printf("Gwater = %e\t Goil = %e\t\n", Fx1.getGcwater(), Fx1.getGcoil());
            pw.flush();
//            
            Fy1 = new SSTKOmegaTF();
            Fy1.setFluid(testO, testW);
            Fy1.setBoundary(bc);
            Fy1.setMesh(Nx, Ny);
            Fy1.application(xold, yold * 1.01, Radia, index);
            pw.printf("Fy1 xold =%e\t yold =%e\t\n", xold, yold * 1.01);
            pw.printf("Gwater = %e\t Goil = %e\t\n", Fy1.getGcwater(), Fy1.getGcoil());
            pw.flush();
//
            F1.application(xold, yold, Radia, index);
            pw.printf("F1 xold =%e\t yold =%e\t\n", xold, yold);
            pw.printf("Gwater = %e\t Goil = %e\t\n", F1.getGcwater(), F1.getGcoil());
            pw.println();
            pw.flush();
//
            F = F1.F(Gwater);
            G = F1.G(Goil);
//            
            Fy = (Fy1.F(Gwater) - F) / (0.01 * yold);
            Gy = (Fy1.G(Goil) - G) / (0.01 * yold);
//            
            Fx = (Fx1.F(Gwater) - F) / (0.01 * xold);
            Gx = (Fx1.G(Goil) - G) / (0.01 * xold);
//            
            xnew = xold + (G * Fy - F * Gy) / (Fx * Gy - Gx * Fy);
            ynew = yold + (F * Gx - G * Fx) / (Fx * Gy - Gx * Fy);
//          
            xnew = Modified(xnew, xold, 1e-3, 10000000);
            ynew = Modified(ynew, yold, 0.1, 0.9);
            xold = xnew;
            yold = ynew;
//            
            errF = abs((F) / Gwater);
            errG = abs((G) / Goil);
        } while (errF >= 1e-3 || errG >= 1e-3);
    }

    public double Modified(double value, double old, double minClip, double maxClip) {
        double temp;
        temp
                = value < minClip ? old * 0.8
                        : value > maxClip ? old * 1.05
                                : value;
        return temp;
    }

    @Override
    public String toString() {
        return getClass().getName() + Title + "/";
    }
}
