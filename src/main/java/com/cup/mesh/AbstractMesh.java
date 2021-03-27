/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.mesh;

import com.cup.boundary.factory.Label;

/**
 * <li>1.网格初始化（目前仅处理了均分网格）</li>
 * <li>2.控制体容积</li>
 * <li>3.控制体各个面积和距离</li>
 *
 * @author winsway
 */
public abstract class AbstractMesh implements Mesh {

    @Override
    public Mesh fineMesh(Mesh cellf) {
        return null;
    }

    /**
     * 单个网格大小
     */
    protected double dx, dy, dz;
    /**
     * 几何域的大小
     */
    protected double LX, LY, LZ;

    /**
     * 网格线的位置
     */
    protected double[] tx, ty, tz;

    /**
     * 计算节点的位置
     */
    protected double[] tPx, tPy, tPz;

    /**
     * 不同方向的网格数目
     */
    protected int NX, NY, NZ;

    /**
     * 同位网格点的数目
     */
    protected int numPx, numPy, numPz;

    /**
     * X方向P的网格线距离；
     */
    protected double[] DXU;
    /**
     * Y方向P的网格线距离；
     */
    protected double[] DYV;
    /**
     * Z方向P的网格线距离；
     */
    protected double[] DZW;
    /**
     * P在X方向两节点之间的距离；
     */
    protected double[] DXP;
    /**
     * P在Y方向两节点之间的距离；
     */
    protected double[] DYP;
    /**
     * P在Z方向两节点之间的距离；
     */
    protected double[] DZP;

    /**
     *
     * @param X
     * @param Y
     * @param Z
     * @param nx
     * @param ny
     * @param nz
     */
    public AbstractMesh(double X, double Y, double Z, int nx, int ny, int nz) {
        this.LX = X;
        this.LY = Y;
        this.LZ = Z;
        this.NX = nx;
        this.NY = ny;
        this.NZ = nz;
        this.numPx = nx + 2;
        this.numPy = ny + 2;
        this.numPz = nz + 2;
        //网格位置    
        this.tx = new double[numPx - 1];
        this.ty = new double[numPy - 1];
        this.tz = new double[numPz - 1];
        this.tPx = new double[numPx];
        this.tPy = new double[numPy];
        this.tPz = new double[numPz];
        dx = X / nx;
        dy = Y / ny;
        dz = Z / nz;
    }

    /**
     *
     * @param X
     * @param Y
     * @param nx
     * @param ny
     */
    public AbstractMesh(double X, double Y, int nx, int ny) {
        this(X, Y, 1, nx, ny, 1);
    }

    /**
     * 初始化网格点
     */
    @Override
    public void blockMesh() {
        throw new UnsupportedOperationException();
    }

    @Override
    public double[] getTx() {
        return tx;
    }

    @Override
    public double[] getTy() {
        return ty;
    }

    @Override
    public double[] getTz() {
        return tz;
    }

    @Override
    public double[] gettPx() {
        return tPx;
    }

    @Override
    public double[] gettPy() {
        return tPy;
    }

    @Override
    public double[] gettPz() {
        return tPz;
    }

    @Override
    public double[] getDXU() {
        return DXU;
    }

    @Override
    public double[] getDYV() {
        return DYV;
    }

    @Override
    public double[] getDZW() {
        return DZW;
    }

    @Override
    public double[] getDXP() {
        return DXP;
    }

    @Override
    public double[] getDYP() {
        return DYP;
    }

    @Override
    public double[] getDZP() {
        return DZP;
    }

    @Override
    public int getNX() {
        return NX;
    }

    @Override
    public int getNY() {
        return NY;
    }

    @Override
    public int getNZ() {
        return NZ;
    }

    @Override
    public int numPx() {
        return numPx;
    }

    @Override
    public int numPy() {
        return numPy;
    }

    @Override
    public int numPz() {
        return numPz;
    }

    @Override
    public double getVol(Label flag) {
        int x, y;
        x = flag.i;
        y = flag.j;
        return DXU[x - 1] * DYV[y - 1] * J(tPx[x], tPy[y]);
    }

    @Override
    public abstract double getVol(int x, int y, int z);

    @Override
    public abstract Mesh clone();

    protected void setDistance(int MaxPx, double[] tx, double[] tPx, double L, double NUM) {
        double dxx = L / (NUM);
        for (int i = 0; i < MaxPx; i++) {
            if (i < (MaxPx - 1)) {
                tx[i] = dxx * i;
            }
            tPx[i] = dxx / 2 + (i - 1) * dxx - Math.min((i - 1), 0.0) * dxx / 2 - (Math.max(i, NUM) - NUM) * dxx / 2;
        }
    }

}
