/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.mesh;

import com.cup.boundary.factory.Label;

/**
 * 网格方法接口
 *
 * @author winsway
 */
public interface Mesh {

    /**
     * 初始化网格点
     */
    public void blockMesh();

    /**
     *
     *
     * @return 得到X方向网格线的坐标
     */
    public double[] getTx();

    /**
     *
     * @return 得到Y方向网格线的坐标
     */
    public double[] getTy();

    public double[] getTz();

    public double[] gettPx();

    public double[] gettPy();

    public double[] gettPz();

    public double[] getDXU();

    public double[] getDYV();

    public double[] getDZW();

    public double[] getDXP();

    public double[] getDYP();

    public double[] getDZP();

    /**
     *
     *
     * @return 得到X方向的网格数目
     */
    public int getNX();

    /**
     *
     *
     * @return 得到Y方向的网格数目
     */
    public int getNY();

    public int getNZ();

    /**
     *
     * @return X方向的点数
     */
    public int numPx();

    /**
     *
     * @return Y方向的点数
     */
    public int numPy();

    /**
     *
     * @return Z方向的点数
     */
    public int numPz();

    /**
     * Jacobi coefficient
     *
     * @param x
     * @param y
     * @return
     */
    public double J(double x, double y);

    /**
     * Metric coefficient
     *
     * @param x
     * @param y
     * @return
     */
    public double alpha(double x, double y);

    /**
     * Metric coefficient
     *
     * @param x
     * @param y
     * @return
     */
    public double gamma(double x, double y);

    /**
     * 判断正交
     *
     * @param x
     * @param y
     * @return
     */
    public double beta(double x, double y);

    /**
     *
     * @param x1
     * @param y1
     * @return get Cartesian coordinates x
     */
    public double realX(double x1, double y1);

    /**
     *
     * @param x1
     * @param y1
     * @return get Cartesian coordinates y
     */
    public double realY(double x1, double y1);

    public double getVol(Label flag);

    public double getVol(int x, int y, int z);

    public Mesh clone();

    public Mesh fineMesh(Mesh cellf);
}
