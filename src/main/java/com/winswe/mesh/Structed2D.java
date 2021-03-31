package com.winswe.mesh;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * 结构网格
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年2月12日 下午3:11:21
 */
public interface Structed2D extends Move {

    /**
     * 返回网格索引
     *
     * @param I 横坐表
     * @param J 纵坐标
     * @return 索引值
     */
    default int getCellIndex(int I, int J) {
        //两边有虚网格
        return (I) + (J) * (this.getNX() + 1);
    }

    /**
     * 网格离散
     */
    public void blockMesh();

    /**
     *
     * @return
     */
    public int getNX();

    /**
     *
     * @return
     */
    public int getNY();

    /**
     *
     * @param x1 计算坐标x1
     * @param x2 计算坐标x2
     * @return 物理坐标系X
     */
    public double X(double x1, double x2);

    /**
     *
     *
     * @param x1 计算坐标x1
     * @param x2 计算坐标x2
     * @return 物理坐标系Y
     */
    public double Y(double x1, double x2);

    /**
     *
     * @return
     */
    public double[] getLineX1();

    /**
     *
     * @return
     */
    public double[] getLineX2();

    /**
     *
     * @return
     */
    public double[] getPointX1();

    /**
     *
     * @return
     */
    public double[] getPointX2();

    /**
     *
     * @param X
     * @param Y
     * @return
     */
    public double getVolume(int X, int Y);

    public double[] getDYV();

    public double[] getDXP();

    public double[] getDXU();

    public double[] getDYP();

    public double alpha(double x1, double x2);

}
