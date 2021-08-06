package com.winswe.mesh;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * structure mesh
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年2月12日 下午3:11:21
 */
public interface Structed2D extends Move {

    /**
     * return mesh index
     *
     * @param I x index
     * @param J y index
     * @return index
     */
    default int getCellIndex(int I, int J) {
        //two side with ghost cell
        return (I) + (J) * (this.getNX() + 1);
    }

    /**
     * mesh discrete
     */
    public void blockMesh();

    /**
     *
     * @return number of X direction cell
     */
    public int getNX();

    /**
     *
     * @return number of Y direction cell
     */
    public int getNY();

    /**
     * physical coordinate X
     *
     * @param x1 computing coordinate x1
     * @param x2 computing coordinate x2
     * @return physical coordinate X
     */
    public double X(double x1, double x2);

    /**
     *
     *
     * @param x1 computing coordinate x1
     * @param x2 computing coordinate x2
     * @return physical coordinate Y
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

    /**
     *
     * @return
     */
    public double[] getDYV();

    /**
     *
     * @return
     */
    public double[] getDXP();

    /**
     *
     * @return
     */
    public double[] getDXU();

    /**
     *
     * @return
     */
    public double[] getDYP();

    /**
     *
     * @param x1
     * @param x2
     * @return
     */
    public double alpha(double x1, double x2);

}
