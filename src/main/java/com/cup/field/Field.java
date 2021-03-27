/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.field;

import com.cup.boundary.Boundary;
import com.cup.field.*;
import com.cup.log.Infor;

/**
 *
 * @author winsway
 */
public interface Field {

    /**
     * Number of rows in the matrix
     *
     * @return
     */
    int numX();

    /**
     * Number of columns in the matrix
     *
     * @return
     */
    int numY();

    /**
     * Number of deep in the matrix
     *
     * @return
     */
    int numZ();

    /**
     * 获取对象场名
     *
     * @return
     */
    public String getName();

    public double[][][] getNewField();

    public double[][][] getOldField();

    public double[][][] getOoldField();

    public double[][][] getTempField();

    public void setNewField(int x, int y, int z, double value);

    /**
     * 新值赋给旧值
     */
    public void copyNew2Old();

    public void copyAtoB(double[][][] A, double[][][] B);

    public Boundary getBoundary();

    public Infor getLog();

    public Field clone();

}
