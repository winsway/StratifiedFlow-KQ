/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.field;

import com.cup.boundary.factory.BoundaryRegion;
import com.cup.boundary.Boundary;
import com.cup.system.FvSolution;
import com.cup.log.Infor;

/**
 *
 * @author winsway
 */
public abstract class AbstractField implements Field {

    /**
     * 场名字
     */
    protected String name;

    /**
     * 边界条件设定
     */
    protected BoundaryRegion bound;

    /**
     * 上上一时步场值
     */
    protected double[][][] ooldTime;

    /**
     * 上一时步场值
     */
    protected double[][][] oldTime;

    /**
     * 当前时步场值
     */
    protected double[][][] newTime;

    /**
     * 求解过程的中间变量场
     */
    protected double[][][] temp;

    /**
     * Number of rows
     */
    protected int numX;

    /**
     * Number of columns
     */
    protected int numY;

    /**
     * Number of deep
     */
    protected int numZ;

    public FvSolution fvsolution;
    public Infor infor = new Infor();

    /**
     * AbstractField构造函数
     *
     * @param numX X方向节点数
     * @param numY Y方向节点数
     * @param numZ Z方向节点数
     * @param name 场名
     */
    protected AbstractField(int numX, int numY, int numZ, String name) {
        this.numX = numX;
        this.numY = numY;
        this.numZ = numZ;
        this.name = name;
        bound = new BoundaryRegion(numX, numY, numZ);
        newTime = new double[numX][numY][numZ];
        oldTime = new double[numX][numY][numZ];
        temp = new double[numX][numY][numZ];

    }

    /**
     * AbstractField构造函数
     *
     * @param numX X方向节点数
     * @param numY Y方向节点数
     * @param numZ Z方向节点数
     */
    protected AbstractField(int numX, int numY, int numZ) {
        this.numX = numX;
        this.numY = numY;
        this.numZ = numZ;
        bound = new BoundaryRegion(numX, numY, numZ);
        newTime = new double[numX][numY][numZ];
        oldTime = new double[numX][numY][numZ];
        temp = new double[numX][numY][numZ];
    }

    /**
     * Constructor for AbstractField, same size as A.The invoking constructor
     * should set this matrix equal the argument matrix
     *
     * @param A
     */
    protected AbstractField(Field A) {
        this(A.numX(), A.numY(), A.numZ());
    }

    @Override
    public int numX() {
        return numX;
    }

    @Override
    public int numY() {
        return numY;
    }

    @Override
    public int numZ() {
        return numZ;
    }

    public boolean is2D() {
        return numZ == 1;
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public void copyNew2Old() {
        this.copyAtoB(newTime, oldTime);
    }

    @Override
    public double[][][] getTempField() {
        return temp;
    }

    @Override
    public double[][][] getOoldField() {
        return ooldTime;
    }

    @Override
    public double[][][] getOldField() {
        return oldTime;
    }

    @Override
    public double[][][] getNewField() {
        return newTime;
    }

    @Override
    public void setNewField(int x, int y, int z, double value) {
        newTime[x][y][z] = value;
    }

    public void setBoundary(int x, int y, int z, int type) {
        bound.setType(x, y, z, type);
    }

    /**
     *
     * @return
     */
    @Override
    public Boundary getBoundary() {
        return bound;
    }

    @Override
    public Infor getLog() {
        return infor;
    }

    @Override
    public void copyAtoB(double[][][] A, double[][][] B) {
        for (int z = 0; z < this.numZ; ++z) {
            for (int j = 0; j < this.numY; ++j) {
                for (int i = 0; i < this.numX; ++i) {
                    B[i][j][z] = A[i][j][z];
                }
            }
        }
    }

    @Override
    public abstract Field clone();

}
