/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.boundary.factory;

import com.cup.boundary.Boundary;

/**
 *
 *
 * @author winsway
 */
public abstract class AbstractBoundary implements Boundary {

    /**
     * 边界条件X方向的点数
     */
    protected int numX;

    /**
     * 边界条件Y方向的点数
     */
    protected int numY;

    /**
     * 边界条件Z方向的点数
     */
    protected int numZ;

    @Override
    public abstract int judgeBoundary(int m);

    @Override
    public abstract int judgeBoundary(int x, int y, int z);

    @Override
    public abstract String getBoundaryTypes();

    @Override
    public abstract int getType(int x, int y, int z);

    protected AbstractBoundary(int X, int Y, int Z) {
        if (X < 0 || Y < 0 || Z < 0) {
            throw new IndexOutOfBoundsException("Matrix size cannot be negative");
        }
        this.numX = X;
        this.numY = Y;
        this.numZ = Z;
    }

}
