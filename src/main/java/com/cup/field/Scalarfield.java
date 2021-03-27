/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.field;

import com.cup.util.Tool;

/**
 * 定义三维标量场
 *
 * @author winsway
 */
public class Scalarfield extends AbstractField {

    public Scalarfield(int numX, int numY, int numZ, String name) {
        super(numX, numY, numZ, name);
    }

    /**
     *
     * @param numX X方向节点数
     * @param numY Y方向节点数
     * @param numZ Z方向节点数
     */
    public Scalarfield(int numX, int numY, int numZ) {
        super(numX, numY, numZ);
    }

    /**
     * 二维场构造函数
     *
     * @param x
     * @param y
     */
    public Scalarfield(int x, int y) {
        this(x, y, 1);
        this.oldTime = new double[x][y][numZ];
        this.newTime = new double[x][y][numZ];
    }

    @Override
    public Scalarfield clone() {
        Scalarfield clone = new Scalarfield(numX, numY, numZ);
        Tool.copyAtoB(this.newTime, clone.newTime);
        Tool.copyAtoB(this.oldTime, clone.oldTime);
        Tool.copyAtoB(this.temp, clone.temp);
        clone.bound = this.bound.clone();

        return clone;
    }

    public double getOldField(int x, int y, int z) {
        return this.oldTime[x][y][z];
    }

    public double getNewField(int x, int y, int z) {
        return this.newTime[x][y][z];
    }

}
