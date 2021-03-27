/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.field;

import com.cup.util.Tool;

/**
 *
 * @author winsw
 */
public class ProF {

    double[][][] F;
    int x, y, z;

    public ProF(int x, int y, int z) {
        this.x = x;
        this.y = y;
        this.z = z;
        F = new double[x][y][z];
    }

    public double[][][] getF() {
        return F;
    }

    @Override
    public ProF clone() throws CloneNotSupportedException {
        ProF clone = new ProF(x, y, z);
        Tool.copyAtoB(this.F, clone.F);
        return clone;
    }

}
