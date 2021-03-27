/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.field;

import com.cup.system.SystemControl;

/**
 * 系数矩阵接口
 *
 * @author winsway
 */
public interface CoeffMatrix {

    public void setAw(int x, int y, int z, double value);

    public void setAe(int x, int y, int z, double value);

    public void setAs(int x, int y, int z, double value);

    public void setAn(int x, int y, int z, double value);

    public void setAb(int x, int y, int z, double value);

    public void setAt(int x, int y, int z, double value);

    public void setAp(int x, int y, int z, double value);

    public void setB(int x, int y, int z, double value);

    public void addB(Field field, double[][][] phi, SystemControl sys);

    public void addB(double[][][] phi);

    public double getAw(int x, int y, int z);

    public double getAe(int x, int y, int z);

    public double getAs(int x, int y, int z);

    public double getAn(int x, int y, int z);

    public double getAb(int x, int y, int z);

    public double getAt(int x, int y, int z);

    public double getAp(int x, int y, int z);

    public double getB(int x, int y, int z);

}
