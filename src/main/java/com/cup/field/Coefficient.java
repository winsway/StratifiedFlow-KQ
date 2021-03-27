/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.field;

import com.cup.boundary.factory.Label;
import com.cup.system.SystemControl;
import com.cup.util.Tool;

/**
 * 构造系数矩阵
 *
 * @author winsway
 */
public class Coefficient extends AbstractCoe {

    public double[][][] ap, aw, ae, as, an, ab, at;
    public double[][][] b;
    double ddt;
    double schemePhi = 0;
    public double Pressure = 0, Other = 0;
    int x, y, z;

    /**
     * 设定非稳态项
     *
     * @param rho 密度
     * @param vol 体积
     * @param dt 时间dt
     * @param var 控制变量“Phi，U，V”
     */
    void setDdt(double rho, double vol, double dt, String var) {
        if (dt == 0 || "Phi".equals(var)) {
            ddt = 0;
        } else {
            ddt = rho * vol / dt;
        }
    }

    /**
     * 构造函数
     *
     * @param x x点的数量
     * @param y y点的数量
     */
    public Coefficient(int x, int y) {
        this(x, y, 3);
    }

    /**
     * 构造函数
     *
     * @param x x点的数量
     * @param y y点的数量
     * @param z z点的数量
     */
    public Coefficient(int x, int y, int z) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.ap = new double[x][y][z];
        this.aw = new double[x][y][z];
        this.ae = new double[x][y][z];
        this.as = new double[x][y][z];
        this.an = new double[x][y][z];
        this.ab = new double[x][y][z];
        this.at = new double[x][y][z];
        this.b = new double[x][y][z];
    }

    @Override
    public void setB(int x, int y, int z, double value) {
        b[x][y][z] = value;
    }

    @Override
    public void setAp(int x, int y, int z, double value) {
        ap[x][y][z] = value;
    }

    @Override
    public void setAt(int x, int y, int z, double value) {
        at[x][y][z] = value;
    }

    @Override
    public void setAb(int x, int y, int z, double value) {
        ab[x][y][z] = value;
    }

    @Override
    public void setAn(int x, int y, int z, double value) {
        an[x][y][z] = value;
    }

    @Override
    public void setAs(int x, int y, int z, double value) {
        as[x][y][z] = value;
    }

    @Override
    public void setAe(int x, int y, int z, double value) {
        ae[x][y][z] = value;
    }

    @Override
    public void setAw(int x, int y, int z, double value) {
        aw[x][y][z] = value;
    }

    @Override
    public void addB(Field field, double[][][] phi, SystemControl sys) {
        Label flag = new Label(field.getName());

    }

    @Override
    public void addB(double[][][] phi) {
        for (int k = 1; k < b[0][0].length - 1; ++k) {
            for (int j = 1; j < b[0].length - 1; ++j) {
                for (int i = 1; i < b.length - 1; ++i) {
                    b[i][j][k] = b[i][j][k] + phi[i][j][k];
                }
            }
        }
    }

    @Override
    public double getAw(int x, int y, int z) {
        return aw[x][y][z];
    }

    @Override
    public double getAe(int x, int y, int z) {
        return ae[x][y][z];
    }

    @Override
    public double getAs(int x, int y, int z) {
        return as[x][y][z];
    }

    @Override
    public double getAn(int x, int y, int z) {
        return an[x][y][z];
    }

    @Override
    public double getAb(int x, int y, int z) {
        return ab[x][y][z];
    }

    @Override
    public double getAt(int x, int y, int z) {
        return at[x][y][z];
    }

    @Override
    public double getAp(int x, int y, int z) {
        return ap[x][y][z];
    }

    @Override
    public double getB(int x, int y, int z) {
        return b[x][y][z];
    }

    @Override
    public Coefficient clone() throws CloneNotSupportedException {
        Coefficient clone = new Coefficient(x, y, z);
        Tool.copyAtoB(this.aw, clone.aw);
        Tool.copyAtoB(this.ae, clone.ae);
        Tool.copyAtoB(this.as, clone.as);
        Tool.copyAtoB(this.an, clone.an);
        Tool.copyAtoB(this.ab, clone.ab);
        Tool.copyAtoB(this.at, clone.at);
        Tool.copyAtoB(this.ap, clone.ap);
        Tool.copyAtoB(this.b, clone.b);
        return clone;
    }

}
