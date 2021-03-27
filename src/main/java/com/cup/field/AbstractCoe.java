/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.field;

/**
 *
 * @author winsway
 */
public abstract class AbstractCoe implements CoeffMatrix {

    protected double ap, aw, ae, as, an, ab, at;
    protected double b;

    protected String name;

    @Override
    public abstract void setB(int x, int y, int z, double value);

    @Override
    public abstract void setAp(int x, int y, int z, double value);

    @Override
    public abstract void setAt(int x, int y, int z, double value);

    @Override
    public abstract void setAb(int x, int y, int z, double value);

    @Override
    public abstract void setAn(int x, int y, int z, double value);

    @Override
    public abstract void setAs(int x, int y, int z, double value);

    @Override
    public abstract void setAe(int x, int y, int z, double value);

    @Override
    public abstract void setAw(int x, int y, int z, double value);

    /**
     * 求和周围系数
     */
    double sumAnb() {
        return (aw + ae + as + an + ab + at);
    }

}
