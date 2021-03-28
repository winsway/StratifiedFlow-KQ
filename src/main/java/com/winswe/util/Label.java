/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.util;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月6日 下午12:36:14
 */
public class Label {

    private int i, j;
    private int west, east, south, north;

    private final int NX, NY;

    public Label(int NX, int NY) {
        this.NX = NX;
        this.NY = NY;
    }

    public void setFlag(int X, int Y) {
        i = X;
        j = Y;
        west = (X == 1) ? 1 : 0;
        east = (X == NX) ? 1 : 0;
        south = (Y == 1) ? 1 : 0;
        north = (Y == NY) ? 1 : 0;
    }

    public int getI() {
        return i;
    }

    public int getJ() {
        return j;
    }

    public int getWest() {
        return west;
    }

    public int getEast() {
        return east;
    }

    public int getSouth() {
        return south;
    }

    public int getNorth() {
        return north;
    }

    /**
     *
     * @return true:at boundary.
     */
    public boolean atBoundary() {
        return (getWest() == 1)
                || (getEast() == 1)
                || (getSouth() == 1)
                || (getNorth() == 1);
    }
}
