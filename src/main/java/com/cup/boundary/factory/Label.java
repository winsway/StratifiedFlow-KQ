/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.boundary.factory;

import com.cup.boundary.Boundary;

/**
 *
 * @author winsway
 */
public class Label {

    public int w, e, s, n, t, b;
    public int i, j, k;
    public String var;

    public Label(String var) {
        this.var = var;
    }

    public void setFlag(int x, int y, int z) {
        this.i = x;
        this.j = y;
        this.k = z;
    }

    public void setFlag(Boundary phi, int i, int j, int k) {
        w = judgeBoundary(phi.getType(i - 1, j, k));
        e = judgeBoundary(phi.getType(i + 1, j, k));
        s = judgeBoundary(phi.getType(i, j - 1, k));
        n = judgeBoundary(phi.getType(i, j + 1, k));
        b = judgeBoundary(phi.getType(i, j, k - 1));
        t = judgeBoundary(phi.getType(i, j, k + 1));
        this.i = i;
        this.j = j;
        this.k = k;
    }

    public void setFlag(Boundary phi, int i, int j) {
        w = judgeBoundary(phi.getType(i - 1, j, 1));
        e = judgeBoundary(phi.getType(i + 1, j, 1));
        s = judgeBoundary(phi.getType(i, j - 1, 1));
        n = judgeBoundary(phi.getType(i, j + 1, 1));
        this.i = i;
        this.j = j;
        k = 1;
    }

    private int judgeBoundary(int m) {
        return (m == 0) ? 0 : 1;
    }

}
