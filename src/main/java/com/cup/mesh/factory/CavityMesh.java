/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.mesh.factory;

import com.cup.boundary.factory.Label;
import com.cup.mesh.AbstractMesh;
import com.cup.mesh.Mesh;

/**
 * 考虑虚拟网格的实现
 *
 * @author winsway
 */
public class CavityMesh extends AbstractMesh {

    public CavityMesh(double X, double Y, double Z, int nx, int ny, int nz) {
        super(X, Y, Z, nx, ny, nz);
    }

    public CavityMesh(double X, double Y, int nx, int ny) {
        this(X, Y, 1, nx, ny, 1);
    }

    @Override
    public Mesh fineMesh(Mesh fine) {
        int nx = fine.getNX();
        int ny = fine.getNY();
        int nz = fine.getNZ();
        Mesh cell = new CavityMesh(this.LX, this.LY, 1, nx / 2, ny / 2, 1);
        cell.blockMesh();
        return cell;
    }

    @Override
    public double getVol(Label flag) {
        int x, y, z;
        x = flag.i;
        y = flag.j;
        z = flag.k;
        return DXU[x] * DYV[y] * J(tPx[x], tPy[y]);
    }

    @Override
    public void blockMesh() {
        //        先定点后定线是FDM的思路
        setBlockMeshX();
        setBlockMeshY();
        setBlockMeshZ();
    }

    public void settPx() {
        var ddx = LX / NX;
        tPx[0] = 0;
        for (int i = 1; i <= NX; i++) {
            tPx[i] = ddx * (i - 0.5);
        }
        tPx[NX + 1] = LX;
    }

    public void setTx() {
        var ddx = LX / NX;
        tx[0] = 0;
        for (int i = 1; i <= NX; i++) {
            tx[i] = ddx * (i) * 1.0;
        }
        tx[tx.length - 1] = LX;
    }

    private void setBlockMeshX() {
        settPx();
        setTx();
        /**
         * 设定X方向网格线之间的距离
         */
        DXU = new double[NX + 2];
        for (int i = 1; i <= NX; ++i) {
            DXU[i] = tx[i] - tx[i - 1];
        }
        DXU[0] = 0;
        DXU[NX + 1] = 0;
        /**
         * 设定X方向P点之间的距离
         */
        DXP = new double[numPx - 1];
        for (int i = 0; i < DXP.length; ++i) {
            DXP[i] = tPx[i + 1] - tPx[i];
        }
        /**
         * X方向设定完成
         */
    }

    public void settPy() {
        var ddy = LY / NY;
        tPy[0] = 0;
        for (int j = 1; j <= NY; j++) {
            tPy[j] = ddy * (j - 0.5);
        }
        tPy[NY + 1] = LY;
    }

    public void setTy() {
        var ddy = LY / NY;
        ty[0] = 0;
        for (int i = 1; i < ty.length - 1; ++i) {
            ty[i] = ddy * i * 1.0;
        }
        ty[ty.length - 1] = LY;
    }

    private void setBlockMeshY() {
        settPy();
        setTy();
        /**
         * 设定Y方向网格线之间的距离
         */
        DYV = new double[NY + 2];
        for (int i = 1; i <= NY; ++i) {
            DYV[i] = ty[i] - ty[i - 1];
        }
        DYV[0] = 0;
        DYV[NY + 1] = 0;
        /**
         * 设定Y方向P点之间的距离
         */
        DYP = new double[numPy - 1];
        for (int i = 0; i < DYP.length; ++i) {
            DYP[i] = tPy[i + 1] - tPy[i];
        }
    }

    private void setBlockMeshZ() {
        setDistance(numPz, tz, tPz, LZ, NZ);
        /**
         * 设定Z方向网格线之间的距离
         */
        DZW = new double[NZ + 2];
        for (int i = 1; i <= NZ; ++i) {
            DZW[i] = tz[i] - tz[i - 1];
        }
        DZW[0] = 0;
        DZW[NZ + 1] = 0;
        /**
         * 设定Z方向P点之间的距离
         */
        DZP = new double[numPz - 1];
        for (int i = 0; i < DZP.length; ++i) {
            DZP[i] = tPz[i + 1] - tPz[i];
        }
        /**
         * Z方向设定完成
         */
    }

    @Override
    public double getVol(int x, int y, int z) {
        return DXU[x] * DYV[y];
    }

    @Override
    public Mesh clone() {
        CavityMesh clone = new CavityMesh(this.LX, this.LY, this.LZ, this.NX, this.NY, this.NZ);
        clone.blockMesh();
        return clone;
    }

    @Override
    public double J(double x, double y) {
        return 1.0;
    }

    @Override
    public double alpha(double x, double y) {
        return 1.0;
    }

    @Override
    public double gamma(double x, double y) {
        return 1.0;
    }

    @Override
    public double beta(double x, double y) {
        return 1.0;
    }

    @Override
    public double realX(double x1, double y1) {
        return x1;
    }

    @Override
    public double realY(double x1, double y1) {
        return y1;
    }
}
