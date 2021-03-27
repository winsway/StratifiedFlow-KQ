package com.cup.mesh.factory;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
import static java.lang.Math.cos;
import static java.lang.Math.cosh;
import static java.lang.Math.sin;
import static java.lang.Math.sinh;
import com.cup.mesh.AbstractMesh;
import static java.lang.Math.PI;

/**
 *
 * 更改DXU,DYV,DWZ的索引
 *
 * @author winsway
 */
public class BiSameV2 extends AbstractMesh {

    /**
     * under-relax factor alpha=0.95
     */
    public double radia;

    public double theta1, theta2, theta3;

    public double a;

    public int interhL;

    private void setBlockMeshX() {
//      网格点的位置
        tPx[0] = -6;
        for (int i = 1; i < tx.length; ++i) {
            tPx[i] = (tx[i] + tx[i - 1]) / 2;
        }
        tPx[numPx - 1] = 6;
        /**
         * 设定X方向网格线之间的距离
         */
        DXU = new double[NX + 2];
        for (int i = 1; i <= NX; ++i) {
            DXU[i] = tx[i] - tx[i - 1];
        }
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

    private void setBlockMeshY() {
        /**
         * 设置Y方向网格点P的坐标
         */
        tPy[0] = theta1;
        for (int i = 1; i < ty.length; ++i) {
            tPy[i] = (ty[i] + ty[i - 1]) / 2;
        }
        tPy[numPy - 1] = theta3;
        /**
         * 设定Y方向网格线之间的距离
         */
        DYV = new double[NY + 2];
        for (int i = 1; i <= NY; ++i) {
            DYV[i] = ty[i] - ty[i - 1];
        }
        /**
         * 设定Y方向P点之间的距离
         */
        DYP = new double[numPy - 1];
        for (int i = 0; i < DYP.length; ++i) {
            DYP[i] = tPy[i + 1] - tPy[i];
        }
        /**
         * Y方向设定完成
         */
    }

    /**
     * initial mesh
     */
    @Override
    public void blockMesh() {
        a = radia * Math.sin(theta1);
        theta3 = theta1 + PI;
//      1.0- hl 气相比例
        double ratio = (theta2 - theta1) / PI;
//        
        int ny = NY / 2;
        interhL = ny;
        this.getLiquidMesh();
        this.getGasMesh();
        this.setNonUniformX(-6, 12, 0, NX, tx);
        this.setBlockMeshX();
        this.setBlockMeshY();
        setDistance(numPz, tz, tPz, LZ, NZ);
    }

    private void getLiquidMesh() {
        setNonUniformMesh(theta2, (theta3 - theta2), interhL, NY, ty);
    }

    private void getGasMesh() {
        int ny = interhL;
        setNonUniformMesh(theta1, (theta2 - theta1), 0, ny, ty);

    }

    private void setNonUniformMesh(double init, double h, int xb, int xe, double[] y) {
        double ALG = 0.95;
        double xi;
        int M = xe - xb;
        for (int j = xb - xb; j <= xe - xb; j++) {
            xi = -1.0 + 2.0 * j / M;
            y[xb + j] = h / 2.0 / ALG
                    * Math.tanh(0.5 * xi
                            * Math.log((1 + ALG) / (1 - ALG)));
        }
        double dh = init - y[xb];
        for (int j = xb; j <= xe; j++) {
            y[j] = y[j] + dh;
            System.out.println("j = " + j + "; y = " + y[j]);
        }
    }

    private void setNonUniformX(double init, double h, int xb, int xe, double[] X) {
        double ALG = 0.95;
        double xi;
        int M = xe - xb;
        double[] x1 = new double[X.length];
        for (int j = xb - xb; j <= xe - xb; j++) {
            xi = -1.0 + 2.0 * j / M;
            x1[xb + j] = h / 2.0 / ALG
                    * Math.tanh(0.5 * xi
                            * Math.log((1 + ALG) / (1 - ALG)));
//            System.out.println("j = " + j + "; x1 = " + x1[j]);
        }
        double h2 = h / 2.0;
        for (int j = xb; j <= (xe - xb) / 2; j++) {
            x1[j] = h2 + x1[j];
//            System.out.println("j = " + j + "; x1 = " + x1[j]);
        }
        for (int j = 0; j <= (xe - xb) / 2; j++) {
            X[j] = -x1[(xe - xb) / 2 - j];
            X[j + (xe - xb) / 2] = x1[j];
        }
        for (int i = 0; i < X.length; i++) {
            System.out.println("I = " + i + "; X = " + X[i]);
        }
    }

    @Override
    public double realX(double x, double y) {
        double temp;
        temp = a * Math.sinh(x) / (Math.cosh(x) - Math.cos(y));
        return temp;
    }

    @Override
    public double realY(double x, double y) {
        double temp;
        temp = a * Math.sin(y) / (Math.cosh(x) - Math.cos(y));
        return temp;
    }

    public BiSameV2(double X, double Y, double Z, int nx, int ny, int nz) {
        super(X, Y, Z, nx, ny, nz);
    }

    public BiSameV2(double X, double Y, int nx, int ny) {
        this(X, Y, 1, nx, ny, 1);
    }

    @Override
    public double J(double x, double y) {
        double c = cosh(x) * cos(y) - 1.0;
        double s = sinh(x) * sin(y);
        return a * a / (c * c + s * s);
    }

    @Override
    public double alpha(double x, double y) {
        double c = cosh(x) * cos(y) - 1.0;
        double s = sinh(x) * sin(y);
        return a * a / (c * c + s * s);
    }

    @Override
    public double gamma(double x, double y) {
        double c = cosh(x) * cos(y) - 1.0;
        double s = sinh(x) * sin(y);
        return a * a / (c * c + s * s);
    }

    @Override
    public double beta(double x, double y) {
        return 0;
    }

    @Override
    public double getVol(int x, int y, int z) {
        return this.getDXU()[x] * this.getDYV()[y]
                * this.J(this.gettPx()[x], this.gettPy()[y]);
    }

    @Override
    public BiSameV2 clone() {
        BiSameV2 clone = new BiSameV2(this.LX, this.LY, this.LZ, this.NX, this.NY, this.NZ);
        clone.a = this.a;
        clone.interhL = this.interhL;
        clone.radia = this.radia;
        clone.theta1 = this.theta1;
        clone.theta2 = this.theta2;
        clone.theta3 = this.theta3;
        clone.blockMesh();
        return clone;
    }

}