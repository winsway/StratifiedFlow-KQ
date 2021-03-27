package com.cup.mesh.factory;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
import com.cup.mesh.AbstractMesh;
import com.cup.mesh.Mesh;
import com.cup.turbulence.DW;
import com.winswe.finiteVolume.geom.Point;
import static java.lang.Math.PI;

/**
 *
 * 1.需要加入坐标数组
 *
 * @author winsway
 */
public class CavityMesh2 extends AbstractMesh {

    /**
     * under-relax factor alpha=0.95
     */
    public double radia;

    public double theta1, theta2, theta3;

    public double a;

    public int interhL;

    /**
     * node从1开始,存储真实物理坐标
     */
    public Point[][][] node, cv;
    /**
     * 和cv对应
     */
    public DW[][][] DN;

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
//        
        int ny = NY / 2;
        interhL = ny;
        this.getLiquidMesh();
        this.getGasMesh();
        this.setNonUniformX(-6, 12, 0, NX, tx);
        this.setBlockMeshX();
        this.setBlockMeshY();
        setDistance(numPz, tz, tPz, LZ, NZ);
//      
        createNode();
        createCV();
        createDN();
    }

    /**
     * 创建节点同时赋值
     */
    private void createNode() {
        node = new Point[NX + 1 + 2][NY + 1 + 2][NZ + 1 + 2];
        double x, y, z;
        for (int K = 1; K <= NZ + 1; K++) {
            for (int J = 1; J <= NY + 1; J++) {
                for (int I = 1; I <= NX + 1; I++) {
                    x = realX(this.getTx()[I - 1], this.getTy()[J - 1]);
                    y = realY(this.getTx()[I - 1], this.getTy()[J - 1]);
                    z = this.getTz()[K - 1];
                    node[I][J][K] = new Point(x, y, z);
                }
            }
        }

    }

    /**
     * 得到网格中心值
     */
    private void createCV() {
        cv = new Point[NX + 2][NY + 2][NZ + 2];
        double x, y, z;
        for (int K = 0; K < numPz; K++) {
            for (int J = 0; J < numPy; J++) {
                for (int I = 0; I < numPx; I++) {
                    x = realX(this.gettPx()[I], this.gettPy()[J]);
                    y = realY(this.gettPx()[I], this.gettPy()[J]);
                    z = this.gettPz()[K];
                    cv[I][J][K] = new Point(x, y, z);
                }
            }
        }
    }

    /**
     * 得到不同网格中心对应的坐标与距离
     */
    private void createDN() {
        DN = new DW[NX + 2][NY + 2][NZ + 2];
        Mesh mesh = this;
        for (int K = 0; K < numPz; K++) {
            for (int J = 0; J < numPy; J++) {
                for (int I = 0; I < numPx; I++) {
                    DN[I][J][K] = new DW();
                    DN[I][J][K].cacD(mesh, cv[I][J][K]);
                }
            }
        }
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
//            System.out.println("j = " + j + "; y = " + y[j]);
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
//        for (int i = 0; i < X.length; i++) {
//           System.out.println("I = " + i + "; X = " + X[i]);
//        }
    }

    @Override
    public double realX(double x, double y) {
        return x;
    }

    @Override
    public double realY(double x, double y) {
        return y;
    }

    public CavityMesh2(double X, double Y, double Z, int nx, int ny, int nz) {
        super(X, Y, Z, nx, ny, nz);
    }

    public CavityMesh2(double X, double Y, int nx, int ny) {
        this(X, Y, 1, nx, ny, 1);
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
    public double getVol(int x, int y, int z) {
        return this.getDXU()[x] * this.getDYV()[y]
                * this.J(this.gettPx()[x], this.gettPy()[y]);
    }

    @Override
    public CavityMesh2 clone() {
        CavityMesh2 clone = new CavityMesh2(this.LX, this.LY, this.LZ, this.NX, this.NY, this.NZ);
        clone.a = this.a;
        clone.interhL = this.interhL;
        clone.radia = this.radia;
        clone.theta1 = this.theta1;
        clone.theta2 = this.theta2;
        clone.theta3 = this.theta3;
        clone.blockMesh();
        return clone;
    }

    @Override
    public Mesh fineMesh(Mesh cellf) {
        CavityMesh2 clone = new CavityMesh2(this.LX, this.LY, this.LZ, this.NX / 2, this.NY / 2, this.NZ);
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
