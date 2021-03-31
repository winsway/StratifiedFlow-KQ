package com.winswe.mesh.factory;

import com.alibaba.fastjson.JSONObject;
import com.winswe.io.IOobject;
import com.winswe.mesh.AbstractMesh;
import com.winswe.mesh.Move;
import com.winswe.mesh.Structed2D;
import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.cos;
import static java.lang.Math.cosh;
import static java.lang.Math.sin;
import static java.lang.Math.sinh;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * 双极坐标网格
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年2月12日 下午3:13:12
 */
public class BipolarCoordianteMesh
        extends AbstractMesh
        implements
        Structed2D, Move {

    /**
     * 单个网格大小
     */
    private final double deltaX1;

    private final double deltaX2;

    /**
     * 网格线的位置
     */
    private final double[] lineX1;

    private final double[] lineX2;

    /**
     * 计算节点的位置
     */
    private final double[] pointX1;

    private final double[] pointX2;

    /**
     * X方向网格线距离； Distance between of two line on X direction
     */
    private double[] distanceLineX1;

    /**
     * Y方向网格线距离； Distance between of two line on Y direction
     */
    private double[] distanceLineX2;

    /**
     * P在X方向两节点之间的距离；
     */
    private double[] distancePointX1;

    /**
     * P在Y方向两节点之间的距离；
     */
    private double[] distancePointX2;

    /**
     * 管径，单位m
     */
    private final double diameter;

    /**
     *
     */
    private double theta1;

    /**
     *
     */
    private final double theta2 = PI;

    /**
     *
     */
    private double theta3;

    /**
     *
     */
    private double a;

    /**
     * 液位高度的网格数
     */
    private final int numberHLMesh;

    /**
     * X方向网格数目
     */
    private final int NX;

    /**
     * Y方向网格数目
     */
    private final int NY;

    /**
     * 几何参数
     */
    private final JSONObject geometricJson;

    public BipolarCoordianteMesh(IOobject ioObject) {
        geometricJson = ioObject.getJsonObject().getJSONObject("geometric");
        //设定几何参数        
        diameter = geometricJson.getDoubleValue("diameter");
        theta1 = acos(1.0 - 2 * 0.5);
        theta3 = theta1 + PI;
        a = diameter / 2 * Math.sin(theta1);
        //设定网格数目
        NX = geometricJson.getIntValue("NX");
        NY = geometricJson.getIntValue("NY");
        numberHLMesh = NY / 2;
        //网格位置 
        this.lineX1 = new double[NX + 1];
        this.lineX2 = new double[NY + 1];
        this.pointX1 = new double[NX + 2];
        this.pointX2 = new double[NY + 2];
        deltaX1 = 12 / NX;
        deltaX2 = PI / NY;
    }

    /**
     * 雅可比因子
     *
     * @param x x坐标
     * @param y y坐标
     * @return 雅可比因子大小
     */
    public double Jacobi(double x, double y) {
        double c = cosh(x) * cos(y) - 1.0;
        double s = sinh(x) * sin(y);
        return a * a / (c * c + s * s);
    }

    /**
     * eta方向的度规系数，放映了eta方向网格的疏密程度
     *
     * @param x x坐标
     * @param y y坐标
     * @return eta方向的度规系数
     */
    @Override
    public double alpha(double x, double y) {
        double c = cosh(x) * cos(y) - 1.0;
        double s = sinh(x) * sin(y);
        return a * a / (c * c + s * s);
    }

    /**
     * xi方向的度规系数，放映xi方向网格线的疏密程度
     *
     * @param x x坐标
     * @param y y坐标
     * @return xi方向的度规系数
     */
    public double gamma(double x, double y) {
        double c = cosh(x) * cos(y) - 1.0;
        double s = sinh(x) * sin(y);
        return a * a / (c * c + s * s);
    }

    /**
     * 反映了物理平面上网格的正交性，网格局部正交时beta=0
     *
     * @param x x坐标
     * @param y y坐标
     * @return 反映了物理平面上网格的正交性，网格局部正交时beta=0
     */
    public double beta(double x, double y) {
        return 0;
    }

    /**
     *
     * @param indexX X方向位置索引
     * @param indexY Y方向位置索引
     * @return 体积大小
     */
    @Override
    public double getVolume(int indexX, int indexY) {
        return distanceLineX1[indexX]
                * distanceLineX2[indexY]
                * Jacobi(pointX1[indexX], pointX2[indexY]);
    }

    @Override
    public void blockMesh() {
        this.setWaterMesh();
        this.setOilMesh();
        this.setNonUniformX(12, 0, NX, lineX1);
        this.setBlockMeshX();
        this.setBlockMeshY();
    }

    /**
     * 设定水相网格
     */
    private void setWaterMesh() {
        setNonUniformMesh(theta2, (theta3 - theta2), numberHLMesh, NY, lineX2);
    }

    /**
     * 设定油相网格
     */
    private void setOilMesh() {
        setNonUniformMesh(theta1, (theta2 - theta1), 0, numberHLMesh, lineX2);
    }

    /**
     * 设置Y方向的网格线
     *
     * @param init 初始值
     * @param h 中间位置
     * @param xb 起点网格数
     * @param xe 终点网格数
     * @param Y y方向的网格线
     */
    private void setNonUniformMesh(double init, double h, int xb, int xe, double[] Y) {
        double ALG = 0.95;
        double xi;
        int M = xe - xb;
        for (int j = xb - xb; j <= xe - xb; j++) {
            xi = -1.0 + 2.0 * j / M;
            Y[xb + j] = h / 2.0 / ALG
                    * Math.tanh(0.5 * xi
                            * Math.log((1 + ALG) / (1 - ALG)));
        }
        double dh = init - Y[xb];
        for (int j = xb; j <= xe; j++) {
            Y[j] = Y[j] + dh;
//            System.out.println("j = " + j + "; y = " + Y[j]);
        }
    }

    /**
     * 设置X方向的网格线
     *
     * @param h X方向宽6*2=12
     * @param xb X方向起点网格数
     * @param xe X方向终点网格数
     * @param X X方向网格线的位置
     */
    private void setNonUniformX(double h, int xb, int xe, double[] X) {
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
//            System.out.println("I = " + i + "; X = " + X[i]);
//        }
    }

    /**
     * 设定X方向的坐标
     */
    private void setBlockMeshX() {

        //网格点的位置
        pointX1[0] = -6;
        for (int i = 1; i < lineX1.length; ++i) {
            pointX1[i] = (lineX1[i] + lineX1[i - 1]) / 2;
        }
        pointX1[NX + 1] = 6;

        //设定X方向网格线之间的距离
        distanceLineX1 = new double[NX + 2];
        for (int i = 1; i <= NX; ++i) {
            distanceLineX1[i] = lineX1[i] - lineX1[i - 1];
        }

        //设定X方向P点之间的距离
        distancePointX1 = new double[NX + 1];
        for (int i = 0; i < distancePointX1.length; ++i) {
            distancePointX1[i] = pointX1[i + 1] - pointX1[i];
        }

        //X方向设定完成
    }

    /**
     * 设定Y方向的坐标
     */
    private void setBlockMeshY() {

        //设置Y方向网格点P的坐标
        pointX2[0] = theta1;
        for (int i = 1; i < lineX2.length; ++i) {
            pointX2[i] = (lineX2[i] + lineX2[i - 1]) / 2;
        }
        pointX2[NY + 1] = theta3;

        //设定Y方向网格线之间的距离
        distanceLineX2 = new double[NY + 2];
        for (int i = 1; i <= NY; ++i) {
            distanceLineX2[i] = lineX2[i] - lineX2[i - 1];
        }

        // 设定Y方向P点之间的距离
        distancePointX2 = new double[NY + 1];
        for (int i = 0; i < distancePointX2.length; ++i) {
            distancePointX2[i] = pointX2[i + 1] - pointX2[i];
        }

        //Y方向设定完成
    }

    @Override
    public void move(double hl) {

        theta1 = acos(1.0 - 2 * hl);
        theta3 = theta1 + PI;
        a = diameter / 2 * Math.sin(theta1);
        this.blockMesh();
    }

    @Override
    public int getNX() {
        return NX;
    }

    @Override
    public int getNY() {
        return NY;
    }

    @Override
    public double X(double x1, double x2) {
        return a * Math.sinh(x1) / (Math.cosh(x1) - Math.cos(x2));
    }

    @Override
    public double Y(double x1, double x2) {
        return a * Math.sin(x2) / (Math.cosh(x1) - Math.cos(x2));
    }

    @Override
    public double[] getPointX1() {
        return pointX1;
    }

    @Override
    public double[] getPointX2() {
        return pointX2;
    }

    @Override
    public double[] getDYV() {
        return distanceLineX2;
    }

    @Override
    public double[] getDXP() {
        return distancePointX1;
    }

    @Override
    public double[] getDXU() {
        return distanceLineX1;
    }

    @Override
    public double[] getDYP() {
        return distancePointX2;
    }

    @Override
    public double[] getLineX1() {
        return lineX1;
    }

    @Override
    public double[] getLineX2() {
        return lineX2;
    }

}
