/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.boundary;

import com.cup.boundary.factory.BoundaryRegion;

/**
 * 基本的边界条件接口。
 *
 * @author winsway
 */
public interface Boundary {

    /**
     * 获取边界条件类型
     *
     * @param x
     * @param y
     * @param z
     * @return
     */
    int getType(int x, int y, int z);

    /**
     * 设定边界条件类型
     *
     * @param x
     * @param y
     * @param z
     * @param value
     */
    void setType(int x, int y, int z, int value);

    /**
     * 得到边界条件类型
     *
     * @return 边界条件类型
     */
    String getBoundaryTypes();

    int judgeBoundary(int x, int y, int z);

    int judgeBoundary(int m);

    enum Type {
        /**
         * 第一类边界条件
         */
        fixedValue,
        /**
         * 零梯度
         */
        zeroGradient,
        /**
         * 空边界，主要是针对二维使用；
         */
        empty;

    }

    /**
     *
     * @param Type 边界条件类型
     * @param Aw 西边界面积
     * @param Vol 控制容积体积
     * @param dxw 到控制界面外点的距离
     * @param gammaw 扩散系数
     * @param Fw 西边界通量
     * @param hfw 换热系数
     * @return 返回附加源项系数部分
     */
    static public double Spad1(int Type, double Aw, double Vol, double dxw,
            double gammaw, double Fw, double hfw) {
        int A1 = 0, A2 = 0, A3 = 0;
        double B1, B2, B3;
        double spad;
        if (Type == 1) {
            A1 = 1;
        }
        if (Type == 2) {
            A2 = 1;
        }
        if (Type == 3) {
            A3 = 1;
        }
        B1 = -(Aw / (dxw / gammaw) + Math.max(Fw, 0)) / Vol;
        B2 = (0 + 0) / Vol;
        B3 = -Aw / Vol * 1 / ((dxw / gammaw) + (1 / hfw));//这个需要进一步处理
        spad = A1 * B1 + A2 * B2 + A3 * B3;
        return spad;
    }

    /**
     *
     * @param Type 类型
     * @param Ae 东边界面积
     * @param Vol 控制容积体积
     * @param dxe 到控制界面外点的距离
     * @param gammae 扩散系数
     * @param Fe 东边界通量
     * @param hfe 换热系数
     * @return 返回附加源项系数部分
     */
    static public double Spad2(int Type, double Ae, double Vol, double dxe,
            double gammae, double Fe, double hfe) {
        int A1 = 0, A2 = 0, A3 = 0;
        double B1, B2, B3;
        double spad;
        if (Type == 1) {
            A1 = 1;
        }
        if (Type == 2) {
            A2 = 1;
        }
        if (Type == 3) {
            A3 = 1;
        }
        B1 = -(Ae / (dxe / gammae) + Math.max(-Fe, 0)) / Vol;
        B2 = (0 + 0) / Vol;
        B3 = -Ae / Vol * 1 / ((dxe / gammae) + (1 / hfe));//需要进一步处理
        spad = A1 * B1 + A2 * B2 + A3 * B3;
        return spad;
    }

    /**
     *
     * @param Type 类型
     * @param Aw 西边界面积
     * @param Vol 控制容积体积
     * @param dxw 到控制界面外点的距离
     * @param gammaw 扩散系数
     * @param Fw 西边界通量
     * @param Tw 西边界定值
     * @param qw 热流密度
     * @param Tfw 西边界环境温度
     * @param hfw 换热系数
     * @return 返回附加源项常数部分
     */
    static public double Scad1(int Type, double Aw, double Vol, double dxw,
            double gammaw, double Fw, double Tw, double qw, double Tfw, double hfw) {
        int A1 = 0, A2 = 0, A3 = 0;
        double B1, B2, B3;
        double scad;
        if (Type == 1) {
            A1 = 1;
        }
        if (Type == 2) {
            A2 = 1;
        }
        if (Type == 3) {
            A3 = 1;
        }
        B1 = (Aw / (dxw / gammaw) + Math.max(Fw, 0)) * Tw / Vol;
        B2 = -(Aw + Math.max(Fw, 0) * dxw / gammaw) * qw / Vol;
        B3 = Aw / Vol * Tfw / ((dxw / gammaw) + (1 / hfw));//有待处理
        scad = A1 * B1 + A2 * B2 + A3 * B3;
        return scad;
    }

    /**
     *
     * @param Type 类型
     * @param Ae 东边界面积
     * @param Vol 控制容积体积
     * @param dxe 到控制界面外点的距离
     * @param gammae 扩散系数
     * @param Fe 东边界通量
     * @param Te 东边界定值
     * @param qe 热流密度
     * @param Tfe 东边界环境温度
     * @param hfe 换热系数
     * @return 返回附加源项常数部分
     */
    static public double Scad2(int Type, double Ae, double Vol, double dxe,
            double gammae, double Fe, double Te, double qe, double Tfe, double hfe) {
        int A1 = 0, A2 = 0, A3 = 0;
        double B1, B2, B3;
        double scad;
        if (Type == 1) {
            A1 = 1;
        }
        if (Type == 2) {
            A2 = 1;
        }
        if (Type == 3) {
            A3 = 1;
        }
        B1 = (Ae / (dxe / gammae) + Math.max(-Fe, 0)) * Te / Vol;
        B2 = (Ae + Math.max(-Fe, 0) * dxe / gammae) * qe / Vol;
        B3 = Ae / Vol * Tfe / ((dxe / gammae) + (1 / hfe));//有待处理
        scad = A1 * B1 + A2 * B2 + A3 * B3;
        return scad;
    }

    /**
     * boundary update
     *
     * @param phi
     * @param bound
     */
    static void bouUp(double[][][] phi, Boundary bound) {
        //X
        int k = 1;
        for (int j = 1; j < phi[0].length - 1; j++) {
            if (bound.getType(0, j, k) == 2) {
                phi[0][j][k] = phi[1][j][k];
            }
            if (bound.getType(phi.length - 1, j, k) == 2) {
                phi[phi.length - 1][j][k] = phi[phi.length - 2][j][k];
            }
        }
        //Y
        for (int i = 1; i < phi.length - 1; i++) {
            if (bound.getType(i, 0, k) == 2) {
                phi[i][0][k] = phi[i][1][k];
            }
            if (bound.getType(i, phi[0].length - 1, k) == 2) {
                phi[i][phi[0].length - 1][k] = phi[i][phi[0].length - 2][k];
            }
        }
    }

    BoundaryRegion getBoundaryRegion();

}
