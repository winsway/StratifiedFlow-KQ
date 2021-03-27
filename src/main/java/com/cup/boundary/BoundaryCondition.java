/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.boundary;

import static java.lang.Math.max;
import com.cup.boundary.factory.AbstractBC;
import com.cup.boundary.factory.RobinBC;

/**
 *
 * @author winsway
 */
public interface BoundaryCondition {

    /**
     *
     * @return
     */
    public AbstractBC getAbstractBC();

    /**
     *
     * @param Type 边界条件类型
     * @param Dw 面扩散系数w
     * @param Fw 面通量w
     * @param dxw W和P之间距离
     * @param vol 体积
     * @param bc 边界条件
     * @return Spad
     */
    static public double Spadw(int Type, double Dw, double Fw, double dxw, double vol, RobinBC bc) {
        if (bc == null || Type == 0 || Type == 7) {
            return 0;
        } else {
            double aw = Dw + max(Fw, 0);
            double a, b;
            a = bc.a;
            b = bc.b;
            double temp = -aw * (b * dxw) / (a + b * dxw) / vol;
            return temp;
        }
    }

    /**
     *
     * @param Type 边界条件类型
     * @param De 面扩散系数e
     * @param Fe 面通量e
     * @param dxe E和P之间距离
     * @param vol 体积
     * @param bc 边界条件
     * @return Spad
     */
    static public double Spade(int Type, double De, double Fe, double dxe, double vol, RobinBC bc) {
        if (bc == null || Type == 0 || Type == 7) {
            return 0;
        } else {
            double ae = De + max(-Fe, 0);
            double a, b;
            a = bc.a;
            b = bc.b;
            double temp = -ae * (b * dxe) / (a + b * dxe) / vol;
            return temp;
        }
    }

    /**
     *
     * @param Type 边界条件类型
     * @param Dw 面扩散系数w
     * @param Fw 面通量w
     * @param dxw W和P之间距离
     * @param vol 体积
     * @param bc 边界条件
     * @return Spad
     */
    static public double Scadw(int Type, double Dw, double Fw, double dxw, double vol, RobinBC bc) {
        if (bc == null || Type == 0 || Type == 7) {
            return 0;
        } else {
            double aw = Dw + max(Fw, 0);
            double a, b, M;
            a = bc.a;
            b = bc.b;
            M = bc.M;
            double temp = aw * (M * dxw) / (a + b * dxw) / vol;
            return temp;
        }
    }

    /**
     *
     * @param Type 边界条件类型
     * @param De 面扩散系数e
     * @param Fe 面通量e
     * @param dxe E和P之间距离
     * @param vol 体积
     * @param bc 边界条件
     * @return Spad
     */
    static public double Scade(int Type, double De, double Fe, double dxe, double vol, RobinBC bc) {
        if (bc == null || Type == 0 || Type == 7) {
            return 0;
        } else {
            double ae = De + max(-Fe, 0);
            double a, b, M;
            a = bc.a;
            b = bc.b;
            M = bc.M;
            double temp = ae * (M * dxe) / (a + b * dxe) / vol;
            return temp;
        }
    }

    /**
     * w系数
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
    static public double Spad1(int Type, double Aw, double Vol, double dxw, double gammaw, double Fw, double hfw) {
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
    static public double Spad2(int Type, double Ae, double Vol, double dxe, double gammae, double Fe, double hfe) {
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
    static public double Scad1(int Type, double Aw, double Vol, double dxw, double gammaw, double Fw, double Tw, double qw, double Tfw, double hfw) {
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
    static public double Scad2(int Type, double Ae, double Vol, double dxe, double gammae, double Fe, double Te, double qe, double Tfe, double hfe) {
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

}
