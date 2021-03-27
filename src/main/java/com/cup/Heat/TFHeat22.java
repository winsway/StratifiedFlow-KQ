/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.Heat;

import com.cup.boundary.factory.RobinBC;
import java.io.FileNotFoundException;

/**
 *
 * @see de Sampaio, P.A.B., J.L.H. Faccini and J. Su, Modelling of stratified
 * gas–liquid two-phase flow in horizontal circular pipes. International Journal
 * of Heat and Mass Transfer, 2008. 51(11-12): p. 2752-2761. 液-液
 * <n>1.双流体温度场计算</n>
 * <n>2.引入对流项的影响</n>
 * <n>3.未引入时间的影响</n>
 * @author 齐雪宇
 */
public final class TFHeat22 {

    TFHeat2[][] X = new TFHeat2[10][2];
    int old = 0, now = 1;
    /**
     * 项目名称
     */
    String Title;
    /**
     * 文件的位置
     */
    String position;
    /**
     * 物性参数
     */
    double[] testO, testW;
    /**
     * 初始压降和液位高度
     */
    double dpdz, hl;
    /**
     * 圆管半径,M
     */
    double Radia;
    /**
     * 设定入口流量，单位是m³/s
     */
    double Qoil, Qwater;
    /**
     * 设定入口的温度，单位都是K
     */
    double TOil, TWater;

    /**
     *
     * @param Qoil
     * @param Qwater
     * @param Radia
     */
    public TFHeat22(double Qoil, double Qwater, double Radia) {
        this.Qoil = Qoil;
        this.Qwater = Qwater;
        this.Radia = Radia;
    }

    /**
     *
     * @param Qoil 入口油的流量m³/s
     * @param Qwater 入口水的流量m³/s
     * @param TOil 入口油的温度K
     * @param TWater 入口水的温度K
     * @param Radia 管径,m
     */
    public TFHeat22(double Qoil, double Qwater, double TOil, double TWater, double Radia) {
        this.Qoil = Qoil;
        this.Qwater = Qwater;
        this.TOil = TOil;
        this.TWater = TWater;
        this.Radia = Radia;
    }

    void setFluid(double[] testO, double[] testW) {
        this.testO = testO;
        this.testW = testW;
    }

    void initial() {
        dpdz = 0;
        hl = 0.5;
        position = position + this.toString();
        RobinBC BC = new RobinBC(0.62, 1200, 1200 * (273.15 + 0), 1.0);
        for (int i = 0; i < X.length; i++) {
            X[i][old] = new TFHeat2(Qoil, Qwater, TOil, TWater, BC, Radia);
            X[i][old].position = position;
            X[i][old].Title = "Xold" + i + "\\";
            X[i][old].setFluid(testO, testW);
            X[i][old].start(dpdz, hl);
            //
            X[i][now] = new TFHeat2(Qoil, Qwater, TOil, TWater, BC, Radia);
            X[i][now].position = position;
            X[i][now].Title = "Xnow" + i + "\\";
            X[i][now].setFluid(testO, testW);
            X[i][now].start(dpdz, hl);
        }
    }

    public void application(int index) throws FileNotFoundException, CloneNotSupportedException {
        initial();
        if (index == 0) {
            scheme0();
        }
        if (index == 1) {
            scheme1();
        }
        if (index == 2) {
            scheme2();
        }
    }

    int index;

    /**
     * 单独测试温度场 瞬态
     *
     * @throws FileNotFoundException
     * @throws CloneNotSupportedException
     */
    void scheme0() throws FileNotFoundException, CloneNotSupportedException {
        dpdz = 10;
        hl = 0.5;
        double dt = 10.0;
        double dz = 0.0;
        index = 2;
        //只进行温度计算
        X[0][now].setTimeConvection(X[0][old], X[0][old], dt, dz);
        X[0][now].application(dpdz, hl, index);
        X[0][old] = X[0][now].clone();
        //开始引入温度计算
        //需要得到上游的信息
        index = 3;
        X[1][now].setTimeConvection(X[1][old], X[0][old], dt, dz);
        X[1][now].application(dpdz, hl * 1.05, index);
        X[1][old] = X[1][now].clone();
    }

    /**
     * 测试温度场瞬态+对流项+网格变化是否有影响
     *
     * @throws FileNotFoundException
     */
    void scheme1() throws FileNotFoundException, CloneNotSupportedException {
        dpdz = 10;
        hl = 0.5;
        double dt = 10.0;
        double dz = 50.0;
        index = 1;
        //只进行速度计算
        X[0][now].setTimeConvection(X[0][old], X[0][old], dt, dz);
        X[0][now].application(dpdz, hl, index);
        X[0][old] = X[0][now].clone();
        //        
        X[1][old].setTimeConvection(X[1][old], X[0][old], dt, dz);
        X[1][old].application(dpdz, hl * 1.05, index);
        //开始引入温度计算
        //需要得到上游的信息
        index = 3;
        X[1][now].setTimeConvection(X[1][old], X[0][old], dt, dz);
        X[1][now].application(dpdz, hl * 0.95, index);
        X[1][old] = X[1][now].clone();
    }

    /**
     * 引入压力和液位修正算法
     *
     * @throws FileNotFoundException
     * @throws CloneNotSupportedException
     */
    void scheme2() throws FileNotFoundException, CloneNotSupportedException {

    }

    @Override
    public String toString() {
        return getClass().getName() + Title;
    }
}
