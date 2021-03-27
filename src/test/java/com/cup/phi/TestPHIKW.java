/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.phi;

import java.io.FileNotFoundException;
import static java.lang.Math.PI;
import org.junit.Test;

/**
 * @see de Sampaio, P.A.B., J.L.H. Faccini and J. Su, Modelling of stratified
 * gas–liquid two-phase flow in horizontal circular pipes. International Journal
 * of Heat and Mass Transfer, 2008. 51(11-12): p. 2752-2761. 液-液
 * <n>整个过程中假设密度不随着温度变化<n>
 * <n>1.双流体温度场计算</n>
 * <n>2.引入对流项的影响</n>
 * <n>3.未引入时间的影响</n>
 * 实验数据：Kumara W A S, Halvorsen B M, Melaaen M C.Particle image velocimetry for
 * characterizing the flow structure of oil–water flow in horizontal and
 * slightly inclined pipes. Chemical Engineering Science. 2010,65(15):4332-4349.
 * @author winsway
 */
public class TestPHIKW {

    @Test
    public void test1() throws FileNotFoundException, CloneNotSupportedException {
        double Radia = 0.056 / 2.0;
        double Qwater = 0.25 * PI * Radia * Radia;
        double Qoil = 0.25 * PI * Radia * Radia;
        double TOil = 300.15, TWater = 300.15;
        double[] porO = {790, 1.64e-3, 0.3, 2000};
        double[] porW = {996, 1.0e-3, 0.6, 4186};
        int index = 2;
        //全局控制输出位置
        String position = "./TestPHIKW/";
        String name = "测试单相kw/";
        //构造函数
        KW2 test = new KW2(Qoil, Qwater, TOil, TWater, Radia);
        test.position = position;
        test.Title = name;
        //
        test.setFluid(porO, porW);
        test.application(index);
    }

}
