/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.Heat;

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
 *
 * @author winsway
 */
public class TestTFHeat22 {

    @Test
    public void test1() throws FileNotFoundException, CloneNotSupportedException {
        double Radia = 0.0243 / 2.0;
        double Qwater = 0.11 * PI * Radia * Radia;
        double Qoil = 0.11 * PI * Radia * Radia;
        double TOil = 300.15, TWater = 300.15;
        double[] porO = {824, 5.0e-3, 0.3, 2000};
        double[] porW = {1000, 1.0e-3, 0.6, 4186};
        int index = 1;
        //全局控制输出位置
        String position = "/D:/winsway/";
        String name = "TF温度耦合+con/";
        //构造函数
        TFHeat22 test = new TFHeat22(Qoil, Qwater, TOil, TWater, Radia);
        test.position = position;
        test.Title = name;
        //
        test.setFluid(porO, porW);
        test.application(index);
    }

}
