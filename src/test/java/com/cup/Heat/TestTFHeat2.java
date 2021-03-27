/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.Heat;

import com.winswe.io.Writer;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;

import static java.lang.Math.PI;
import org.junit.Test;

/**
 * @see de Sampaio, P.A.B., J.L.H. Faccini and J. Su, Modelling of stratified
 * gas–liquid two-phase flow in horizontal circular pipes. International Journal
 * of Heat and Mass Transfer, 2008. 51(11-12): p. 2752-2761. 液-液
 * <n>1.双流体温度场计算</n>
 * <n>2.引入对流项的影响</n>
 * <n>3.未引入时间的影响</n>
 *
 * @author winsway
 */
public class TestTFHeat2 {

    @Test
    public void test1() throws FileNotFoundException {
        double dpdz = 252.0, yint = 0.50;
        double Radia = 0.0243 / 2.0;
        double Qwater = 0.22 * PI * Radia * Radia;
        double Qoil = 0.22 * PI * Radia * Radia;
        double[] porO = {996, 8.6e-4, 0.6, 4186};
        double[] porW = {996, 8.6e-4, 0.6, 4186};
//      全局控制输出位置
        String position = "/D:/winsway/";
        String name = "TF温度耦合/";
//      设定dpdz,hl;
        double xold = dpdz, yold = yint;
//      构造函数
        TFHeat2 test = new TFHeat2(Qoil, Qwater, Radia);
        test.position = position;
        test.Title = name;
//      设定输出位置
        String dirName = test.position + test.toString();
        String title = "pressure and liquid high";
        Writer.createDir(dirName);
        String fileName = dirName + title + "." + "txt";
        Writer.createFile(fileName);
//      设定文件的输出
        PrintWriter pw
                = new PrintWriter(new OutputStreamWriter(new FileOutputStream(fileName)));
        pw.printf("Qwater = %e\t Qoil = %e\t\n", Qwater, Qoil);
        pw.flush();

        System.out.println("Qwater = " + Qwater + " Qoil = " + Qoil);
        System.out.println("F1   " + "xold =" + xold + " yold =" + yold);
        pw.printf("F1 xold =%e\t yold =%e\t\n", xold, yold);
        pw.flush();
        //
        test.setFluid(porO, porW);
        test.application(xold, yold, 6);
        pw.close();
    }

}
