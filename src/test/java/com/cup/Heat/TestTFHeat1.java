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
 * 测试温度场，没有问题
 *
 * @author winsway
 */
public class TestTFHeat1 {

    @Test
    public void test1() throws FileNotFoundException {
        double dpdz = 252.0, yint = 0.50;
        double Radia = 0.0243 / 2.0;
        double Qwater = 0.22 * PI * Radia * Radia;
        double Qoil = 0.22 * PI * Radia * Radia;
        double[] porO = {996, 8.6e-4, 0.6, 4186};
        double[] porW = {996, 8.6e-4, 0.6, 4186};
        String position = "./";
        String name = "TF温度耦合/";
//              
        double xold = dpdz, yold = yint;
//
        TFHeat1 test = new TFHeat1(Qoil, Qwater, Radia);
        test.position = position;
        test.Title = name;
//       
        String dirName = test.position + test.toString();
        String title = "pressure and liquid high";
        Writer.createDir(dirName);
        String fileName = dirName + title + "." + "txt";
        Writer.createFile(fileName);
//        
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
