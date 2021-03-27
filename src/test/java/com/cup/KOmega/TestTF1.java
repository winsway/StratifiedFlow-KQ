/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.KOmega;

import com.winswe.io.Writer;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import org.junit.Test;

/**
 *
 * @author winsw
 */
public class TestTF1 {

    /**
     * 测试单相双流体层流方程
     *
     * @throws FileNotFoundException
     */
    @Test
    public void test1() throws FileNotFoundException {
        double dpdz = 252.0, yint = 0.50;
        double Radia = 0.0243 / 2.0;
        double[] porO = {1000, 1e-3, 0.1260, 2000};
        double[] porW = {1000, 1e-3, 0.1260, 2000};
        komegaTFOWAngeli1 test = new komegaTFOWAngeli1(0, 0, Radia);
        test.position = "/D:/winsway/";
        test.Title = "层流/";
        test.setFluid(porO, porW);
        test.application(dpdz, yint, 3);
    }

}
