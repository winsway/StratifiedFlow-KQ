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
 * 测试湍流+界面摩擦 测试刘楠楠数据
 * <n>1.采用新的方式更改界面更新方法，避免计算速度降低，同时考虑界面摩擦阻力</n>
 *
 * @author winsway
 */
public class TestTF4Wave1 {

    @Test
    public void test1() throws FileNotFoundException {
        String position, Title;
        position = "/D:/winsway/";
        Title = "湍流+界面摩擦+刘楠楠1/";
        double dpdz = 315.0, yint = 0.40;
        double Radia = 0.040 / 2.0;
        double Qwater = 0.3 * PI * Radia * Radia;
        double Qoil = 0.7 * PI * Radia * Radia;
        double[] porO = {824, 5.0e-3, 0.126, 2000};
        double[] porW = {1000, 1.0e-3, 0.6, 4186};
//        
        double F, G;
        double Fy, Fx, Gy, Gx;
        double xold = dpdz, yold = yint;
        double xnew, ynew;
        double errF, errG;
//
        komegaTFOWAngeli2Wave1 test = new komegaTFOWAngeli2Wave1(Qoil, Qwater, Radia);
        test.position = position;
        test.Title = Title;
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
        do {
            System.out.println("Qwater = " + Qwater + " Qoil = " + Qoil);
            System.out.println("F1   " + "xold =" + xold + " yold =" + yold);
            pw.printf("F1 xold =%e\t yold =%e\t\n", xold, yold);
            pw.flush();
            komegaTFOWAngeli2Wave1 F1 = new komegaTFOWAngeli2Wave1(Qoil, Qwater, Radia);
            F1.position = position;
            F1.Title = Title;
            F1.setFluid(porO, porW);
            F1.application(xold, yold, 5);
//            
            F = F1.F(Qwater);
            G = F1.G(Qoil);
            pw.printf("Qcwater = %e\t Qcoil = %e\t\n", F + Qwater, G + Qoil);
            pw.flush();
//            
            System.out.println("Fx1   " + "xold * 1.01 =" + xold * 1.01 + " yold =" + yold);
            pw.printf("Fx1 xold * 1.01 =%e\t yold =%e\t\n", xold * 1.01, yold);
            pw.flush();
            komegaTFOWAngeli2Wave1 Fx1 = new komegaTFOWAngeli2Wave1(Qoil, Qwater, Radia);
            Fx1.position = position;
            Fx1.Title = Title;
            Fx1.setFluid(porO, porW);
            Fx1.application(xold * 1.01, yold, 5);
//            
            System.out.println("Fy1   " + " xold =" + xold + " yold * 1.01 =" + yold * 1.01);
            pw.printf("Fy1 xold =%e\t yold * 1.01 =%e\t\n", xold, yold * 1.01);
            pw.flush();
            komegaTFOWAngeli2Wave1 Fy1 = new komegaTFOWAngeli2Wave1(Qoil, Qwater, Radia);
            Fy1.position = position;
            Fy1.Title = Title;
            Fy1.setFluid(porO, porW);
            Fy1.application(xold, yold * 1.01, 5);
//            
            Fy = (Fy1.F(Qwater) - F) / (0.01 * yold);
            Gy = (Fy1.G(Qoil) - G) / (0.01 * yold);
            Fx = (Fx1.F(Qwater) - F) / (0.01 * xold);
            Gx = (Fx1.G(Qoil) - G) / (0.01 * xold);
            xnew = xold + (G * Fy - F * Gy) / (Fx * Gy - Gx * Fy);
            ynew = yold + (F * Gx - G * Fx) / (Fx * Gy - Gx * Fy);
            xold = xnew;
            yold = ynew;
            System.out.println("xnew = " + xnew + " ynew = " + ynew);
            errF = abs((F) / Qwater);
            errG = abs((G) / Qoil);
            System.out.println("errF = " + errF + " errG = " + errG);
        } while (errF >= 1e-3 || errG >= 1e-3);
        pw.close();
    }

}
