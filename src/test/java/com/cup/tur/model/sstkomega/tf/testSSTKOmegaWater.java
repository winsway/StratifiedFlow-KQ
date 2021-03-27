/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.tur.model.sstkomega.tf;

import com.cup.boundary.factory.RobinBC;
import java.io.FileNotFoundException;
import static java.lang.Math.PI;
import org.junit.Test;

/**
 * 测试K-Omega模型LiuNanNan油水两相
 *
 * @author winsway
 */
public class testSSTKOmegaWater {

    @Test
    public void test() throws FileNotFoundException, CloneNotSupportedException {
        String position = "/D:/winsway/", Title = "测试water";
//
        double Radia = 0.0243 / 2.0;
        double Um = 0.411522633744856;
//        
        double lamOil = 0.5, lamWater = 0.5;
        double Qwater = Um * lamWater * PI * Radia * Radia;
        double Qoil = Um * lamOil * PI * Radia * Radia;
        double Toil = 273.15 + 30, Twater = 273.15 + 30;
        double[] porO = {824, 5.00e-3, 0.3, 2000};
        double[] porW = {1000, 1.00e-3, 0.6, 4186};
        //        
        int index = 0;
        //应该引入反算Re的压降;   
        int W = 0, K = 1, Omega = 2, T = 3;
        RobinBC[] bc = new RobinBC[4];
        bc[W] = new RobinBC(0.0, 1.0, 0.0, 1.0);
        bc[K] = new RobinBC(0.0, 1.0, 0.0, 1.0);
        bc[Omega] = new RobinBC(0.0, 1.0, 0.0, 1.0);
        bc[T] = new RobinBC(59.0, 15.0, 375.0, 1.0);
        //        
        Control test = new Control(position, Title);
        test.setFluid(porW, porW);
        test.setBoundary(bc);
        test.setMesh(50, 50);
        test.setQ(Qoil, Qwater);
        test.setT(Toil, Twater);
        test.setRadia(Radia);
        test.application(index);
    }

}