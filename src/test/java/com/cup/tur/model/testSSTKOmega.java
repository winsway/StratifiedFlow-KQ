/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.tur.model;

import com.cup.boundary.factory.RobinBC;
import java.io.FileNotFoundException;
import static java.lang.Math.PI;
import org.junit.Test;

/**
 * 测试SSTK-Omega模型的单相
 *
 * @author winsway
 */
public class testSSTKOmega {

    @Test
    public void test() throws FileNotFoundException {
        double Radia = 0.056 / 2.0;
        String position = "./", Title = "测试水单相";
        double Um = 0.5;
        double lamOil = 0.5, lamWater = 0.5;
        double Qwater = Um * lamWater * PI * Radia * Radia;
        double Qoil = Um * lamOil * PI * Radia * Radia;
        double[] porO = {790, 1.64e-3, 0.3, 2000};
        double[] porW = {996, 1.0e-3, 0.6, 4186};
        double dpdz = 200, hl = 0.5;
        //        
        int index = 1;
        //应该引入反算Re的压降;   
        int W = 0, K = 1, Omega = 2, T = 3;
        RobinBC[] bc = new RobinBC[4];
        bc[W] = new RobinBC(0.0, 1.0, 0.0, 1.0);
        bc[K] = new RobinBC(0.0, 1.0, 0.0, 1.0);
        bc[Omega] = new RobinBC(0.0, 1.0, 0.0, 1.0);
        bc[T] = new RobinBC(59.0, 15.0, 375.0, 1.0);
        //        
        SSTKOmega test = new SSTKOmega(position, Title);
        test.setFluid(porW, porW);
        test.setBoundary(bc);
        test.application(dpdz, hl, Radia, index);

    }

}
