/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.tur.model;

import com.cup.boundary.factory.RobinBC;
import java.io.File;
import java.io.FileNotFoundException;
import static java.lang.Math.PI;
import org.junit.Test;

/**
 * 测试K-Omega模型的单相
 *
 * @author winsway
 */
public class testKOmega {

    @Test
    public void test() throws FileNotFoundException {
        String position = "./", Title = "测试水单相Re=10000Mesh=50";
//
        double Radia = 0.0243 / 2.0;
        double dpdz = 141.585970949543, hl = 0.5;
        double Um = 0.5;
//        
        double lamOil = 0.5, lamWater = 0.5;
        double Qwater = Um * lamWater * PI * Radia * Radia;
        double Qoil = Um * lamOil * PI * Radia * Radia;
        double[] porO = {790, 1.64e-3, 0.3, 2000};
        double[] porW = {1000, 1.0e-3, 0.6, 4186};
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
        KOmega test = new KOmega(position, Title);
        test.setFluid(porW, porW);
        test.setBoundary(bc);
        test.application(dpdz, hl, Radia, index);
    }

    @Test
    public void test1() throws FileNotFoundException {
        String position = "./", Title = "测试水单相Re=7776Mesh=200";
//
        double Radia = 0.0243 / 2.0;
        double dpdz = 88.43995444799775, hl = 0.5;
        double Um = 0.32;
//        
        double lamOil = 0.5, lamWater = 0.5;
        double Qwater = Um * lamWater * PI * Radia * Radia;
        double Qoil = Um * lamOil * PI * Radia * Radia;
        double[] porO = {790, 1.64e-3, 0.3, 2000};
        double[] porW = {1000, 1.0e-3, 0.6, 4186};
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
        KOmega test = new KOmega(position, Title);
        test.setFluid(porW, porW);
        test.setBoundary(bc);
        test.application(dpdz, hl, Radia, index);
    }

}