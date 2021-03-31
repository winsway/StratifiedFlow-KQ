/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.solver;

import org.junit.Test;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月31日 下午8:49:13
 */
public class TwoPhaseSolverTest {

    @Test
    public void waterFlow() {
        System.out.println("readConfigure");
        String path = "./tutorials/case/twophase";
        String caseName = "laminar";
        TwoPhaseSolver twoPhaseSolver = new TwoPhaseSolver(path, caseName);
        twoPhaseSolver.readConfigure();
        //for mesh
        twoPhaseSolver.createMesh();
        twoPhaseSolver.outPutMesh();
        //for field
        twoPhaseSolver.createField();
        twoPhaseSolver.startIteration();
        twoPhaseSolver.outPutFields();
    }

}
