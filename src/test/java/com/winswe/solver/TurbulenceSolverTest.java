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
 * @date 2021年3月7日 下午5:32:44
 */
public class TurbulenceSolverTest {

    @Test
    public void waterKEpsilonFlow() {
        System.out.println("readConfigure");
        String path = "./tutorials/case/single/turbulence";
        String caseName = "kepsilon";
        TurbulenceSolver turbulenceSolver = new TurbulenceSolver(path, caseName);
        turbulenceSolver.readConfigure();
        //for mesh
        turbulenceSolver.createMesh();
        turbulenceSolver.outPutMesh();
        //for field
        turbulenceSolver.createField();
        turbulenceSolver.solve();
        turbulenceSolver.outPutFields();
    }

    @Test
    public void waterKOmegaFlow() {
        System.out.println("readConfigure");
        String path = "./tutorials/case/single/turbulence";
        String caseName = "komega";
        TurbulenceSolver turbulenceSolver = new TurbulenceSolver(path, caseName);
        turbulenceSolver.readConfigure();
        //for mesh
        turbulenceSolver.createMesh();
        turbulenceSolver.outPutMesh();
        //for field
        turbulenceSolver.createField();
        turbulenceSolver.solve();
        turbulenceSolver.outPutFields();
    }

    @Test
    public void waterSSTKOmegaFlow() {
        System.out.println("readConfigure");
        String path = "./tutorials/case/single/turbulence";
        String caseName = "sstkomega";
        TurbulenceSolver turbulenceSolver = new TurbulenceSolver(path, caseName);
        turbulenceSolver.readConfigure();
        //for mesh
        turbulenceSolver.createMesh();
        turbulenceSolver.outPutMesh();
        //for field
        turbulenceSolver.createField();
        turbulenceSolver.solve();
        turbulenceSolver.outPutFields();
    }
}
