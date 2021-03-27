/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.mesh.factory;

import com.winswe.io.IOobject;
import com.winswe.mesh.Structed2D;
import org.junit.Test;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 */
public class BipolarCoordianteMeshTest {

    public BipolarCoordianteMeshTest() {
    }

    /**
     * Test of gamma method, of class BipolarCoordianteMesh.
     */
    @Test
    public void testGamma() {
        System.out.println("gamma");
        double x = 0.0;
        double y = 0.0;
        BipolarCoordianteMesh instance = null;
        double expResult = 0.0;
        double result = instance.gamma(x, y);

    }

    /**
     * Test of blockMesh method, of class BipolarCoordianteMesh.
     */
    @Test
    public void testBlockMesh() {
        System.out.println("blockMesh");
        String path = "./tutorials/mesh";
        String caseName = "BipolarCoordianteMesh";
        IOobject instance = new IOobject(path, caseName);
        System.out.println("读取配置文件");
        instance.readConfigure();
        System.out.println("生成网格");
        Structed2D mesh = new BipolarCoordianteMesh(instance);
        mesh.blockMesh();
        System.out.println("输出网格文件");
        instance.outPutMesh(mesh);
    }

}
