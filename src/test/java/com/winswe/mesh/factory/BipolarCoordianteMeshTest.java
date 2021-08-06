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
 * test bipolar coordinate mesh
 *
 * @author winswe <halo.winswe@gmail.com>
 */
public class BipolarCoordianteMeshTest {

    /**
     * Test of blockMesh method, of class BipolarCoordianteMesh.
     */
    @Test
    public void testBlockMesh() {
        System.out.println("blockMesh");
        String path = "./tutorials/mesh";
        String caseName = "BipolarCoordianteMesh";
        IOobject instance = new IOobject(path, caseName);
        System.out.println("read configure file");
        instance.readConfigure();
        System.out.println("generate mesh");
        Structed2D mesh = new BipolarCoordianteMesh(instance);
        mesh.blockMesh();
        System.out.println("output mesh file");
        instance.outPutMesh(mesh);
    }

    @Test
    public void forDifferentLiqudLevel() {
        System.out.println("blockMesh");
        String path = "./tutorials/mesh";
        String caseName = "for different liquid level";
        IOobject instance = new IOobject(path, caseName);
        System.out.println("read configure file");
        instance.readConfigure();
        System.out.println("generate mesh");
        Structed2D mesh = new BipolarCoordianteMesh(instance, 0.3);
        mesh.blockMesh();
        System.out.println("output mesh file");
        instance.outPutMesh(mesh);
    }

}
