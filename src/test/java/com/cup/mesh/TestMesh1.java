/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.mesh;

import com.winswe.io.Writer;
import com.cup.mesh.factory.BiSameV2;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import static java.lang.Math.PI;
import static java.lang.Math.acos;
import org.junit.Test;

/**
 * 给出物理平面和计算平面的网格
 *
 * @author winsway
 */
public class TestMesh1 {

    public String position, Title;

    @Test
    public void BiSameV2() throws FileNotFoundException {
        Mesh mesh;
        BiSameV2 test = new BiSameV2(12, Math.PI, 100, 100);
        test.radia = 0.125;
        double D = test.radia * 2.0;
        double hl = 0.5;
        hl = D * hl;
        double gamma = acos(1.0 - 2 * hl / D);
        test.theta1 = gamma;
        test.theta2 = PI;
        mesh = test;
        mesh.blockMesh();
        TestMesh1 print = new TestMesh1();
        print.position = "/D:/winsway/";
        print.Title = "Mesh";
        print.RESULT(mesh);
    }

    public void RESULT(Mesh mesh) throws FileNotFoundException {
        String dirName = position + this.toString() + "/";
        Writer.createDir(dirName);
        String fileName = dirName + "Mesh1" + "." + "DAT";
        Writer.createFile(fileName);
        PrintWriter RES1
                = new PrintWriter(new OutputStreamWriter(new FileOutputStream(fileName)));
        RES1.println("Title=" + "\"" + "Field" + "\"");
        String var = "Variables=\"X\",\"Y\"";
        RES1.println(var);
//        
        RES1.println("Zone"
                + " I=" + (mesh.getNX() + 2)
                + " J=" + (mesh.getNY() + 2)
                + " F=POINT");
        for (int J = 0; J <= mesh.getNY() + 1; ++J) {
            for (int I = 0; I <= mesh.getNX() + 1; ++I) {
                RES1.printf("%16.6E\t", mesh.realX(mesh.gettPx()[I], mesh.gettPy()[J]));
                RES1.printf("%16.6E\t", mesh.realY(mesh.gettPx()[I], mesh.gettPy()[J]));
                RES1.println();
            }
        }
        RES1.close();
//        
        fileName = dirName + "Mesh2" + "." + "DAT";
        PrintWriter RES2
                = new PrintWriter(new OutputStreamWriter(new FileOutputStream(fileName)));
        RES2.println("Title=" + "\"" + "Field" + "\"");
        var = "Variables=\"X\",\"Y\"";
        RES2.println(var);
//        
        RES2.println("Zone"
                + " I=" + (mesh.getNX() + 2)
                + " J=" + (mesh.getNY() + 2)
                + " F=POINT");
        for (int J = 0; J <= mesh.getNY() + 1; ++J) {
            for (int I = 0; I <= mesh.getNX() + 1; ++I) {
                RES2.printf("%16.6E\t", mesh.gettPx()[I]);
                RES2.printf("%16.6E\t", mesh.gettPy()[J]);
                RES2.println();
            }
        }
        RES2.close();
    }

    @Override
    public String toString() {
        return getClass().getName() + Title;
    }
}
