/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.io;

import com.winswe.io.Writer;
import com.cup.field.Field;
import com.cup.mesh.Mesh;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;

/**
 *
 * @author winsway
 */
public final class Output {

    Field U;
    Mesh mesh;
    String dir, title;
    double[][][] UL;

    public Output(Field U, Mesh mesh, String dir, String title) {
        this.U = U;
        this.mesh = mesh;
        this.dir = dir;
        this.title = title;
        this.output1();
    }

    public Output(double[][][] U, Mesh mesh, String dir, String title) {
        this.UL = U;
        this.mesh = mesh;
        this.dir = dir;
        this.title = title;
        this.output2();
    }

    public void output1() {
        try {
            writerMesh(dir,
                    "dat",
                    title,
                    U.getNewField(),
                    mesh.gettPx(),
                    mesh.gettPy(),
                    mesh.gettPz());
        } catch (FileNotFoundException ex) {
            System.out.println("Can not output field file, please check!");
        }
    }

    public void output2() {
        try {
            writerMesh(dir,
                    "dat",
                    title,
                    UL,
                    mesh.gettPx(),
                    mesh.gettPy(),
                    mesh.gettPz());
        } catch (FileNotFoundException ex) {
            System.out.println("Can not output field file, please check!");
        }
    }

    void writerMesh(String dirName, String formate, String titlename,
            double[][][] phi, double[] x, double[] y, double[] z)
            throws FileNotFoundException {
        //Create direction
        Writer.createDir(dirName);
        //Create file
        String str = String.format("%s", titlename);
        String fileName = dirName + str + "." + formate;
        Writer.createFile(fileName);
        inputField(fileName, titlename, phi, x, y, z);
    }

    void inputField(String fileName, String titlename,
            double[][][] phi, double[] x, double[] y, double[] z)
            throws FileNotFoundException {
        PrintWriter pw
                = new PrintWriter(new OutputStreamWriter(new FileOutputStream(fileName)));
        //TITLE
        pw.println("Title=" + "\"" + titlename + "\"");
        String var = "Variables=\"X\",\"Y\",\"Z\"";
        pw.println(var + "\"" + titlename + "\"");
        pw.println("Zone"
                + " I=" + phi.length
                + " J=" + phi[0].length
                + " K=" + (phi[0][0].length - 2)
                + " F=POINT");
//
        for (int k = 1; k < z.length - 1; ++k) {
            for (int j = 0; j < y.length; ++j) {
                for (int i = 0; i < x.length; ++i) {
                    pw.printf("%16.6f\t", mesh.realX(x[i], y[j]));
                    pw.printf("%16.6f\t", mesh.realY(x[i], y[j]));
                    pw.printf("%16.6f\t", z[k]);
                    pw.printf("%16.6f\t", phi[i][j][k]);
                    pw.println();
                }
            }
        }
        pw.close();
    }

}
