/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.io;

//import com.cup.field.Scalarfield;
/**
 *
 * @author winsway
 */
public class DisplayField {

    static public void displayField2D(double[][] temp) {
        System.out.println(temp.toString());
        for (int i = 0; i < temp.length; i++) {
            for (int j = 0; j < temp[i].length; j++) {
                System.out.printf("%8.4f   ", temp[i][j]);
            }
            System.out.println("");
        }
    }

    static public void displayField3D(double[][][] temp) {
        System.out.println(temp.toString());
        for (int i = 0; i < temp.length; i++) {
            for (int j = 0; j < temp[i].length; j++) {
                System.out.printf("%8.4f   ", temp[i][j][1]);
            }
            System.out.println("");
        }
    }

    static public void displayField3D(double[][][] temp, String name) {
        System.out.println(name);
        for (int i = 0; i < temp.length; i++) {
            for (int j = 0; j < temp[i].length; j++) {
                System.out.printf("%8.4f   ", temp[i][j][1]);
            }
            System.out.println("");
        }
    }

}
