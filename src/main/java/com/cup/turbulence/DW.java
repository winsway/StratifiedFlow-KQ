/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.turbulence;

import com.cup.mesh.Mesh;
import com.cup.util.Tool;
import com.winswe.finiteVolume.geom.Point;

/**
 *
 * @author winsw
 */
public class DW {

    /**
     * 索引
     */
    int x, y, z;
    /**
     * 距离
     */
    double distance;

    public void cacD(Mesh mesh, Point p) {
        double d1;
        double d2 = 1e30;
        for (int J = 0; J <= mesh.getNY() + 1; J++) {
            for (int I = 0; I <= mesh.getNX() + 1; I++) {
                if (I == 0 || I == (mesh.getNX() + 1) || J == 0 || J == (mesh.getNY() + 1)) {
                    d1 = Tool.distance(
                            mesh.realX(mesh.gettPx()[I], mesh.gettPy()[J]),
                            mesh.realY(mesh.gettPx()[I], mesh.gettPy()[J]),
                            p.x, p.y
                    );
                    if (d1 < d2) {
                        x = I;
                        y = J;
                        distance = d1;
                        d2 = d1;
                    }
                }
            }
        }
    }

    public void cacD(Mesh mesh, int x1, int y1) {
        double d1;
        double d2 = 1e30;
        for (int J = 0; J <= mesh.getNY() + 1; J++) {
            for (int I = 0; I <= mesh.getNX() + 1; I++) {
                if (I == 0 || I == (mesh.getNX() + 1) || J == 0 || J == (mesh.getNY() + 1)) {
                    d1 = Tool.distance(
                            mesh.realX(mesh.gettPx()[I], mesh.gettPy()[J]),
                            mesh.realY(mesh.gettPx()[I], mesh.gettPy()[J]),
                            mesh.realX(mesh.gettPx()[x1], mesh.gettPy()[y1]),
                            mesh.realY(mesh.gettPx()[x1], mesh.gettPy()[y1])
                    );
                    if (d1 < d2) {
                        x = I;
                        y = J;
                        distance = d1;
                        d2 = d1;
                    }
                }
            }
        }
    }

    public int getX() {
        return x;
    }

    public int getY() {
        return y;
    }

    public int getZ() {
        return z;
    }

    public double getDistance() {
        return distance;
    }

}
