/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.turbulence;

import com.cup.mesh.Mesh;
import com.cup.util.Tool;

/**
 * 计算到壁面的距离，得到对应的点的坐标。
 *
 * @author winsway
 */
public class Distance {

    int x, y, z;
    double distance;

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

    public void cacDint(Mesh mesh, double YP, int x1, int y1, String var
    ) {
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
                    if (J == 0 && var == "Water") {
                        d1 = d1 + YP;
                    }
                    if (J == (mesh.getNY() + 1) && var == "Oil") {
                        d1 = d1 + YP;
                    }
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

    /**
     *
     * @return DN
     */
    public double getDistance() {
        return distance;
    }

}
