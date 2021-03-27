package com.winswe.finiteVolume.geom;

import com.winswe.math.basic.var.Vector;
import static java.lang.Math.sqrt;

public class Point {

    public final double x, y, z;

    public Point(Point temp) {
        this.x = temp.x;
        this.y = temp.y;
        this.z = temp.z;
    }

    public Point(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public Vector toVector() {
        return new Vector(x, y, z);
    }

    public double distance(Point other) {
        double dx = other.x - this.x;
        double dy = other.y - this.y;
        double dz = other.z - this.z;

        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    public Point midPoint(Point A, Point B) {
        return new Point(
                0.5 * (A.x + B.x),
                0.5 * (A.y + B.y),
                0.5 * (A.z + B.z)
        );
    }

    @Override
    public String toString() {
        return "Point{"
                + "x=" + x
                + ", y=" + y
                + ", z=" + z
                + '}';
    }
}
