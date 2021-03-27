package com.winswe.mesh.element;

import com.winswe.math.basic.var.Vector;

public class Shape {
    public final double volume;
    public final Vector centroid;

    public Shape(double volume, Vector centroid) {
        this.volume = volume;
        this.centroid = centroid;
    }

    @Override
    public String toString() {
        return "Shape{" +
                "volume=" + volume +
                ", centroid=" + centroid +
                '}';
    }
}
