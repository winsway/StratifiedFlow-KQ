package com.winswe.mesh.element;

import com.winswe.finiteVolume.geom.Point;
import com.winswe.math.basic.var.Vector;

import java.util.ArrayList;

public class Node {

    public final Vector p;
    public final ArrayList<Cell> neighbors;
    private int index;

    public Node(double x, double y, double z) {
        index = -1;
        p = new Vector(x, y, z);
        this.neighbors = new ArrayList<>();
    }

    public Node(Vector p) {
        this(
                p.comtValue(0),
                p.comtValue(1),
                p.comtValue(2)
        );
    }

    /**
     *
     * @param index 设定索引
     */
    public void setIndex(int index) {
        if (this.index == -1) {
            this.index = index;
        } else {
            throw new IllegalStateException("The index can be set only once.");
        }
    }

    public int index() {
        return this.index;
    }

    /**
     *
     * @return 转化为Point
     */
    public Point location() {
        return p.toPoint();
    }

    @Override
    public String toString() {
        return "Node{"
                + "x=" + p.comtValue(0)
                + ", y=" + p.comtValue(1)
                + ", z=" + p.comtValue(2)
                + '}';
    }
}
