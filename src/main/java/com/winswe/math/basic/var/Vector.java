package com.winswe.math.basic.var;

import com.winswe.finiteVolume.geom.Point;

/**
 * 矢量
 *
 * @author winswe <halo.winswe@gmail.com>
 */
public class Vector {

    private final static double EPS = 1e-12;

    public static final Vector ZERO = new Vector(0, 0, 0);
    public static final Vector ONE = new Vector(1, 1, 1);
    public static final String[] name = {"x", "y", "z"};

    public double[] value = new double[3];

    public Vector(double x, double y, double z) {
        value[0] = x;
        value[1] = y;
        value[2] = z;
    }

    public Vector(Point from, Point to) {
        this(to.x - from.x,
                to.y - from.y,
                to.z - from.z);
    }

    public Point toPoint() {
        return new Point(
                comtValue(0),
                comtValue(1),
                comtValue(2)
        );
    }

    /**
     * 在本身基础上相加
     *
     * @param other
     */
    public void addSelf(Vector other) {
        setComValue(0, comtValue(0) + other.comtValue(0));
        setComValue(1, comtValue(1) + other.comtValue(1));
        setComValue(2, comtValue(2) + other.comtValue(2));
    }

    public Vector add(Vector other) {
        return new Vector(
                comtValue(0) + other.comtValue(0),
                comtValue(1) + other.comtValue(1),
                comtValue(2) + other.comtValue(2)
        );
    }

    /**
     * 在本身基础上相加
     *
     * @param other
     */
    public void subSelf(Vector other) {
        setComValue(0, comtValue(0) - other.comtValue(0));
        setComValue(1, comtValue(1) - other.comtValue(1));
        setComValue(2, comtValue(2) - other.comtValue(2));
    }

    public Vector sub(Vector other) {
        return new Vector(
                comtValue(0) - other.comtValue(0),
                comtValue(1) - other.comtValue(1),
                comtValue(2) - other.comtValue(2)
        );
    }

    /**
     *
     * @param scalar 自身乘上标量值
     */
    public void multSelf(double scalar) {
        setComValue(0, comtValue(0) * scalar);
        setComValue(1, comtValue(1) * scalar);
        setComValue(2, comtValue(2) * scalar);
    }

    public Vector mult(double scalar) {
        return new Vector(
                comtValue(0) * scalar,
                comtValue(1) * scalar,
                comtValue(2) * scalar
        );
    }

    public double magSqr() {
        return comtValue(0) * comtValue(0)
                + comtValue(1) * comtValue(1)
                + comtValue(2) * comtValue(2);
    }

    public double mag() {
        return Math.sqrt(magSqr());
    }

    public Vector unit() {
        double mag = mag();
        if (mag < EPS) {
            throw new ArithmeticException("Divide by zero.");
        }

        return new Vector(
                comtValue(0) / mag,
                comtValue(1) / mag,
                comtValue(2) / mag
        );
    }

    public double dotDouble(Vector other) {
        return (this.comtValue(0) * other.comtValue(0)
                + this.comtValue(1) * other.comtValue(1)
                + this.comtValue(2) * other.comtValue(2));
    }

    @Override
    public String toString() {
        return "("
                + comtValue(0) + " "
                + comtValue(1) + " "
                + comtValue(2)
                + ")";
    }

    @Override
    public boolean equals(final Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof Vector)) {
            return false;
        }
        Vector other = (Vector) obj;
        if (this.comtValue(0) != other.comtValue(0)
                && this.comtValue(1) != other.comtValue(1)
                && this.comtValue(2) != other.comtValue(2)) {
            return false;
        }
        return true;
    }

    public String componentNames(int ncmpt) {
        return name[ncmpt];
    }

    public double comtValue(int ncmpt) {
        return value[ncmpt];
    }

    public void setComValue(int ncmpt, double value) {
        if (ncmpt <= 2) {
            this.value[ncmpt] = value;
        } else {
            throw new UnsupportedOperationException("分量设定错误!");
        }
    }

    public Vector cmptMultiply(Vector arg0, Vector arg1) {
        return new Vector(
                arg0.comtValue(0) * arg1.comtValue(0),
                arg0.comtValue(1) * arg1.comtValue(1),
                arg0.comtValue(2) * arg1.comtValue(2)
        );
    }

    public void set(Vector arg1) {
        this.set(
                arg1.comtValue(0),
                arg1.comtValue(1),
                arg1.comtValue(2)
        );
    }

    /**
     * 设定值
     *
     * @param x
     * @param y
     * @param z
     */
    public void set(double x, double y, double z) {
        setComValue(0, x);
        setComValue(1, y);
        setComValue(2, z);
    }

    /**
     *
     * @return new Vector
     * @throws CloneNotSupportedException
     */
    @Override
    public Vector clone() throws CloneNotSupportedException {
        return new Vector(value[0], value[1], value[2]);
    }

    public String className() {
        return "vector";
    }

    public Vector cross(Vector arg1) {
        Vector arg0 = this;
        double xComp
                = arg0.comtValue(1) * arg1.comtValue(2)
                - arg1.comtValue(1) * arg0.comtValue(2);
        double yComp
                = arg1.comtValue(0) * arg0.comtValue(2)
                - arg0.comtValue(0) * arg1.comtValue(2);
        double zComp
                = arg0.comtValue(0) * arg1.comtValue(1)
                - arg1.comtValue(0) * arg0.comtValue(1);
        return new Vector(xComp, yComp, zComp);
    }

}
