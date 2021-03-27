package com.winswe.finiteVolume.geom;

import com.winswe.math.basic.var.Vector;

/**
 * An interface for defining various VTK Geometry.
 */
public interface Geometry {

    /**
     * Ordered points as defined by the VTK documentation
     *
     * @return An array of <code>Point</code> objects following the sequence as
     * defined by the VTK geometry.
     */
    Point[] points();

    /**
     *
     * @return 长度
     */
    double length();

    /**
     *
     * @return 面积
     */
    double area();

    /**
     *
     * @return 体积
     */
    double volume();

    /**
     *
     * @return 质心坐标
     */
    Point centroid();

    /**
     *
     * @return 单位法向量
     */
    Vector unitNormal();
}
