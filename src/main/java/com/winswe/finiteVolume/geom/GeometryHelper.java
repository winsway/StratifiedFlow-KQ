package com.winswe.finiteVolume.geom;

import com.winswe.math.basic.var.Vector;

/**
 * Helper methods for calculation of geometry properties. Look at my repository
 * on GitHub:
 * <a href="https://github.com/heySourabh/GeometryCalculations">
 * https://github.com/heySourabh/GeometryCalculations
 * </a>
 * for further details and derivations.
 * 为了进行数值计算，建议将实体平移到接近原点的位置，以减少由于大数平方引起的数值误差。
 *
 * @author Sourabh Bhat (heySourabh@gmail.com)
 */
public class GeometryHelper {

    private static class BoundingBox {

        final double xMin, xMax, yMin, yMax, zMin, zMax;

        private BoundingBox(
                double xMin, double xMax,
                double yMin, double yMax,
                double zMin, double zMax) {
            this.xMin = xMin;
            this.xMax = xMax;
            this.yMin = yMin;
            this.yMax = yMax;
            this.zMin = zMin;
            this.zMax = zMax;
        }

        private Point midPoint() {
            return new Point(
                    (xMin + xMax) * 0.5,
                    (yMin + yMax) * 0.5,
                    (zMin + zMax) * 0.5);
        }
    }

    public static class TriGeom {

        final Point[] points;

        public TriGeom(Point p0, Point p1, Point p2) {
            points = new Point[]{p0, p1, p2};
        }

        private Point[] points() {
            return points;
        }
    }

    private static BoundingBox boundingBox(TriGeom... triangles) {
        double xMin = triangles[0].points()[0].x;
        double xMax = triangles[0].points()[0].x;
        double yMin = triangles[0].points()[0].y;
        double yMax = triangles[0].points()[0].y;
        double zMin = triangles[0].points()[0].z;
        double zMax = triangles[0].points()[0].z;
        /*计算得到几何坐标的最大和最小值*/
        for (TriGeom triangle : triangles) {
            for (int j = 0; j < 3; j++) {
                Point p = triangle.points()[j];
                xMin = Math.min(xMin, p.x);
                xMax = Math.max(xMax, p.x);
                yMin = Math.min(yMin, p.y);
                yMax = Math.max(yMax, p.y);
                zMin = Math.min(zMin, p.z);
                zMax = Math.max(zMax, p.z);
            }
        }
        return new BoundingBox(xMin, xMax, yMin, yMax, zMin, zMax);
    }

    /**
     * 平移三角形
     *
     * @param tri 原始三角
     * @param dx 平移x
     * @param dy 平移y
     * @param dz 平移z
     * @return 三角形的新坐标
     */
    private static TriGeom translateTriangle(
            TriGeom tri,
            double dx, double dy, double dz) {
        Point[] p = tri.points();
        return new TriGeom(
                new Point(p[0].x + dx, p[0].y + dy, p[0].z + dz),
                new Point(p[1].x + dx, p[1].y + dy, p[1].z + dz),
                new Point(p[2].x + dx, p[2].y + dy, p[2].z + dz));
    }

    /**
     *
     * @param triangles
     * @param dx
     * @param dy
     * @param dz
     * @return 平移三角形的位置
     */
    private static TriGeom[] translateTriangles(
            TriGeom[] triangles,
            double dx, double dy, double dz) {
        TriGeom[] newTriangles = new TriGeom[triangles.length];
        for (int i = 0; i < triangles.length; i++) {
            newTriangles[i] = translateTriangle(triangles[i], dx, dy, dz);
        }
        return newTriangles;
    }

    /**
     * 三角形体积计算=三角形面积矢量的x分量*质心的x值
     *
     * @param tri 三角形
     * @return 三角形面积乘以质心
     */
    private static double volumeUnder(TriGeom tri) {
        Point[] p = tri.points();
        //三角形质心        
        double xc = (p[0].x + p[1].x + p[2].x) / 3.0;
        Vector vab = new Vector(p[0], p[1]);
        Vector vac = new Vector(p[0], p[2]);
        // x-component of cross product叉乘后的x分量
        double ax
                = vab.comtValue(1) * vac.comtValue(2)
                - vac.comtValue(1) * vab.comtValue(2);
        //矩形的面积x分量大小乘以0.5
        return ax * xc * 0.5;
    }

    private static Vector scaledCentroidUnder(TriGeom tri) {
        Point[] p = tri.points();
        Point p0 = p[0];
        Point p1 = p[1];
        Point p2 = p[2];

        Vector v01 = new Vector(p0, p1);
        Vector v02 = new Vector(p0, p2);

        Vector areaVector = v01.cross(v02);

        double xc
                = areaVector.comtValue(0) * (p0.x * p0.x + p1.x * p1.x + p2.x * p2.x
                + p0.x * p1.x + p0.x * p2.x + p1.x * p2.x);
        double yc
                = areaVector.comtValue(1) * (p0.y * p0.y + p1.y * p1.y + p2.y * p2.y
                + p0.y * p1.y + p0.y * p2.y + p1.y * p2.y);
        double zc
                = areaVector.comtValue(2) * (p0.z * p0.z + p1.z * p1.z + p2.z * p2.z
                + p0.z * p1.z + p0.z * p2.z + p1.z * p2.z);

        return new Vector(xc, yc, zc);
    }

    /**
     *
     * @param triangles 三角形数组
     * @return 体积大小
     */
    public static double signedVolume(TriGeom[] triangles) {
        double volume = 0.0;
        for (TriGeom t : triangles) {
            volume += volumeUnder(t);
        }
        return volume;
    }

    /**
     * Calculates the volume of a solid formed by a set of triangles. The
     * triangles must have points either all-clockwise or all-anti-clockwise,
     * looking from outside of the solid.
     * 计算由一组三角形组成的固体的体积。从实体外部看，三角形必须有全顺时针或全逆时针的点。
     *
     * @param triangles array of triangles 三角形的数组
     * @return non-negative volume enclosed by the set of triangles
     * 由一组三角形包围的非负体积
     */
    public static double volume(TriGeom[] triangles) {
        /*得到最大最小值的中点值*/
        Point translateBy = boundingBox(triangles).midPoint();
        /*为了进行数值计算，建议将实体平移到接近原点的位置，以减少由于大数平方引起的数值误差。*/
        TriGeom[] newTriangles = translateTriangles(triangles,
                -translateBy.x, -translateBy.y, -translateBy.z);

        return Math.abs(signedVolume(newTriangles));
    }

    /**
     * Calculates the centroid point of a solid formed by a set of triangles.
     * The triangles must have points either all-clockwise or
     * all-anti-clockwise, looking from outside of the solid. 计算由一组三角形形成的实体的质心点。
     * 从实体的外部看，三角形必须具有顺时针或逆时针的点。
     *
     * @param triangles array of triangles 三角形数组
     * @return centroid point of solid region enclosed by the set of triangles
     * 三角形集合包围的实体区域的质心点
     */
    public static Point centroid(TriGeom[] triangles) {
        Point translateBy = boundingBox(triangles).midPoint();
        /*为了进行数值计算，建议将实体平移到接近原点的位置，以减少由于大数平方引起的数值误差。*/
        TriGeom[] newTriangles = translateTriangles(triangles,
                -translateBy.x, -translateBy.y, -translateBy.z);

        double vol = signedVolume(newTriangles);
        Vector centroidPos = new Vector(0, 0, 0);

        for (TriGeom t : newTriangles) {
            centroidPos = centroidPos.add(scaledCentroidUnder(t));
        }
        centroidPos = centroidPos.mult(1.0 / vol / 24.0);
        /*平移回去了*/
        centroidPos = centroidPos.add(
                new Vector(translateBy.x, translateBy.y, translateBy.z)
        );
        return new Point(
                centroidPos.comtValue(0),
                centroidPos.comtValue(1),
                centroidPos.comtValue(2)
        );
    }
}
