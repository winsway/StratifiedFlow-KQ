/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.boundary;

import static java.lang.Math.max;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月5日 下午9:16:53
 */
public interface BoundaryCondition1 {

    /**
     *
     * @param Type boundary condition type
     * @param Dw diffusion coefficient of face on boundary
     * @param Fw face flux on boundary
     * @param dxw distance from face centroid to inter point P
     * @param volume volume
     * @param bc boundary condition parameter
     * @return boundary value
     */
    static public double spadWest(
            int Type,
            double Dw,
            double Fw,
            double dxw,
            double volume,
            RobinBC bc) {
        if (bc == null || Type == 0 || Type == 7) {
            return 0;
        } else {
            double aw = Dw + max(Fw, 0);
            double a, b;
            a = bc.getA();
            b = bc.getB();
            double temp = -aw * (b * dxw) / (a + b * dxw) / volume;
            return temp;
        }
    }

    /**
     *
     * @param Type boundary condition type
     * @param Ds diffusion coefficient of face on boundary
     * @param Fs face flux on boundary
     * @param dxs distance from face centroid to inter point P
     * @param volume volume
     * @param bc boundary condition parameter
     * @return boundary value
     */
    static public double spadSouth(
            int Type,
            double Ds,
            double Fs,
            double dxs,
            double volume,
            RobinBC bc) {
        if (bc == null || Type == 0 || Type == 7) {
            return 0;
        } else {
            double aw = Ds + max(Fs, 0);
            double a, b;
            a = bc.getA();
            b = bc.getB();
            double temp = -aw * (b * dxs) / (a + b * dxs) / volume;
            return temp;
        }
    }

    /**
     *
     * @param Type boundary condition type
     * @param De diffusion coefficient of face on boundary
     * @param Fe face flux on boundary
     * @param dxe distance from face centroid to inter point P
     * @param volume volume
     * @param bc boundary condition parameter
     * @return boundary value
     */
    static public double SpadEast(
            int Type,
            double De,
            double Fe,
            double dxe,
            double volume,
            RobinBC bc) {
        if (bc == null || Type == 0 || Type == 7) {
            return 0;
        } else {
            double ae = De + max(-Fe, 0);
            double a, b;
            a = bc.getA();
            b = bc.getB();
            double temp = -ae * (b * dxe) / (a + b * dxe) / volume;
            return temp;
        }
    }

    /**
     *
     * @param Type boundary condition type
     * @param Dn diffusion coefficient of face on boundary
     * @param Fn face flux on boundary
     * @param dxn distance from face centroid to inter point P
     * @param volume volume
     * @param bc boundary condition parameter
     * @return boundary value
     */
    static public double SpadNorth(
            int Type,
            double Dn,
            double Fn,
            double dxn,
            double volume,
            RobinBC bc) {
        if (bc == null || Type == 0 || Type == 7) {
            return 0;
        } else {
            double ae = Dn + max(-Fn, 0);
            double a, b;
            a = bc.getA();
            b = bc.getB();
            double temp = -ae * (b * dxn) / (a + b * dxn) / volume;
            return temp;
        }
    }

    /**
     *
     * @param Type boundary condition type
     * @param Dw diffusion coefficient of face on boundary
     * @param Fw flux of face on boundary
     * @param dxw distance from face centroid to inter point P
     * @param volume volume
     * @param bc boundary condition
     * @return
     */
    static public double scadWest(
            int Type,
            double Dw,
            double Fw,
            double dxw,
            double volume,
            RobinBC bc
    ) {
        if (bc == null || Type == 0 || Type == 7) {
            return 0;
        } else {
            double aw = Dw + max(Fw, 0);
            double a, b, M;
            a = bc.getA();
            b = bc.getB();
            M = bc.getM();
            double temp = aw * (M * dxw) / (a + b * dxw) / volume;
            return temp;
        }
    }

    /**
     *
     * @param Type boundary condition type
     * @param Ds diffusion coefficient of face on boundary
     * @param Fs flux of face on boundary
     * @param dxs distance from face centroid to inter point P
     * @param volume volume
     * @param bc boundary condition
     * @return
     */
    static public double scadSouth(
            int Type,
            double Ds,
            double Fs,
            double dxs,
            double volume,
            RobinBC bc
    ) {
        if (bc == null || Type == 0 || Type == 7) {
            return 0;
        } else {
            double aw = Ds + max(Fs, 0);
            double a, b, M;
            a = bc.getA();
            b = bc.getB();
            M = bc.getM();
            double temp = aw * (M * dxs) / (a + b * dxs) / volume;
            return temp;
        }
    }

    /**
     * @param Type boundary condition type
     * @param De diffusion coefficient of face on boundary
     * @param Fe face flux on boundary
     * @param dxe distance from face centroid to inter point P
     * @param volume volume
     * @param bc boundary condition parameter
     * @return boundary value
     *
     */
    static public double scadEast(
            int Type,
            double De,
            double Fe,
            double dxe,
            double volume,
            RobinBC bc
    ) {
        if (bc == null || Type == 0 || Type == 7) {
            return 0;
        } else {
            double ae = De + max(-Fe, 0);
            double a, b, M;
            a = bc.getA();
            b = bc.getB();
            M = bc.getM();
            double temp = ae * (M * dxe) / (a + b * dxe) / volume;
            return temp;
        }
    }

    /**
     * @param Type boundary condition type
     * @param Dn diffusion coefficient of face on boundary
     * @param Fn face flux on boundary
     * @param dxn distance from face centroid to inter point P
     * @param volume volume
     * @param bc boundary condition parameter
     * @return boundary value
     *
     */
    static public double scadNorth(
            int Type,
            double Dn,
            double Fn,
            double dxn,
            double volume,
            RobinBC bc
    ) {
        if (bc == null || Type == 0 || Type == 7) {
            return 0;
        } else {
            double ae = Dn + max(-Fn, 0);
            double a, b, M;
            a = bc.getA();
            b = bc.getB();
            M = bc.getM();
            double temp = ae * (M * dxn) / (a + b * dxn) / volume;
            return temp;
        }
    }
}
