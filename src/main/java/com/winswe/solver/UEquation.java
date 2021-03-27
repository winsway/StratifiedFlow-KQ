/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.solver;

import com.alibaba.fastjson.JSONObject;
import static com.cup.util.Flux.faceValue;

import com.winswe.boundary.BoundaryCondition;
import com.winswe.field.VolScalarField;
import com.winswe.io.IOobject;
import com.winswe.matrix.Matrix;
import com.winswe.matrix.solve.MatrixSolver;
import com.winswe.matrix.solve.SolverPerformance;
import com.winswe.mesh.Structed2D;
import com.winswe.util.Label;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月1日 上午4:45:34
 */
public class UEquation {

    private final double dpdz;
    private final Structed2D mesh;
    private final VolScalarField U;
    private final VolScalarField mu;
    private final VolScalarField rho;
    private final IOobject ioObject;
    private final double URF;
    private final MatrixSolver matrixSolver;
    private final SolverPerformance solverPerformance;
    private final JSONObject jSONObject;
    private final Matrix coeW;

    /**
     *
     * @param dpdz pressure gradient
     * @param mesh mesh
     * @param U velocity
     * @param mu effective dynamic viscosity
     * @param density density
     * @param ioObject io object
     */
    public UEquation(
            double dpdz,
            Structed2D mesh,
            VolScalarField U,
            VolScalarField mu,
            VolScalarField density,
            IOobject ioObject
    ) {
        this.dpdz = dpdz;
        this.mesh = mesh;
        this.U = U;
        this.mu = mu;
        this.rho = density;
        this.ioObject = ioObject;
        this.coeW = new Matrix(mesh.getNX(), mesh.getNY());
        //        
        jSONObject = ioObject.
                getJsonObject().
                getJSONObject("Uequation");
        this.URF = jSONObject.getDoubleValue("Relax Factor");

        solverPerformance = new SolverPerformance(
                "Uequation",
                jSONObject.getIntValue("Sub Loop Number"),
                jSONObject.getDoubleValue("Convergence Criterion")
        );

        matrixSolver = MatrixSolver.factory(
                jSONObject.getString("Solve Name"),
                coeW,
                U.getFI(),
                mesh,
                solverPerformance);

    }

    /**
     * discrete equation
     */
    public void discrete() {
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double gamw, game, gams, gamn;
        double Spad, Scad, Sp, Sc;
        Label flag = new Label(mesh.getNX(), mesh.getNY());
        int IJW, IJE, IJS, IJN, IJ;

        coeW.intialMatrix();

        for (int Y = 1; Y <= mesh.getNY(); ++Y) {
            for (int X = 1; X <= mesh.getNX(); ++X) {

                flag.setFlag(X, Y);
                double volume = mesh.getVolume(X, Y);

                //INDEX 
                IJ = mesh.getCellIndex(X, Y);
                IJW = mesh.getCellIndex(X - 1, Y);
                IJE = mesh.getCellIndex(X + 1, Y);
                IJS = mesh.getCellIndex(X, Y - 1);
                IJN = mesh.getCellIndex(X, Y + 1);

                Fw = 0;
                Fe = 0;
                Fs = 0;
                Fn = 0;

                gamw
                        = faceValue(
                                mu.getFI()[IJW], mu.getFI()[IJ] - mu.getFI()[IJW],
                                mesh.getPointX1()[X - 1], mesh.getLineX1()[X - 1],
                                mesh.getPointX1()[X]
                        );
                game
                        = faceValue(
                                mu.getFI()[IJ], mu.getFI()[IJE] - mu.getFI()[IJ],
                                mesh.getPointX1()[X], mesh.getLineX1()[X],
                                mesh.getPointX1()[X + 1]
                        );
                gams
                        = faceValue(
                                mu.getFI()[IJS], mu.getFI()[IJ] - mu.getFI()[IJS],
                                mesh.getPointX2()[Y - 1], mesh.getLineX2()[Y - 1],
                                mesh.getPointX2()[Y]
                        );
                gamn
                        = faceValue(
                                mu.getFI()[IJ], mu.getFI()[IJN] - mu.getFI()[IJ],
                                mesh.getPointX2()[Y], mesh.getLineX2()[Y],
                                mesh.getPointX2()[Y + 1]
                        );

                Dw = -gamw * mesh.getDYV()[Y] / mesh.getDXP()[X - 1];
                De = -game * mesh.getDYV()[Y] / mesh.getDXP()[X];
                Ds = -gams * mesh.getDXU()[X] / mesh.getDYP()[Y - 1];
                Dn = -gamn * mesh.getDXU()[X] / mesh.getDYP()[Y];

                coeW.getAW()[IJ] = (1 - flag.getWest()) * (Math.max(Fw, 0.0) + Dw);
                coeW.getAE()[IJ] = (1 - flag.getEast()) * (Math.max(-Fe, 0.0) + De);
                coeW.getAS()[IJ] = (1 - flag.getSouth()) * (Math.max(Fs, 0.0) + Ds);
                coeW.getAN()[IJ] = (1 - flag.getNorth()) * (Math.max(-Fn, 0.0) + Dn);

                //setting add source.
                Sp = 0;
                Spad
                        = BoundaryCondition.spadWest(
                                Dw, Fw, mesh.getDXP()[X - 1], volume,
                                U.getBoundaryConditionParameter()[IJW])
                        + BoundaryCondition.SpadEast(
                                De, Fe, mesh.getDXP()[X], volume,
                                U.getBoundaryConditionParameter()[IJE])
                        + BoundaryCondition.spadSouth(
                                Ds, Fs, mesh.getDYP()[Y - 1], volume,
                                U.getBoundaryConditionParameter()[IJS])
                        + BoundaryCondition.spadSouth(
                                Dn, Fn, mesh.getDYP()[Y], volume,
                                U.getBoundaryConditionParameter()[IJN]);

                //get the coefficient of ap
                coeW.getAP()[IJ]
                        = -(coeW.getAW()[IJ] + coeW.getAE()[IJ]
                        + coeW.getAS()[IJ] + coeW.getAN()[IJ]
                        + (Fe - Fw) + (Fn - Fs)
                        - (Sp + Spad) * volume);

                Sc = dpdz;

                Scad
                        = BoundaryCondition.scadWest(
                                Dw, Fw, mesh.getDXP()[X - 1], volume,
                                U.getBoundaryConditionParameter()[IJW])
                        + BoundaryCondition.scadEast(
                                De, Fe, mesh.getDXP()[X], volume,
                                U.getBoundaryConditionParameter()[IJE])
                        + BoundaryCondition.scadSouth(
                                Ds, Fs, mesh.getDYP()[Y - 1], volume,
                                U.getBoundaryConditionParameter()[IJS])
                        + BoundaryCondition.scadNorth(
                                Dn, Fn, mesh.getDYP()[Y], volume,
                                U.getBoundaryConditionParameter()[IJN]);

                coeW.getS()[IJ] = (Sc + Scad) * volume;
//                System.out.println(" coefficietn debug");
                //under relax
                coeW.getS()[IJ]
                        = coeW.getS()[IJ]
                        + (URFFI - 1.0) * coeW.getAP()[IJ]
                        * U.getFI()[IJ];
                coeW.getAP()[IJ] = coeW.getAP()[IJ] * URFFI;
            }
        }
    }

    /**
     * solve equation
     */
    public void solve() {
        matrixSolver.solve();
    }

    public SolverPerformance getSolverPerformance() {
        return solverPerformance;
    }

}
