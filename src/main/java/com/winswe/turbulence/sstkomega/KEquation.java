/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.turbulence.sstkomega;

import com.alibaba.fastjson.JSONObject;
import com.winswe.turbulence.Turbulence;
import static com.cup.util.Flux.faceValue;
import com.winswe.boundary.BoundaryCondition;
import com.winswe.field.VolScalarField;
import com.winswe.io.IOobject;
import com.winswe.matrix.Matrix;
import com.winswe.matrix.solve.MatrixSolver;
import com.winswe.matrix.solve.SolverPerformance;
import com.winswe.mesh.Structed2D;
import com.winswe.util.Label;
import static java.lang.Math.max;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月25日 下午9:01:19
 */
public class KEquation {

    private final Structed2D mesh;
    private final VolScalarField k;
    private final VolScalarField omega;
    private final VolScalarField velocity;
    private final VolScalarField mu;
    private final VolScalarField mut;
    private final VolScalarField rho;
    private final IOobject ioObject;
    private final double URF;
    private final MatrixSolver matrixSolver;
    private final SolverPerformance solverPerformance;
    private final JSONObject jSONObject;
    private final Matrix coe;
    private final VolScalarField S;
    private final VolScalarField F1SST;
    private final VolScalarField GEN;

    public KEquation(
            Structed2D mesh,
            VolScalarField k,
            VolScalarField omega,
            VolScalarField velocity,
            VolScalarField mu,
            VolScalarField mut,
            VolScalarField density,
            IOobject ioObject,
            VolScalarField S,
            VolScalarField F1SST,
            VolScalarField GEN
    ) {
        this.mesh = mesh;
        this.k = k;
        this.omega = omega;
        this.velocity = velocity;
        this.mu = mu;
        this.mut = mut;
        this.rho = density;
        this.ioObject = ioObject;
        this.coe = new Matrix(mesh.getNX(), mesh.getNY());
        this.S = S;
        this.F1SST = F1SST;
        this.GEN = GEN;
        //        
        jSONObject = ioObject.
                getJsonObject().
                getJSONObject("K Equation");
        this.URF = jSONObject.getDoubleValue("Relax Factor");

        solverPerformance = new SolverPerformance(
                "Kequation",
                jSONObject.getIntValue("Sub Loop Number"),
                jSONObject.getDoubleValue("Convergence Criterion")
        );

        matrixSolver = MatrixSolver.factory(
                jSONObject.getString("Solve Name"),
                coe,
                k.getFI(),
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
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;

        double sigma = 0;

        Label flag = new Label(mesh.getNX(), mesh.getNY());

        int IJW, IJE, IJS, IJN, IJ;

        coe.intialMatrix();

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

                sigma = SSTKOmegaModel.getPhi(
                        SSTKOmegaModel.sigmak1,
                        SSTKOmegaModel.sigmak2,
                        F1SST.getFI()[IJ]
                );

                mup = (mu.getFI()[IJ] + mut.getFI()[IJ] * sigma);
                muw = (mu.getFI()[IJW] + mut.getFI()[IJW] * sigma);
                mue = (mu.getFI()[IJE] + mut.getFI()[IJE] * sigma);
                mus = (mu.getFI()[IJS] + mut.getFI()[IJS] * sigma);
                mun = (mu.getFI()[IJN] + mut.getFI()[IJN] * sigma);

                gamw
                        = faceValue(
                                muw, mup - muw,
                                mesh.getPointX1()[X - 1], mesh.getLineX1()[X - 1],
                                mesh.getPointX1()[X]
                        );
                game
                        = faceValue(
                                mup, mue - mup,
                                mesh.getPointX1()[X], mesh.getLineX1()[X],
                                mesh.getPointX1()[X + 1]
                        );
                gams
                        = faceValue(
                                mus, mup - mus,
                                mesh.getPointX2()[Y - 1], mesh.getLineX2()[Y - 1],
                                mesh.getPointX2()[Y]
                        );
                gamn
                        = faceValue(
                                mup, mun - mup,
                                mesh.getPointX2()[Y], mesh.getLineX2()[Y],
                                mesh.getPointX2()[Y + 1]
                        );

                Dw = -gamw * mesh.getDYV()[Y] / mesh.getDXP()[X - 1];
                De = -game * mesh.getDYV()[Y] / mesh.getDXP()[X];
                Ds = -gams * mesh.getDXU()[X] / mesh.getDYP()[Y - 1];
                Dn = -gamn * mesh.getDXU()[X] / mesh.getDYP()[Y];

                coe.getAW()[IJ] = (1 - flag.getWest()) * (Math.max(Fw, 0.0) + Dw);
                coe.getAE()[IJ] = (1 - flag.getEast()) * (Math.max(-Fe, 0.0) + De);
                coe.getAS()[IJ] = (1 - flag.getSouth()) * (Math.max(Fs, 0.0) + Ds);
                coe.getAN()[IJ] = (1 - flag.getNorth()) * (Math.max(-Fn, 0.0) + Dn);

                //setting add source.
                Sp = -SSTKOmegaModel.betaStar * rho.getFI()[IJ] * omega.getFI()[IJ];
                Spad
                        = BoundaryCondition.spadWest(
                                Dw, Fw, mesh.getDXP()[X - 1], volume,
                                k.getBoundaryConditionParameter()[IJW])
                        + BoundaryCondition.spadEast(
                                De, Fe, mesh.getDXP()[X], volume,
                                k.getBoundaryConditionParameter()[IJE])
                        + BoundaryCondition.spadSouth(
                                Ds, Fs, mesh.getDYP()[Y - 1], volume,
                                k.getBoundaryConditionParameter()[IJS])
                        + BoundaryCondition.spadNorth(
                                Dn, Fn, mesh.getDYP()[Y], volume,
                                k.getBoundaryConditionParameter()[IJN]);

                //get the coefficient of ap
                coe.getAP()[IJ]
                        = -(coe.getAW()[IJ] + coe.getAE()[IJ]
                        + coe.getAS()[IJ] + coe.getAN()[IJ]
                        + (Fe - Fw) + (Fn - Fs) - Spad * volume)
                        - Sp * volume;

                //<editor-fold>
                if (flag.atBoundary()) {
                    double viss = mu.getFI()[IJ];
                    double yPlus, distance = 0;
                    double tau = 0;
                    double CK = SSTKOmegaModel.getCK(k.getFI()[IJ]);
                    double viscWall, visW;
                    if (flag.getWest() == 1) {
                        distance = Turbulence.distanceBetweenTwoPoints(
                                mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.X(mesh.getPointX1()[X - 1], mesh.getPointX2()[Y]),
                                mesh.Y(mesh.getPointX1()[X - 1], mesh.getPointX2()[Y]));
                        yPlus = Turbulence.yPlusSubLay(rho.getFI()[IJ], mu.getFI()[IJ], CK, distance);
                        viscWall = Turbulence.getViscosityWall(yPlus, mu.getFI()[IJ]);
                        visW = max(mu.getFI()[IJ], viscWall);
                        if (yPlus >= SSTKOmegaModel.CTRANS) {
                            viss = visW;
                        }
                        tau
                                = Turbulence.getTauWall(
                                        viss,
                                        velocity.getFI()[IJ],
                                        velocity.getFI()[IJW],
                                        distance);
                    }
                    if (flag.getEast() == 1) {
                        distance = Turbulence.distanceBetweenTwoPoints(
                                mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.X(mesh.getPointX1()[X + 1], mesh.getPointX2()[Y]),
                                mesh.Y(mesh.getPointX1()[X + 1], mesh.getPointX2()[Y]));
                        yPlus = Turbulence.yPlusSubLay(rho.getFI()[IJ], mu.getFI()[IJ], CK, distance);
                        viscWall = Turbulence.getViscosityWall(yPlus, mu.getFI()[IJ]);
                        visW = max(mu.getFI()[IJ], viscWall);
                        if (yPlus >= SSTKOmegaModel.CTRANS) {
                            viss = visW;
                        }
                        tau
                                = Turbulence.getTauWall(
                                        viss,
                                        velocity.getFI()[IJ],
                                        velocity.getFI()[IJE],
                                        distance);
                    }
                    if (flag.getSouth() == 1) {
                        distance = Turbulence.distanceBetweenTwoPoints(
                                mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y - 1]),
                                mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y - 1]));
                        yPlus = Turbulence.yPlusSubLay(rho.getFI()[IJ], mu.getFI()[IJ], CK, distance);
                        viscWall = Turbulence.getViscosityWall(yPlus, mu.getFI()[IJ]);
                        visW = max(mu.getFI()[IJ], viscWall);
                        if (yPlus >= SSTKOmegaModel.CTRANS) {
                            viss = visW;
                        }
                        tau
                                = Turbulence.getTauWall(
                                        viss,
                                        velocity.getFI()[IJ],
                                        velocity.getFI()[IJS],
                                        distance);
                    }
                    if (flag.getNorth() == 1) {
                        distance = Turbulence.distanceBetweenTwoPoints(
                                mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y + 1]),
                                mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y + 1]));
                        yPlus = Turbulence.yPlusSubLay(rho.getFI()[IJ], mu.getFI()[IJ], CK, distance);
                        viscWall = Turbulence.getViscosityWall(yPlus, mu.getFI()[IJ]);
                        visW = max(mu.getFI()[IJ], viscWall);
                        if (yPlus >= SSTKOmegaModel.CTRANS) {
                            viss = visW;
                        }
                        tau = Turbulence.getTauWall(
                                viss,
                                velocity.getFI()[IJ],
                                velocity.getFI()[IJN],
                                distance);
                    }
                    GEN.getFI()[IJ] = Turbulence.getProduceTerm(tau, k.getFI()[IJ], distance);
                } else {
                    GEN.getFI()[IJ] = Math.min(
                            mut.getFI()[IJ] * S.getFI()[IJ],
                            10.0
                            * SSTKOmegaModel.betaStar
                            * rho.getFI()[IJ]
                            * omega.getFI()[IJ]
                            * k.getFI()[IJ]
                    );
                }
                //</editor-fold>

                Sc = GEN.getFI()[IJ];

                Scad
                        = BoundaryCondition.scadWest(
                                Dw, Fw, mesh.getDXP()[X - 1], volume,
                                k.getBoundaryConditionParameter()[IJW])
                        + BoundaryCondition.scadEast(
                                De, Fe, mesh.getDXP()[X], volume,
                                k.getBoundaryConditionParameter()[IJE])
                        + BoundaryCondition.scadSouth(
                                Ds, Fs, mesh.getDYP()[Y - 1], volume,
                                k.getBoundaryConditionParameter()[IJS])
                        + BoundaryCondition.scadNorth(
                                Dn, Fn, mesh.getDYP()[Y], volume,
                                k.getBoundaryConditionParameter()[IJN]);

                coe.getS()[IJ] = (Sc + Scad) * volume;
//                System.out.println(" coefficietn debug");
                //under relax
                coe.getS()[IJ]
                        = coe.getS()[IJ]
                        + (URFFI - 1.0) * coe.getAP()[IJ]
                        * k.getFI()[IJ];
                coe.getAP()[IJ] = coe.getAP()[IJ] * URFFI;
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
