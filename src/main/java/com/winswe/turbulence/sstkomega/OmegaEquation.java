/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.turbulence.sstkomega;

import com.alibaba.fastjson.JSONObject;
import static com.cup.util.Flux.faceValue;
import com.winswe.boundary.BoundaryCondition;
import com.winswe.field.VolScalarField;
import com.winswe.io.IOobject;
import com.winswe.matrix.Matrix;
import com.winswe.matrix.solve.MatrixSolver;
import com.winswe.matrix.solve.SolverPerformance;
import com.winswe.mesh.Structed2D;
import com.winswe.turbulence.Turbulence;
import com.winswe.util.Label;
import static java.lang.Math.max;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月25日 下午9:01:38
 */
public class OmegaEquation {

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

    public OmegaEquation(
            Structed2D mesh,
            VolScalarField k,
            VolScalarField epsilon,
            VolScalarField U,
            VolScalarField mu,
            VolScalarField mut,
            VolScalarField density,
            IOobject ioObject,
            VolScalarField S
    ) {
        this.mesh = mesh;
        this.k = k;
        this.omega = epsilon;
        this.velocity = U;
        this.mu = mu;
        this.mut = mut;
        this.rho = density;
        this.ioObject = ioObject;
        this.coe = new Matrix(mesh.getNX(), mesh.getNY());
        this.S = S;
        //        
        jSONObject = ioObject.
                getJsonObject().
                getJSONObject("Omega Equation");
        this.URF = jSONObject.getDoubleValue("Relax Factor");

        solverPerformance = new SolverPerformance(
                "Omega Equation",
                jSONObject.getIntValue("Sub Loop Number"),
                jSONObject.getDoubleValue("Convergence Criterion")
        );

        matrixSolver = MatrixSolver.factory(
                jSONObject.getString("Solve Name"),
                coe,
                epsilon.getFI(),
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

                mup = (mu.getFI()[IJ] + mut.getFI()[IJ] / KOmegaModel.sigmaOmega);
                muw = (mu.getFI()[IJW] + mut.getFI()[IJW] / KOmegaModel.sigmaOmega);
                mue = (mu.getFI()[IJE] + mut.getFI()[IJE] / KOmegaModel.sigmaOmega);
                mus = (mu.getFI()[IJS] + mut.getFI()[IJS] / KOmegaModel.sigmaOmega);
                mun = (mu.getFI()[IJN] + mut.getFI()[IJN] / KOmegaModel.sigmaOmega);

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
                Sp
                        = -KOmegaModel.beta * rho.getFI()[IJ]
                        * omega.getFI()[IJ];

                Spad
                        = BoundaryCondition.spadWest(
                                Dw, Fw, mesh.getDXP()[X - 1], volume,
                                omega.getBoundaryConditionParameter()[IJW])
                        + BoundaryCondition.spadEast(
                                De, Fe, mesh.getDXP()[X], volume,
                                omega.getBoundaryConditionParameter()[IJE])
                        + BoundaryCondition.spadSouth(
                                Ds, Fs, mesh.getDYP()[Y - 1], volume,
                                omega.getBoundaryConditionParameter()[IJS])
                        + BoundaryCondition.spadNorth(
                                Dn, Fn, mesh.getDYP()[Y], volume,
                                omega.getBoundaryConditionParameter()[IJN]);

                //get the coefficient of ap
                coe.getAP()[IJ]
                        = -(coe.getAW()[IJ] + coe.getAE()[IJ]
                        + coe.getAS()[IJ] + coe.getAN()[IJ]
                        + (Fe - Fw) + (Fn - Fs) - Spad * volume)
                        - Sp * volume;

                Sc
                        = KOmegaModel.alpha * (mut.getFI()[IJ] * S.getFI()[IJ])
                        * omega.getFI()[IJ] / k.getFI()[IJ];

                Scad
                        = BoundaryCondition.scadWest(
                                Dw, Fw, mesh.getDXP()[X - 1], volume,
                                omega.getBoundaryConditionParameter()[IJW])
                        + BoundaryCondition.scadEast(
                                De, Fe, mesh.getDXP()[X], volume,
                                omega.getBoundaryConditionParameter()[IJE])
                        + BoundaryCondition.scadSouth(
                                Ds, Fs, mesh.getDYP()[Y - 1], volume,
                                omega.getBoundaryConditionParameter()[IJS])
                        + BoundaryCondition.scadNorth(
                                Dn, Fn, mesh.getDYP()[Y], volume,
                                omega.getBoundaryConditionParameter()[IJN]);

                coe.getS()[IJ] = (Sc + Scad) * volume;
            }
        }

        double distance = 1E-30;
        for (int Y = 1; Y <= mesh.getNY(); ++Y) {
            for (int X = 1; X <= mesh.getNX(); ++X) {

                flag.setFlag(X, Y);

                //INDEX 
                IJ = mesh.getCellIndex(X, Y);

                if (flag.atBoundary()) {
                    if (flag.getWest() == 1) {
                        distance = Turbulence.distanceBetweenTwoPoints(
                                mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.X(mesh.getPointX1()[X - 1], mesh.getPointX2()[Y]),
                                mesh.Y(mesh.getPointX1()[X - 1], mesh.getPointX2()[Y]));

                    }
                    if (flag.getEast() == 1) {
                        distance = Turbulence.distanceBetweenTwoPoints(
                                mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.X(mesh.getPointX1()[X + 1], mesh.getPointX2()[Y]),
                                mesh.Y(mesh.getPointX1()[X + 1], mesh.getPointX2()[Y]));

                    }
                    if (flag.getSouth() == 1) {
                        distance = Turbulence.distanceBetweenTwoPoints(
                                mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y - 1]),
                                mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y - 1]));

                    }
                    if (flag.getNorth() == 1) {
                        distance = Turbulence.distanceBetweenTwoPoints(
                                mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                                mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y + 1]),
                                mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y + 1]));

                    }

                    double value
                            = Math.sqrt(max(0, k.getFI()[IJ]))
                            / (KOmegaModel.Cu25 * Turbulence.CAPPA * distance);
//ED(IJP)=SQRT(MAX(ZERO,TE(IJP)))/(CMU25*CAPPA*DN(IW))
//6.0 * mu[X][Y][Z] / (rho[X][Y][Z] * beta * DN * DN);
                    coe.getAP()[IJ] = 1.0;
                    coe.getAW()[IJ] = 0;//w
                    coe.getAE()[IJ] = 0;//e
                    coe.getAS()[IJ] = 0;//s
                    coe.getAN()[IJ] = 0;//n
                    coe.getS()[IJ] = value;
                }

                coe.getS()[IJ]
                        = coe.getS()[IJ]
                        + (URFFI - 1.0) * coe.getAP()[IJ]
                        * omega.getFI()[IJ];
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
