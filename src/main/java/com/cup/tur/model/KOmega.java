/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.tur.model;

import static com.cup.boundary.BoundaryCondition.Scade;
import static com.cup.boundary.BoundaryCondition.Scadw;
import static com.cup.boundary.BoundaryCondition.Spade;
import static com.cup.boundary.BoundaryCondition.Spadw;
import com.cup.boundary.factory.Label;
import com.cup.boundary.factory.RobinBC;
import com.cup.field.Coefficient;
import com.cup.field.Field;
import com.cup.field.ProF;
import com.cup.field.Scalarfield;
import com.winswe.io.Writer;
import static com.cup.log.Residual.getRes;
import com.cup.mesh.Mesh;
import com.cup.mesh.factory.BiSameV2;
import com.cup.system.ControlDict;
import com.cup.system.FvSolution;
import com.cup.system.SystemControl;
import com.cup.turbulence.Turbulence;
import static com.cup.turbulence.Turbulence.getTAU;
import static com.cup.turbulence.Turbulence.setK;
import static com.cup.turbulence.Turbulence.setOmega;
import static com.cup.util.Flux.faceValue;
import com.cup.util.MatrixSolver;
import com.cup.util.Tool;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.max;

/**
 *
 * @author winsway
 */
public class KOmega {

    FvSolution fvsolution = new FvSolution();
    ControlDict controldict = new ControlDict();
    SystemControl sys;
//    
    ProF[] pro = new ProF[6];
    int S = 0, rho = 1, mu = 2, mut = 3, lambda = 4, cp = 5;
//    
    Field[] F = new Field[4];
    int W = 0, K = 1, Omega = 2, T = 3;
//    
    Mesh mesh;
    BiSameV2 bm;
//    
    double[] SOR = {0.95, 0.9, 0.9, 0.9};
    /**
     * 需要引入边界条件的基础参数
     */
    RobinBC[] bc;
    /**
     * 物性参数
     */
    double[] testO, testW;
//    
    PrintWriter fileWriter;
//
    double dpdz, hl, R;
//
    double[][][] phi;

    public KOmega(String position, String Title) {
        this.position = position;
        this.Title = Title;
    }

    public void setBoundary(RobinBC[] bc) {
        this.bc = bc;
    }

    public void setFluid(double[] testO, double[] testW) {
        this.testO = testO;
        this.testW = testW;
    }

    public void application(double dpdz, double hl, double R, int I) throws FileNotFoundException {
        this.dpdz = dpdz;
        this.hl = hl;
        this.R = R;
        createGrid();
        createFields();
        initField();
        setBoundary();
        sys = new SystemControl(mesh, controldict, fvsolution);
        Solver(I);
        System.out.println("求解完成！");
    }

    public void createGrid() {
        bm = new BiSameV2(12, Math.PI, 50, 50);
        bm.radia = R;
        double D = bm.radia * 2.0;
        hl = D * hl;
        double gamma = acos(1.0 - 2 * hl / D);
        bm.theta1 = gamma;
        bm.theta2 = PI;
        mesh = bm;
        mesh.blockMesh();
    }

    void createFields() {
        for (int i = 0; i < F.length; i++) {
            F[i] = new Scalarfield(bm.numPx(), bm.numPy(), bm.numPz());
        }
        for (int i = 0; i < pro.length; i++) {
            pro[i] = new ProF(bm.numPx(), bm.numPy(), bm.numPz());
        }
        phi = new double[bm.numPx()][bm.numPy()][bm.numPz()];
        System.out.println("场创建完成!");
    }

    void initField() {
        for (int K = 0; K <= bm.getNZ() + 1; K++) {
            for (int J = 0; J <= bm.getNY() + 1; J++) {
                for (int I = 0; I <= bm.getNX() + 1; I++) {
                    if (J <= bm.interhL) {
                        phi[I][J][K] = -1;
                        pro[rho].getF()[I][J][K] = testO[0];
                        pro[mu].getF()[I][J][K] = testO[1];
                        pro[lambda].getF()[I][J][K] = testO[2];
                        pro[cp].getF()[I][J][K] = testO[3];
                    } else {
                        phi[I][J][K] = 1;
                        pro[rho].getF()[I][J][K] = testW[0];
                        pro[mu].getF()[I][J][K] = testW[1];
                        pro[lambda].getF()[I][J][K] = testW[2];
                        pro[cp].getF()[I][J][K] = testW[3];
                    }
                }
            }
        }
    }

    void setBoundary() {
        int xw = 0, xe = bm.getNX() + 1;
        int yw = 0, ye = bm.getNY() + 1;
        //<editor-fold>
        //X边界
        for (int Z = 0; Z <= bm.getNZ() + 1; Z++) {
            for (int Y = 0; Y <= bm.getNY() + 1; Y++) {
                F[W].getBoundary().getBoundaryRegion().setType(xw, Y, Z, 1, bc[W]);
                F[K].getBoundary().getBoundaryRegion().setType(xw, Y, Z, 1, bc[K]);
                F[Omega].getBoundary().getBoundaryRegion().setType(xw, Y, Z, 1, bc[Omega]);
                F[T].getBoundary().getBoundaryRegion().setType(xw, Y, Z, 3, bc[T]);
                //xe
                F[W].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 1, bc[W]);
                F[K].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 1, bc[K]);
                F[Omega].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 1, bc[Omega]);
                F[T].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 3, bc[T]);
            }
        }
        //Y边界
        for (int Z = 0; Z <= bm.getNZ() + 1; Z++) {
            for (int X = 0; X <= bm.getNX() + 1; X++) {
                F[W].getBoundary().getBoundaryRegion().setType(X, yw, Z, 1, bc[W]);
                F[K].getBoundary().getBoundaryRegion().setType(X, yw, Z, 1, bc[K]);
                F[Omega].getBoundary().getBoundaryRegion().setType(X, yw, Z, 1, bc[Omega]);
                F[T].getBoundary().getBoundaryRegion().setType(X, yw, Z, 3, bc[T]);
                //ye
                F[W].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, bc[W]);
                F[K].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, bc[K]);
                F[Omega].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, bc[Omega]);
                F[T].getBoundary().getBoundaryRegion().setType(X, ye, Z, 3, bc[T]);
            }
        }
    }

    double[] RES = new double[4];

    void Solver(int I) throws FileNotFoundException {
        if (I == 1) {
            scheme1();
        }
    }

    double GravityZ = 0;

    void scheme1() throws FileNotFoundException {
        int count = 0;
        RESULT(count);
        System.out.println("计算速度场W");
        cacVelocity();
        System.out.println("初始化湍流参数！");
        setTurPar();
        do {
            count++;
            Turbulence.cacScenter(F[W], mesh, pro[S].getF());
            RES[W] = cacVelocity();
            cacTur();
            updateProF();
        } while (getMaxRES() > 1e-4 || count < 500);
        RESULT(count);
    }

    /**
     * 计算速度场
     *
     * @return 最大残差
     */
    double cacVelocity() {
        double resw;
        System.out.println("计算速度场W");
        cacW(
                dpdz,
                mesh, F[W],
                pro[mu].getF(), pro[rho].getF(), pro[mut].getF(), sys
        );
        resw = getRes(F[W].getNewField(), F[W].getOldField());
        F[W].copyNew2Old();
        return resw;
    }

    /**
     * 计算湍流场
     *
     * @return 最大残差
     */
    double cacTur() {
        System.out.println("计算湍流动能!");
        cacK(mesh, F[K],
                F[Omega].getNewField(), F[W].getNewField(),
                pro[mu].getF(), pro[rho].getF(), pro[mut].getF(), pro[S].getF(), sys);
        RES[K] = getRes(F[K].getNewField(), F[K].getOldField());
        F[K].copyNew2Old();
        System.out.println("计算湍流动能的比耗散率!");
        cacOmega(mesh, F[Omega],
                F[K].getNewField(), F[W].getNewField(),
                pro[mu].getF(), pro[rho].getF(), pro[mut].getF(), pro[S].getF(), sys);
        RES[Omega] = getRes(F[Omega].getNewField(), F[Omega].getOldField());
        F[Omega].copyNew2Old();
        return RES[Omega];
    }

    /**
     * 更新物性参数
     */
    void updateProF() {
        for (int Z = 1; Z <= mesh.getNZ(); ++Z) {
            for (int Y = 1; Y <= mesh.getNY(); ++Y) {
                for (int X = 1; X <= mesh.getNX(); ++X) {
                    pro[mut].getF()[X][Y][Z]
                            = cacMut(
                                    pro[rho].getF()[X][Y][Z],
                                    F[K].getNewField()[X][Y][Z],
                                    F[Omega].getNewField()[X][Y][Z]
                            );
                }
            }
        }

//        double value;
//        double CK, DN, YPL, VISCW, VISW;
//        //<editor-fold>
//        //X边界
//        for (int Z = 1; Z <= bm.getNZ(); Z++) {
//            for (int Y = 1; Y <= bm.getNY(); Y++) {
//                int X = 1;
//                DN = Tool.distance(
//                        mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y]),
//                        mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y]),
//                        mesh.realX(mesh.gettPx()[X - 1], mesh.gettPy()[Y]),
//                        mesh.realY(mesh.gettPx()[X - 1], mesh.gettPy()[Y]));
//                CK = Turbulence.getCK(F[K].getOldField()[X][Y][Z]);
//                YPL = Turbulence.getYPL(pro[rho].getF()[X][Y][Z], pro[mu].getF()[X][Y][Z], CK, DN);
//                VISCW = Turbulence.getVISCW(YPL, pro[mu].getF()[X][Y][Z]);
//                VISW = max(pro[mu].getF()[X][Y][Z], VISCW);
//                value = max(VISW - pro[mu].getF()[X][Y][Z], 0);
//                pro[mut].getF()[X][Y][Z] = value;
//                //xe
//                X = mesh.getNX();
//                DN = Tool.distance(
//                        mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y]),
//                        mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y]),
//                        mesh.realX(mesh.gettPx()[X + 1], mesh.gettPy()[Y]),
//                        mesh.realY(mesh.gettPx()[X + 1], mesh.gettPy()[Y]));
//                CK = Turbulence.getCK(F[K].getOldField()[X][Y][Z]);
//                YPL = Turbulence.getYPL(pro[rho].getF()[X][Y][Z], pro[mu].getF()[X][Y][Z], CK, DN);
//                VISCW = Turbulence.getVISCW(YPL, pro[mu].getF()[X][Y][Z]);
//                VISW = max(pro[mu].getF()[X][Y][Z], VISCW);
//                value = max(VISW - pro[mu].getF()[X][Y][Z], 0);
//                pro[mut].getF()[X][Y][Z] = value;
//            }
//        }
//        //Y边界
//        for (int Z = 1; Z <= bm.getNZ(); Z++) {
//            for (int X = 1; X <= bm.getNX(); X++) {
//                int Y = 1;
//                DN = Tool.distance(
//                        mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y]),
//                        mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y]),
//                        mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y - 1]),
//                        mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y - 1]));
//                CK = Turbulence.getCK(F[K].getOldField()[X][Y][Z]);
//                YPL = Turbulence.getYPL(pro[rho].getF()[X][Y][Z], pro[mu].getF()[X][Y][Z], CK, DN);
//                VISCW = Turbulence.getVISCW(YPL, pro[mu].getF()[X][Y][Z]);
//                VISW = max(pro[mu].getF()[X][Y][Z], VISCW);
//                value = max(VISW - pro[mu].getF()[X][Y][Z], 0);
//                pro[mut].getF()[X][Y][Z] = value;
//                //xe
//                Y = mesh.getNY();
//                DN = Tool.distance(
//                        mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y]),
//                        mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y]),
//                        mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y + 1]),
//                        mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y + 1]));
//                CK = Turbulence.getCK(F[K].getOldField()[X][Y][Z]);
//                YPL = Turbulence.getYPL(pro[rho].getF()[X][Y][Z], pro[mu].getF()[X][Y][Z], CK, DN);
//                VISCW = Turbulence.getVISCW(YPL, pro[mu].getF()[X][Y][Z]);
//                VISW = max(pro[mu].getF()[X][Y][Z], VISCW);
//                value = max(VISW - pro[mu].getF()[X][Y][Z], 0);
//                pro[mut].getF()[X][Y][Z] = value;
//            }
//        }
    }

    /**
     * 设置初始的湍流参数k,omega,mut
     */
    void setTurPar() {
        double value;
        for (int Z = 1; Z <= mesh.getNZ(); ++Z) {
            for (int Y = 1; Y <= mesh.getNY(); ++Y) {
                for (int X = 1; X <= mesh.getNX(); ++X) {
                    value = setK(F[W].getNewField()[X][Y][Z]);
                    F[K].setNewField(X, Y, Z, value);
//                    
                    value = setOmega(F[K].getNewField()[X][Y][Z], R);
                    F[Omega].setNewField(X, Y, Z, value);
//                    
                    pro[mut].getF()[X][Y][Z]
                            = cacMut(
                                    pro[rho].getF()[X][Y][Z],
                                    F[K].getNewField()[X][Y][Z],
                                    F[Omega].getNewField()[X][Y][Z]
                            );
                }
            }
        }
    }

    /**
     * 计算湍流粘度
     *
     * @param rho
     * @param K
     * @param omega
     * @return 湍流粘度
     */
    public double cacMut(double rho, double K, double omega) {
        double alpha = 1.0;
        double temp = alpha * rho * K / (omega + 1e-30);
        return temp;
    }

    /**
     *
     *
     * @return 得到最大残差
     */
    double getMaxRES() {
        double max = 0;
        System.out.println("RESW = " + RES[W]);
        System.out.println("RESK = " + RES[K]);
        System.out.println("RESOmega = " + RES[Omega]);
        System.out.println("REST = " + RES[T]);
        for (int i = 0; i < RES.length; i++) {
            max = RES[i] > max ? RES[i] : max;
        }
        return RES[W];
    }

    /**
     * 计算速度场
     *
     * @param dpdz
     */
    void cacW(double dpdz,
            Mesh mesh, Field F,
            double[][][] mu, double[][][] rho,
            double[][][] Mut, SystemControl sys
    ) {
        double URF = SOR[this.W];
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;
        double alpha = 1.0;
        Label flag = new Label("W");
        int i = mesh.numPx();
        int j = mesh.numPy();
        Coefficient coeW = new Coefficient(i, j);
        for (int Z = 1; Z <= mesh.getNZ(); ++Z) {
            for (int Y = 1; Y <= mesh.getNY(); ++Y) {
                for (int X = 1; X <= mesh.getNX(); ++X) {
                    //<editor-fold>
                    flag.setFlag(F.getBoundary(), X, Y);
                    double vol = mesh.getVol(X, Y, Z);
                    Fw = 0;
                    Fe = 0;
                    Fs = 0;
                    Fn = 0;
//
                    mup = (mu[X][Y][Z] + Mut[X][Y][Z] / alpha);
                    muw = (mu[X - 1][Y][Z] + Mut[X - 1][Y][Z] / alpha);
                    mue = (mu[X + 1][Y][Z] + Mut[X + 1][Y][Z] / alpha);
                    mus = (mu[X][Y - 1][Z] + Mut[X][Y - 1][Z] / alpha);
                    mun = (mu[X][Y + 1][Z] + Mut[X][Y + 1][Z] / alpha);
//得到gamma
                    gamw = faceValue(muw, mup - muw, mesh.gettPx()[X - 1], mesh.getTx()[X - 1], mesh.gettPx()[X]);
                    game = faceValue(mup, mue - mup, mesh.gettPx()[X], mesh.getTx()[X], mesh.gettPx()[X + 1]);
                    gams = faceValue(mus, mup - mus, mesh.gettPy()[Y - 1], mesh.getTy()[Y - 1], mesh.gettPy()[Y]);
                    gamn = faceValue(mup, mun - mup, mesh.gettPy()[Y], mesh.getTy()[Y], mesh.gettPy()[Y + 1]);
//GET d
                    Dw = gamw * mesh.getDYV()[Y] / mesh.getDXP()[X - 1];
                    De = game * mesh.getDYV()[Y] / mesh.getDXP()[X];
                    Ds = gams * mesh.getDXU()[X] / mesh.getDYP()[Y - 1];
                    Dn = gamn * mesh.getDXU()[X] / mesh.getDYP()[Y];
//系数
                    coeW.aw[X][Y][Z] = (1 - flag.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeW.ae[X][Y][Z] = (1 - flag.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeW.as[X][Y][Z] = (1 - flag.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeW.an[X][Y][Z] = (1 - flag.n) * (Math.max(-Fn, 0.0) + Dn);//n
//设定附加源项
                    Sp = 0;
                    Spad
                            = Spadw(F.getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, mesh.getDXP()[X - 1], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Spade(F.getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, mesh.getDXP()[X], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Spadw(F.getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, mesh.getDYP()[Y - 1], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Spade(F.getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, mesh.getDYP()[Y], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
//得到ap系数
                    coeW.ap[X][Y][Z]
                            = coeW.aw[X][Y][Z] + coeW.ae[X][Y][Z]
                            + coeW.as[X][Y][Z] + coeW.an[X][Y][Z]
                            + (Fe - Fw) + (Fn - Fs) - (Sp + Spad) * vol;
//得到源项
                    Sc = dpdz - rho[X][Y][Z] * GravityZ;
                    Scad
                            = Scadw(F.getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, mesh.getDXP()[X - 1], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Scade(F.getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, mesh.getDXP()[X], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Scadw(F.getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, mesh.getDYP()[Y - 1], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Scade(F.getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, mesh.getDYP()[Y], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
                    coeW.b[X][Y][Z] = (Sc + Scad) * vol;
                    //</editor-fold>
                    //低松弛
                    coeW.b[X][Y][Z]
                            = coeW.b[X][Y][Z] + (URFFI - 1.0) * coeW.ap[X][Y][Z]
                            * F.getOldField()[X][Y][Z];
                    coeW.ap[X][Y][Z] = coeW.ap[X][Y][Z] * URFFI;
                }
            }
        }
        F.getLog().rl = MatrixSolver.gaussSeidel(F.getNewField(), coeW, sys);
        F.getLog().log("W");
    }

    /**
     * 计算湍流动能
     *
     * @param mesh 网格
     * @param F 场
     * @param Omega
     * @param W
     * @param mu
     * @param rho
     * @param Mut
     * @param S
     * @param sys
     */
    void cacK(Mesh mesh, Field F,
            double[][][] Omega, double[][][] W,
            double[][][] mu, double[][][] rho,
            double[][][] Mut, double[][][] S,
            SystemControl sys) {
        double URF = SOR[K];
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;
//        
        double sigma2 = 0.5, betastar = 0.09;
        double GEN = 0;
//        
        int i = mesh.numPx();
        int j = mesh.numPy();
        Coefficient coeW = new Coefficient(i, j);
        Label flag = new Label("K");
        for (int Z = 1; Z <= mesh.getNZ(); ++Z) {
            for (int Y = 1; Y <= mesh.getNY(); ++Y) {
                for (int X = 1; X <= mesh.getNX(); ++X) {
                    flag.setFlag(F.getBoundary(), X, Y);
                    double vol = mesh.getVol(X, Y, Z);
                    Fw = 0;
                    Fe = 0;
                    Fs = 0;
                    Fn = 0;
//                    
                    mup = (mu[X][Y][Z] + Mut[X][Y][Z] / sigma2);
                    muw = (mu[X - 1][Y][Z] + Mut[X - 1][Y][Z] / sigma2);
                    mue = (mu[X + 1][Y][Z] + Mut[X + 1][Y][Z] / sigma2);
                    mus = (mu[X][Y - 1][Z] + Mut[X][Y - 1][Z] / sigma2);
                    mun = (mu[X][Y + 1][Z] + Mut[X][Y + 1][Z] / sigma2);
//
                    gamw = faceValue(muw, mup - muw, mesh.gettPx()[X - 1], mesh.getTx()[X - 1], mesh.gettPx()[X]);
                    game = faceValue(mup, mue - mup, mesh.gettPx()[X], mesh.getTx()[X], mesh.gettPx()[X + 1]);
                    gams = faceValue(mus, mup - mus, mesh.gettPy()[Y - 1], mesh.getTy()[Y - 1], mesh.gettPy()[Y]);
                    gamn = faceValue(mup, mun - mup, mesh.gettPy()[Y], mesh.getTy()[Y], mesh.gettPy()[Y + 1]);
//                
                    Dw = gamw * mesh.getDYV()[Y] / mesh.getDXP()[X - 1];
                    De = game * mesh.getDYV()[Y] / mesh.getDXP()[X];
                    Ds = gams * mesh.getDXU()[X] / mesh.getDYP()[Y - 1];
                    Dn = gamn * mesh.getDXU()[X] / mesh.getDYP()[Y];
//                    
                    coeW.aw[X][Y][Z] = (1 - flag.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeW.ae[X][Y][Z] = (1 - flag.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeW.as[X][Y][Z] = (1 - flag.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeW.an[X][Y][Z] = (1 - flag.n) * (Math.max(-Fn, 0.0) + Dn);//n
//
                    Sp = -betastar * rho[X][Y][Z] * Omega[X][Y][Z];
                    Spad
                            = Spadw(F.getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, mesh.getDXP()[X - 1], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Spade(F.getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, mesh.getDXP()[X], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Spadw(F.getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, mesh.getDYP()[Y - 1], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Spade(F.getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, mesh.getDYP()[Y], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
//得到ap系数
                    coeW.ap[X][Y][Z]
                            = coeW.aw[X][Y][Z] + coeW.ae[X][Y][Z]
                            + coeW.as[X][Y][Z] + coeW.an[X][Y][Z]
                            + (Fe - Fw) + (Fn - Fs) - (Sp + Spad) * vol;
//                    
                    //<editor-fold>
                    if ((flag.w == 1) || (flag.e == 1) || (flag.s == 1) || (flag.n == 1)) {
                        double VISS = mu[X][Y][Z];
                        double YPL, DN = 0;
                        double TAU = 0;
                        double CTRANS = 11.63;
                        double CK = Turbulence.getCK(F.getOldField()[X][Y][Z]);
                        double VISCW, VISW;
                        if (flag.w == 1) {
                            DN = Tool.distance(
                                    mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realX(mesh.gettPx()[X - 1], mesh.gettPy()[Y]),
                                    mesh.realY(mesh.gettPx()[X - 1], mesh.gettPy()[Y]));
                            YPL = Turbulence.getYPL(rho[X][Y][Z], mu[X][Y][Z], CK, DN);
                            VISCW = Turbulence.getVISCW(YPL, mu[X][Y][Z]);
                            VISW = max(mu[X][Y][Z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = getTAU(VISS, W[X][Y][Z], W[X - 1][Y][Z], DN);
                        }
                        if (flag.e == 1) {
                            DN = Tool.distance(
                                    mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realX(mesh.gettPx()[X + 1], mesh.gettPy()[Y]),
                                    mesh.realY(mesh.gettPx()[X + 1], mesh.gettPy()[Y]));
                            YPL = Turbulence.getYPL(rho[X][Y][Z], mu[X][Y][Z], CK, DN);
                            VISCW = Turbulence.getVISCW(YPL, mu[X][Y][Z]);
                            VISW = max(mu[X][Y][Z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = getTAU(VISS, W[X][Y][Z], W[X + 1][Y][Z], DN);
                        }
                        if (flag.s == 1) {
                            DN = Tool.distance(
                                    mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y - 1]),
                                    mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y - 1]));
                            YPL = Turbulence.getYPL(rho[X][Y][Z], mu[X][Y][Z], CK, DN);
                            VISCW = Turbulence.getVISCW(YPL, mu[X][Y][Z]);
                            VISW = max(mu[X][Y][Z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = getTAU(VISS, W[X][Y][Z], W[X][Y - 1][Z], DN);
                        }
                        if (flag.n == 1) {
                            DN = Tool.distance(
                                    mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y + 1]),
                                    mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y + 1]));
                            YPL = Turbulence.getYPL(rho[X][Y][Z], mu[X][Y][Z], CK, DN);
                            VISCW = Turbulence.getVISCW(YPL, mu[X][Y][Z]);
                            VISW = max(mu[X][Y][Z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = getTAU(VISS, W[X][Y][Z], W[X][Y + 1][Z], DN);
                        }
                        GEN = Turbulence.getGEN(TAU, F.getOldField()[X][Y][Z], DN);
                    } else {
                        GEN = (Mut[X][Y][Z] * S[X][Y][Z]);
                    }
                    //</editor-fold>
//
                    Sc = GEN;
                    Scad
                            = Scadw(F.getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, mesh.getDXP()[X - 1], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Scade(F.getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, mesh.getDXP()[X], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Scadw(F.getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, mesh.getDYP()[Y - 1], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Scade(F.getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, mesh.getDYP()[Y], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
                    coeW.b[X][Y][Z] = (Sc + Scad) * vol;
//
                    coeW.b[X][Y][Z]
                            = coeW.b[X][Y][Z] + (URFFI - 1.0) * coeW.ap[X][Y][Z]
                            * F.getOldField()[X][Y][Z];
                    coeW.ap[X][Y][Z] = coeW.ap[X][Y][Z] * URFFI;
                }
            }
        }
        F.getLog().rl = MatrixSolver.gaussSeidel(F.getNewField(), coeW, sys);
        F.getLog().log("K");
    }

    /**
     * 计算湍流动能的比能量耗散率
     *
     * @param mesh 网格
     * @param F 场
     * @param K
     * @param W
     * @param mu
     * @param rho
     * @param Mut
     * @param S
     * @param sys
     */
    void cacOmega(Mesh mesh, Field F,
            double[][][] K, double[][][] W,
            double[][][] mu, double[][][] rho,
            double[][][] Mut, double[][][] S,
            SystemControl sys) {
        double URF = SOR[Omega];
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;
//        
        double sigma1 = 0.5, alpha = 0.555, beta = 0.075;
        int i = mesh.numPx();
        int j = mesh.numPy();
        Coefficient coeW = new Coefficient(i, j);
        Label flag = new Label("Omega");
        for (int Z = 1; Z <= mesh.getNZ(); ++Z) {
            for (int Y = 1; Y <= mesh.getNY(); ++Y) {
                for (int X = 1; X <= mesh.getNX(); ++X) {
                    flag.setFlag(F.getBoundary(), X, Y);
                    double vol = mesh.getVol(X, Y, Z);
                    Fw = 0;
                    Fe = 0;
                    Fs = 0;
                    Fn = 0;
//                    
                    mup = (mu[X][Y][Z] + Mut[X][Y][Z] / sigma1);
                    muw = (mu[X - 1][Y][Z] + Mut[X - 1][Y][Z] / sigma1);
                    mue = (mu[X + 1][Y][Z] + Mut[X + 1][Y][Z] / sigma1);
                    mus = (mu[X][Y - 1][Z] + Mut[X][Y - 1][Z] / sigma1);
                    mun = (mu[X][Y + 1][Z] + Mut[X][Y + 1][Z] / sigma1);
//                    
                    gamw = faceValue(muw, mup - muw, mesh.gettPx()[X - 1], mesh.getTx()[X - 1], mesh.gettPx()[X]);
                    game = faceValue(mup, mue - mup, mesh.gettPx()[X], mesh.getTx()[X], mesh.gettPx()[X + 1]);
                    gams = faceValue(mus, mup - mus, mesh.gettPy()[Y - 1], mesh.getTy()[Y - 1], mesh.gettPy()[Y]);
                    gamn = faceValue(mup, mun - mup, mesh.gettPy()[Y], mesh.getTy()[Y], mesh.gettPy()[Y + 1]);
//                  
                    Dw = gamw * mesh.getDYV()[Y] / mesh.getDXP()[X - 1];
                    De = game * mesh.getDYV()[Y] / mesh.getDXP()[X];
                    Ds = gams * mesh.getDXU()[X] / mesh.getDYP()[Y - 1];
                    Dn = gamn * mesh.getDXU()[X] / mesh.getDYP()[Y];
//                    
                    coeW.aw[X][Y][Z] = (1 - flag.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeW.ae[X][Y][Z] = (1 - flag.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeW.as[X][Y][Z] = (1 - flag.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeW.an[X][Y][Z] = (1 - flag.n) * (Math.max(-Fn, 0.0) + Dn);//n
//               
                    Sp = -F.getOldField()[X][Y][Z] * rho[X][Y][Z] * beta;
                    Spad
                            = Spadw(F.getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, mesh.getDXP()[X - 1], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Spade(F.getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, mesh.getDXP()[X], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Spadw(F.getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, mesh.getDYP()[Y - 1], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Spade(F.getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, mesh.getDYP()[Y], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
//得到ap系数
                    coeW.ap[X][Y][Z]
                            = coeW.aw[X][Y][Z] + coeW.ae[X][Y][Z]
                            + coeW.as[X][Y][Z] + coeW.an[X][Y][Z]
                            + (Fe - Fw) + (Fn - Fs) - (Sp + Spad) * vol;
//                    
                    Sc
                            = alpha * (Mut[X][Y][Z] * S[X][Y][Z])
                            * F.getOldField()[X][Y][Z] / K[X][Y][Z];
                    Scad
                            = Scadw(F.getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, mesh.getDXP()[X - 1], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Scade(F.getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, mesh.getDXP()[X], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Scadw(F.getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, mesh.getDYP()[Y - 1], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Scade(F.getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, mesh.getDYP()[Y], vol,
                                    F.getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
                    coeW.b[X][Y][Z] = (Sc + Scad) * vol;
                }
            }
        }
//      
        double DN = 1E-30;
        for (int Z = 1; Z <= mesh.getNZ(); ++Z) {
            for (int Y = 1; Y <= mesh.getNY(); ++Y) {
                for (int X = 1; X <= mesh.getNX(); ++X) {
                    flag.setFlag(F.getBoundary(), X, Y);
                    if ((flag.w == 1) || (flag.e == 1) || (flag.s == 1) || (flag.n == 1)) {
                        if (flag.w == 1) {
                            DN = Tool.distance(
                                    mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realX(mesh.gettPx()[X - 1], mesh.gettPy()[Y]),
                                    mesh.realY(mesh.gettPx()[X - 1], mesh.gettPy()[Y]));
                        }
                        if (flag.e == 1) {
                            DN = Tool.distance(
                                    mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realX(mesh.gettPx()[X + 1], mesh.gettPy()[Y]),
                                    mesh.realY(mesh.gettPx()[X + 1], mesh.gettPy()[Y]));
                        }
                        if (flag.s == 1) {
                            DN = Tool.distance(
                                    mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y - 1]),
                                    mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y - 1]));
                        }
                        if (flag.n == 1) {
                            DN = Tool.distance(
                                    mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y]),
                                    mesh.realX(mesh.gettPx()[X], mesh.gettPy()[Y + 1]),
                                    mesh.realY(mesh.gettPx()[X], mesh.gettPy()[Y + 1]));
                        }
                        coeW.ap[X][Y][Z] = 1.0;
                        coeW.aw[X][Y][Z] = 0;//w
                        coeW.ae[X][Y][Z] = 0;//e
                        coeW.as[X][Y][Z] = 0;//s
                        coeW.an[X][Y][Z] = 0;//n
                        coeW.b[X][Y][Z]
                                = 6.0 * mu[X][Y][Z] / (rho[X][Y][Z] * beta * DN * DN);
                    }
                    //under-relaxization  
                    coeW.b[X][Y][Z]
                            = coeW.b[X][Y][Z] + (URFFI - 1.0) * coeW.ap[X][Y][Z]
                            * F.getOldField()[X][Y][Z];
                    coeW.ap[X][Y][Z] = coeW.ap[X][Y][Z] * URFFI;
                }
            }
        }
//求解
        F.getLog().rl = MatrixSolver.gaussSeidel(F.getNewField(), coeW, sys);
        F.getLog().log("Omega");
    }

    /**
     * 输出文件控制
     *
     * @param count 迭代次数
     * @throws FileNotFoundException
     */
    void RESULT(int count) throws FileNotFoundException {
        String dirName = position + this.toString() + "/";
        Writer.createDir(dirName);
        String fileName = dirName + "RESULT" + count + "." + "DAT";
        Writer.createFile(fileName);
        fileWriter = new PrintWriter(fileName);

        fileWriter.println("Title=" + "\"" + "Field" + "\"");
        String var = "Variables=\"X\",\"Y\",\"W\",\"K\",\"Omega\",\"T\",\"S\",\"phi\",\"mut\"";
        fileWriter.println(var);
//        
        fileWriter.println("Zone"
                + " I=" + (mesh.getNX() + 2)
                + " J=" + (mesh.getNY() + 2)
                + " F=POINT");
        int K = 1;
        for (int J = 0; J <= mesh.getNY() + 1; ++J) {
            for (int I = 0; I <= mesh.getNX() + 1; ++I) {
                fileWriter.printf("%16.6E\t", mesh.realX(mesh.gettPx()[I], mesh.gettPy()[J]));
                fileWriter.printf("%16.6E\t", mesh.realY(mesh.gettPx()[I], mesh.gettPy()[J]));
                fileWriter.printf("%16.6E\t", F[W].getNewField()[I][J][K]);
                fileWriter.printf("%16.6E\t", F[K].getNewField()[I][J][K]);
                fileWriter.printf("%16.6E\t", F[Omega].getNewField()[I][J][K]);
                fileWriter.printf("%16.6E\t", F[T].getNewField()[I][J][K]);
                fileWriter.printf("%16.6E\t", pro[S].getF()[I][J][K]);
                fileWriter.printf("%16.6E\t", phi[I][J][K]);
                fileWriter.printf("%16.6E\t", pro[mut].getF()[I][J][K]);
                fileWriter.println();
            }
        }
        fileWriter.close();
    }

    String position, Title;

    @Override
    public String toString() {
        return getClass().getName() + Title;
    }
}
