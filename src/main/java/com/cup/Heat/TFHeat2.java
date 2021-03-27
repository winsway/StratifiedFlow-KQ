/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.Heat;

import com.cup.system.SystemControl;
import com.cup.field.Scalarfield;
import com.cup.field.Field;
import com.cup.util.MatrixSolver;
import static com.cup.boundary.BoundaryCondition.Scade;
import static com.cup.boundary.BoundaryCondition.Scadw;
import static com.cup.boundary.BoundaryCondition.Spade;
import static com.cup.boundary.BoundaryCondition.Spadw;
import com.cup.boundary.factory.Label;
import com.cup.boundary.factory.RobinBC;
import com.cup.system.ControlDict;
import com.cup.system.FvSolution;
import com.cup.field.Coefficient;
import com.cup.field.ProF;
import com.winswe.io.Writer;
import static com.cup.log.Residual.getRes;
import com.cup.mesh.Mesh;
import com.cup.mesh.factory.BiSameV2;
import com.cup.mesh.factory.BiSingleV2;
import static com.cup.util.Flux.faceValue;
import com.cup.util.Tool;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.max;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @see de Sampaio, P.A.B., J.L.H. Faccini and J. Su, Modelling of stratified
 * gas–liquid two-phase flow in horizontal circular pipes. International Journal
 * of Heat and Mass Transfer, 2008. 51(11-12): p. 2752-2761. 液-液
 * <n>1.双流体温度场计算</n>
 * <n>2.引入对流项的影响</n>
 * <n>3.未引入时间的影响</n>
 * @author 齐雪宇
 */
public class TFHeat2 implements Cloneable {

    FvSolution fvsolution = new FvSolution();
    ControlDict controldict = new ControlDict();
    SystemControl sysOil, sysWater;
    double dpdz, hl;
    double Qso, Qsw;
    double R, g = 0;
    ProF[] proW = new ProF[6], proO = new ProF[6];
    int S = 0, rho = 1, mu = 2, mut = 3, lambda = 4, cp = 5;
    Field[] FW = new Field[4], FO = new Field[4];
    Coefficient[] coeW = new Coefficient[4], coeO = new Coefficient[4];
    int W = 0, K = 1, Omega = 2, T = 3;
    Mesh mesh;
    BiSameV2 bm;
    BiSingleV2 OilMesh, WaterMesh;
    double[] SOR = {0.95, 0.9, 0.9, 0.9};
    TFHeat2 Xold, X_1old;
    double DRT, dz;
    /**
     * 需要引入一个逻辑变量控制时间和对流
     */
    Boolean time = Boolean.FALSE, convec = Boolean.FALSE;
    double[][] GAMW, GAMO;
    double[][] TimeW, TimeO;

    /**
     * 需要引入边界条件的基础参数
     */
    RobinBC bc;

    /**
     * 输出信息
     */
    String position, Title;
    /**
     * 物性参数
     */
    double[] testO;
    double[] testW;

    void setTimeConvection(TFHeat2 xold, TFHeat2 x_1old, double dt, double dz) {
        this.Xold = xold;
        this.X_1old = x_1old;
        this.dz = dz;
        if (dt != 0) {
            time = Boolean.TRUE;
            DRT = 1 / dt;
        } else {
            DRT = 0;
        }
        if (dz != 0) {
            convec = Boolean.TRUE;
        }

    }

    void setFluid(double[] testO, double[] testW) {
        this.testO = testO;
        this.testW = testW;
    }

    public void application(double dpdz, double hl, int I) throws FileNotFoundException {
        this.dpdz = dpdz;
        this.hl = hl;
        createGrid();
        createFields();
        initField();
        setBoundary();
        sysOil = new SystemControl(OilMesh, controldict, fvsolution);
        sysWater = new SystemControl(WaterMesh, controldict, fvsolution);
        Solver(I);
        System.out.println("求解完成！");
    }

    void start(double dpdz, double hl) {
        this.dpdz = dpdz;
        this.hl = hl;
        createGrid();
        createFields();
        initField();
        setBoundary();
        sysOil = new SystemControl(OilMesh, controldict, fvsolution);
        sysWater = new SystemControl(WaterMesh, controldict, fvsolution);
        System.out.println("初始化！");
    }

    void createFields() {
        for (int i = 0; i < FW.length; i++) {
            FW[i] = new Scalarfield(WaterMesh.numPx(), WaterMesh.numPy(), WaterMesh.numPz());
            FO[i] = new Scalarfield(OilMesh.numPx(), OilMesh.numPy(), OilMesh.numPz());
        }
        for (int i = 0; i < coeW.length; i++) {
            coeW[i] = new Coefficient(WaterMesh.numPx(), WaterMesh.numPy(), WaterMesh.numPz());
            coeO[i] = new Coefficient(OilMesh.numPx(), OilMesh.numPy(), OilMesh.numPz());
        }
        for (int i = 0; i < 6; i++) {
            proW[i] = new ProF(WaterMesh.numPx(), WaterMesh.numPy(), WaterMesh.numPz());
            proO[i] = new ProF(OilMesh.numPx(), OilMesh.numPy(), OilMesh.numPz());
        }
        System.out.println("场创建完成!");
    }

    void initField() {
//        油相设定初场
//<editor-fold>
        for (int K = 0; K <= OilMesh.getNZ() + 1; K++) {
            for (int I = 0; I <= OilMesh.getNX() + 1; I++) {
                for (int J = 0; J <= OilMesh.getNY() + 1; J++) {
                    FO[T].getNewField()[I][J][K] = TOil;
                    proO[rho].getF()[I][J][K] = testO[0];
                    proO[mu].getF()[I][J][K] = testO[1];
                    proO[lambda].getF()[I][J][K] = testO[2];
                    proO[cp].getF()[I][J][K] = testO[3];
                }
            }
        }
//</editor-fold>
//        水相设定初场
//<editor-fold>
        for (int K = 0; K <= WaterMesh.getNZ() + 1; K++) {
            for (int I = 0; I <= WaterMesh.getNX() + 1; I++) {
                for (int J = 0; J <= WaterMesh.getNY() + 1; J++) {
                    FW[T].getNewField()[I][J][K] = TWater;
                    proW[rho].getF()[I][J][K] = testW[0];
                    proW[mu].getF()[I][J][K] = testW[1];
                    proW[lambda].getF()[I][J][K] = testW[2];
                    proW[cp].getF()[I][J][K] = testW[3];
                }
            }
        }
//</editor-fold>
    }

    /**
     * 设定边界条件
     */
    void setBoundary() {
        RobinBC[] temp = new RobinBC[5];
        int xw = 0, xe = WaterMesh.getNX() + 1;
        int yw = 0, ye = WaterMesh.getNY() + 1;
        //X边界
        //<editor-fold>
        for (int Z = 0; Z <= WaterMesh.getNZ() + 1; Z++) {
            for (int Y = 0; Y <= WaterMesh.getNY() + 1; Y++) {
                //xw
                temp[W] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[W].getBoundary().getBoundaryRegion().setType(xw, Y, Z, 1, temp[W]);
                temp[K] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[K].getBoundary().getBoundaryRegion().setType(xw, Y, Z, 1, temp[K]);
                temp[Omega] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[Omega].getBoundary().getBoundaryRegion().setType(xw, Y, Z, 1, temp[Omega]);
                temp[T] = new RobinBC(59.0, 15.0, 375.0, 1.0);
                FW[T].getBoundary().getBoundaryRegion().setType(xw, Y, Z, 3, bc);
                //xe
                temp[W] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[W].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 1, temp[W]);
                temp[K] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[K].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 1, temp[K]);
                temp[Omega] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[Omega].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 1, temp[Omega]);
                temp[T] = new RobinBC(59.0, 15.0, 375.0, 1.0);
                FW[T].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 3, bc);
            }
        }
        //</editor-fold>
        //Y边界
        //<editor-fold>
        for (int Z = 0; Z <= WaterMesh.getNZ() + 1; Z++) {
            for (int X = 0; X <= WaterMesh.getNX() + 1; X++) {
                //yw
                temp[W] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[W].getBoundary().getBoundaryRegion().setType(X, yw, Z, 1, temp[W]);
                temp[K] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[K].getBoundary().getBoundaryRegion().setType(X, yw, Z, 1, temp[K]);
                temp[Omega] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[Omega].getBoundary().getBoundaryRegion().setType(X, yw, Z, 1, temp[Omega]);
                temp[T] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[T].getBoundary().getBoundaryRegion().setType(X, yw, Z, 1, bc);
                //ye
                temp[W] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[W].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, temp[W]);
                temp[K] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[K].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, temp[K]);
                temp[Omega] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[Omega].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, temp[Omega]);
                temp[T] = new RobinBC(59.0, 15.0, 375.0, 1.0);
                FW[T].getBoundary().getBoundaryRegion().setType(X, ye, Z, 3, bc);
            }
        }
        //</editor-fold>
        xw = 0;
        xe = OilMesh.getNX() + 1;
        yw = 0;
        ye = OilMesh.getNY() + 1;
        //X边界
        //<editor-fold>
        for (int Z = 0; Z <= OilMesh.getNZ() + 1; Z++) {
            for (int Y = 0; Y <= OilMesh.getNY() + 1; Y++) {
                //xw
                temp[W] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[W].getBoundary().getBoundaryRegion().setType(xw, Y, Z, 1, temp[W]);
                temp[K] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[K].getBoundary().getBoundaryRegion().setType(xw, Y, Z, 1, temp[K]);
                temp[Omega] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[Omega].getBoundary().getBoundaryRegion().setType(xw, Y, Z, 1, temp[Omega]);
                temp[T] = new RobinBC(59.0, 15.0, 375.0, 1.0);
                FO[T].getBoundary().getBoundaryRegion().setType(xw, Y, Z, 3, bc);
                //xe
                temp[W] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[W].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 1, temp[W]);
                temp[K] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[K].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 1, temp[K]);
                temp[Omega] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[Omega].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 1, temp[Omega]);
                temp[T] = new RobinBC(59.0, 15.0, 375.0, 1.0);
                FO[T].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 3, bc);
            }
        }
        //</editor-fold>
        //Y边界
        //<editor-fold>
        for (int Z = 0; Z <= OilMesh.getNZ() + 1; Z++) {
            for (int X = 0; X <= OilMesh.getNX() + 1; X++) {
                //yw
                temp[W] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[W].getBoundary().getBoundaryRegion().setType(X, yw, Z, 1, temp[W]);
                temp[K] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[K].getBoundary().getBoundaryRegion().setType(X, yw, Z, 1, temp[K]);
                temp[Omega] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[Omega].getBoundary().getBoundaryRegion().setType(X, yw, Z, 1, temp[Omega]);
                temp[T] = new RobinBC(59.0, 15.0, 375.0, 1.0);
                FO[T].getBoundary().getBoundaryRegion().setType(X, yw, Z, 3, bc);
                //ye
                temp[W] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[W].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, temp[W]);
                temp[K] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[K].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, temp[K]);
                temp[Omega] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[Omega].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, temp[Omega]);
                temp[T] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[T].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, bc);
            }
        }
        //</editor-fold>
        System.out.println("边界设定完成!");
    }

    public boolean createGrid() {
        bm = new BiSameV2(12, Math.PI, 100, 100);
        bm.radia = R;
        double D = bm.radia * 2.0;
        hl = D * hl;
        double gamma = acos(1.0 - 2 * hl / D);
        bm.theta1 = gamma;
        bm.theta2 = PI;
        mesh = bm;
        mesh.blockMesh();
//       设定两相网格
//       油相
        System.out.println("begin OilMesh");
        OilMesh = new BiSingleV2(12, bm.theta2 - bm.theta1,
                bm.getNX(), bm.interhL);
        OilMesh.theta1 = bm.theta1;
        OilMesh.theta2 = bm.theta2;
        OilMesh.a = bm.a;
        OilMesh.radia = bm.radia;
        OilMesh.blockMesh();
//       水相
        System.out.println("begin WaterMesh");
        WaterMesh = new BiSingleV2(12, bm.theta3 - bm.theta2,
                bm.getNX(), bm.getNY() - bm.interhL);
        WaterMesh.theta1 = bm.theta2;
        WaterMesh.theta2 = bm.theta3;
        WaterMesh.a = bm.a;
        WaterMesh.radia = bm.radia;
        WaterMesh.blockMesh();
        return true;
    }

    public boolean Solver(int i) throws FileNotFoundException {
        time();
        convection();
        if (i == 1) {
            scheme1();
        }
        if (i == 2) {
            scheme2();
        }
        if (i == 3) {
            scheme3();
        }
        return true;
    }

    double[] RES = new double[4];

    /**
     * 单相双流体层流方程求解
     */
    void scheme1() throws FileNotFoundException {
        int count = 0;
        double Qw, Qo;
        RESULT(count);
//        进入大循环
        do {
            count++;
            updateProF();
            RES[W] = cacVelocity();
        } while (getMaxRES() > 1e-3 || count < 500);
        Qw = Tool.getTotalVol(FW[W], WaterMesh);
        Qo = Tool.getTotalVol(FO[W], OilMesh);
        System.out.println("计算总流量=Qw" + Qw);
        System.out.println("计算总流量=Qo" + Qo);
        System.out.println("总流量=" + (Qw + Qo));
        RESULT(count);
    }

    /**
     * 双流体方程+温度场
     */
    void scheme2() throws FileNotFoundException {
        int count = 0;
        RESULT(count);
        do {
            count++;
            updateProF();
            RES[T] = cacTemperature();
        } while (RES[T] > 1e-3 && count < 10000);
        RESULT(count);
    }

    /**
     * 速度场方程+温度场
     *
     * @throws FileNotFoundException
     */
    void scheme3() throws FileNotFoundException {
        int count = 0;
        RESULT(count);
//        进入大循环
        do {
            count++;
            updateProF();
            RES[W] = cacVelocity();
            RES[T] = cacTemperature();
        } while (getMaxRES() > 1e-3 || count < 500);
        RESULT(count);
    }

    void updateProF() {
        for (int Z = 1; Z <= WaterMesh.getNZ(); ++Z) {
            for (int Y = 1; Y <= WaterMesh.getNY(); ++Y) {
                for (int X = 1; X <= WaterMesh.getNX(); ++X) {
                    //更新粘度
                    updateMu(X, Y, Z);
                    //更新定压比热容;
                    updateCp(X, Y, Z);
                }
            }
        }
    }

    double cacVelocity() {
        double resw = 0, reso = 0;
        System.out.println("计算速度场W");
        cacW(dpdz, WaterMesh, FW[W], proW[mu].getF(),
                proW[rho].getF(), proW[mut].getF(), sysWater);
        cacW(dpdz, OilMesh, FO[W], proO[mu].getF(),
                proO[rho].getF(), proO[mut].getF(), sysOil);
        updateInterfaceW();
        resw = getRes(FW[W].getNewField(), FW[W].getOldField());
        reso = getRes(FO[W].getNewField(), FO[W].getOldField());
        FW[W].copyNew2Old();
        FO[W].copyNew2Old();
        System.out.println("resw=" + resw + " reso=" + reso);
        return max(resw, reso);
    }

    double cacTemperature() {
        double resw, reso;
        System.out.println("计算温度场");
//      <editor-fold>
        cacT(OilMesh, FO[T], proO[lambda].getF(), proO[cp].getF(),
                proO[rho].getF(), proO[mut].getF(), GAMO, TimeO,
                sysOil);
        cacT(WaterMesh, FW[T], proW[lambda].getF(), proW[cp].getF(),
                proW[rho].getF(), proW[mut].getF(), GAMW, TimeW,
                sysWater);
        updateInterfaceT();
        resw = getRes(FW[T].getNewField(), FW[T].getOldField());
        reso = getRes(FO[T].getNewField(), FO[T].getOldField());
        FW[T].copyNew2Old();
        FO[T].copyNew2Old();
//     </editor-fold>
        System.out.println("resw=" + resw + " reso=" + reso);
        return max(resw, reso);
    }

    double getMaxRES() {
        double max = RES[0];
        System.out.println("RESW = " + RES[W]);
        System.out.println("RESK = " + RES[K]);
        System.out.println("RESOmega = " + RES[Omega]);
        System.out.println("REST = " + RES[T]);
        for (int i = 0; i < RES.length; i++) {
            max = RES[i] > max ? RES[i] : max;
        }
        return max;
    }

    /**
     * 计算速度场
     *
     * @param dpdz
     */
    void cacW(double dpdz,
            Mesh mesh, Field W,
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
        Label flagW = new Label("W");
        int i = mesh.numPx();
        int j = mesh.numPy();
        Coefficient coeW = new Coefficient(i, j);
        for (int Z = 1; Z <= WaterMesh.getNZ(); ++Z) {
            for (int Y = 1; Y <= WaterMesh.getNY(); ++Y) {
                for (int X = 1; X <= WaterMesh.getNX(); ++X) {
                    //<editor-fold>
                    flagW.setFlag(W.getBoundary(), X, Y);
                    double vol
                            = mesh.getDXU()[X] * mesh.getDYV()[Y]
                            * mesh.J(mesh.gettPx()[X], mesh.gettPy()[Y]);
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
                    coeW.aw[X][Y][Z] = (1 - flagW.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeW.ae[X][Y][Z] = (1 - flagW.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeW.as[X][Y][Z] = (1 - flagW.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeW.an[X][Y][Z] = (1 - flagW.n) * (Math.max(-Fn, 0.0) + Dn);//n
//设定附加源项
                    Sp = 0;
                    Spad
                            = Spadw(W.getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, mesh.getDXP()[X - 1], vol,
                                    W.getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Spade(W.getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, mesh.getDXP()[X], vol,
                                    W.getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Spadw(W.getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, mesh.getDYP()[Y - 1], vol,
                                    W.getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Spade(W.getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, mesh.getDYP()[Y], vol,
                                    W.getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
//得到ap系数
                    coeW.ap[X][Y][Z]
                            = coeW.aw[X][Y][Z] + coeW.ae[X][Y][Z]
                            + coeW.as[X][Y][Z] + coeW.an[X][Y][Z]
                            + (Fe - Fw) + (Fn - Fs) - (Sp + Spad) * vol;
//得到源项
                    Sc = dpdz - rho[X][Y][Z] * g * Math.sin(0);
                    Scad
                            = Scadw(W.getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, mesh.getDXP()[X - 1], vol,
                                    W.getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Scade(W.getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, mesh.getDXP()[X], vol,
                                    W.getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Scadw(W.getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, mesh.getDYP()[Y - 1], vol,
                                    W.getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Scade(W.getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, mesh.getDYP()[Y], vol,
                                    W.getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
                    coeW.b[X][Y][Z] = (Sc + Scad) * vol;
                    //</editor-fold>
//under-relaxization
                    coeW.b[X][Y][Z]
                            = coeW.b[X][Y][Z] + (URFFI - 1.0) * coeW.ap[X][Y][Z]
                            * W.getOldField()[X][Y][Z];
                    coeW.ap[X][Y][Z] = coeW.ap[X][Y][Z] * URFFI;
                }
            }
        }
        W.getLog().rl = MatrixSolver.gaussSeidel(W.getNewField(), coeW, sys);
        W.getLog().log("W");
    }

    /**
     * 计算温度场
     *
     * @param dpdz
     */
    void cacT(Mesh mesh, Field T,
            double[][][] lambda, double[][][] cp,
            double[][][] rho, double[][][] Mut,
            double[][] GAM, double[][] Time,
            SystemControl sys
    ) {
        double URF = SOR[this.T];
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;
        double simgaT = 0.7;
        Label flagT = new Label("T");
        int i = mesh.numPx();
        int j = mesh.numPy();
        Coefficient coeT = new Coefficient(i, j);
        for (int Z = 1; Z <= WaterMesh.getNZ(); ++Z) {
            for (int Y = 1; Y <= WaterMesh.getNY(); ++Y) {
                for (int X = 1; X <= WaterMesh.getNX(); ++X) {
                    //<editor-fold>
                    flagT.setFlag(T.getBoundary(), X, Y);
                    double vol = mesh.getVol(X, Y, Z);
                    Fw = 0;
                    Fe = 0;
                    Fs = 0;
                    Fn = 0;
//
                    mup = (lambda[X][Y][Z] + Mut[X][Y][Z] * cp[X][Y][Z] / simgaT);
                    muw = (lambda[X - 1][Y][Z] + Mut[X - 1][Y][Z] * cp[X - 1][Y][Z] / simgaT);
                    mue = (lambda[X + 1][Y][Z] + Mut[X + 1][Y][Z] * cp[X + 1][Y][Z] / simgaT);
                    mus = (lambda[X][Y - 1][Z] + Mut[X][Y - 1][Z] * cp[X][Y - 1][Z] / simgaT);
                    mun = (lambda[X][Y + 1][Z] + Mut[X][Y + 1][Z] * cp[X][Y + 1][Z] / simgaT);
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
                    coeT.aw[X][Y][Z] = (1 - flagT.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeT.ae[X][Y][Z] = (1 - flagT.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeT.as[X][Y][Z] = (1 - flagT.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeT.an[X][Y][Z] = (1 - flagT.n) * (Math.max(-Fn, 0.0) + Dn);//n
//设定附加源项
                    Sp = 0;
                    Spad
                            = Spadw(T.getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, mesh.getDXP()[X - 1], vol,
                                    T.getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Spade(T.getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, mesh.getDXP()[X], vol,
                                    T.getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Spadw(T.getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, mesh.getDYP()[Y - 1], vol,
                                    T.getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Spade(T.getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, mesh.getDYP()[Y], vol,
                                    T.getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
//得到ap系数
                    coeT.ap[X][Y][Z]
                            = DRT * vol * rho[X][Y][Z] * cp[X][Y][Z]
                            + coeT.aw[X][Y][Z] + coeT.ae[X][Y][Z]
                            + coeT.as[X][Y][Z] + coeT.an[X][Y][Z]
                            + (Fe - Fw) + (Fn - Fs) - (Sp + Spad) * vol;
//得到源项
                    Sc = GAM[X][Y];
                    Scad
                            = Scadw(T.getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, mesh.getDXP()[X - 1], vol,
                                    T.getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Scade(T.getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, mesh.getDXP()[X], vol,
                                    T.getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Scadw(T.getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, mesh.getDYP()[Y - 1], vol,
                                    T.getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Scade(T.getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, mesh.getDYP()[Y], vol,
                                    T.getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
                    coeT.b[X][Y][Z]
                            = Time[X][Y] * vol + (Sc + Scad) * vol;
                    //</editor-fold>
//under-relaxization
                    coeT.b[X][Y][Z]
                            = coeT.b[X][Y][Z] + (URFFI - 1.0) * coeT.ap[X][Y][Z]
                            * T.getOldField()[X][Y][Z];
                    coeT.ap[X][Y][Z] = coeT.ap[X][Y][Z] * URFFI;
                }
            }
        }
        T.getLog().rl = MatrixSolver.gaussSeidel(T.getNewField(), coeT, sys);
        T.getLog().log("T");
    }

    void time() {
        TimeW = new double[WaterMesh.numPx()][WaterMesh.numPy()];
        TimeO = new double[OilMesh.numPx()][OilMesh.numPy()];
        if (time) {
            int Z = 1;
            //<editor-fold>
            for (int Y = 1; Y <= WaterMesh.getNY(); ++Y) {
                for (int X = 1; X <= WaterMesh.getNX(); ++X) {
                    TimeW[X][Y]
                            = Xold.proW[rho].getF()[X][Y][Z]
                            * Xold.proW[cp].getF()[X][Y][Z]
                            * Xold.FW[T].getNewField()[X][Y][Z]
                            * DRT;
                }
            }
            //</editor-fold>
            //<editor-fold>
            for (int Y = 1; Y <= OilMesh.getNY(); ++Y) {
                for (int X = 1; X <= OilMesh.getNX(); ++X) {
                    TimeO[X][Y]
                            = Xold.proO[rho].getF()[X][Y][Z]
                            * Xold.proO[cp].getF()[X][Y][Z]
                            * Xold.FO[T].getNewField()[X][Y][Z]
                            * DRT;
                }
            }
            //</editor-fold>
        }
    }

    void convection() {
        GAMW = new double[WaterMesh.numPx()][WaterMesh.numPy()];
        GAMO = new double[OilMesh.numPx()][OilMesh.numPy()];
        if (convec) {
            double phiC, phiU, Ac, Au;
            int Z = 1;
            //<editor-fold>
            for (int Y = 1; Y <= WaterMesh.getNY(); ++Y) {
                for (int X = 1; X <= WaterMesh.getNX(); ++X) {
                    phiC
                            = Xold.proW[rho].getF()[X][Y][Z]
                            * Xold.proW[cp].getF()[X][Y][Z]
                            * Xold.FW[W].getNewField()[X][Y][Z]
                            * Xold.FW[T].getNewField()[X][Y][Z];
                    phiU
                            = X_1old.proW[rho].getF()[X][Y][Z]
                            * X_1old.proW[cp].getF()[X][Y][Z]
                            * X_1old.FW[W].getNewField()[X][Y][Z]
                            * X_1old.FW[T].getNewField()[X][Y][Z];
                    Ac
                            = Xold.WaterMesh.getVol(X, Y, Z);
                    Au
                            = X_1old.WaterMesh.getVol(X, Y, Z);
                    GAMW[X][Y] = getCon(phiC, phiU, 1.0, 1.0, dz);
                }
            }
            //</editor-fold>
            //<editor-fold>
            for (int Y = 1; Y <= OilMesh.getNY(); ++Y) {
                for (int X = 1; X <= OilMesh.getNX(); ++X) {
                    phiC
                            = Xold.proO[rho].getF()[X][Y][Z]
                            * Xold.proO[cp].getF()[X][Y][Z]
                            * Xold.FO[W].getNewField()[X][Y][Z]
                            * Xold.FO[T].getNewField()[X][Y][Z];
                    phiU
                            = X_1old.proO[rho].getF()[X][Y][Z]
                            * X_1old.proO[cp].getF()[X][Y][Z]
                            * X_1old.FO[W].getNewField()[X][Y][Z]
                            * X_1old.FO[T].getNewField()[X][Y][Z];
                    Ac
                            = Xold.OilMesh.getVol(X, Y, Z);
                    Au
                            = X_1old.OilMesh.getVol(X, Y, Z);
                    GAMO[X][Y] = getCon(phiC, phiU, 1.0, 1.0, dz);
                }
            }
            //</editor-fold>
        }
    }

    double getCon(double phiC, double phiU, double Ac, double Au, double dz) {
        return -(phiC - Au / Ac * phiU) / dz;
    }

    /**
     * 引入温度因素改变粘度
     */
    void updateMu(int x, int y, int z) {

    }

    /**
     * 更新定压比热容
     */
    void updateCp(int x, int y, int z) {

    }

    /**
     * 更新界面速度和边界条件
     */
    void updateInterfaceW() {
        for (int x = 1; x <= mesh.getNX(); ++x) {
            double yOil = Tool.distance(
                    mesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY() + 1]),
                    mesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY() + 1]),
                    mesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY()]),
                    mesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY()]));
            double yWater = Tool.distance(
                    mesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[1]),
                    mesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[1]),
                    mesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[0]),
                    mesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[0]));
            double w = (proW[mu].getF()[x][1][1] + proW[mut].getF()[x][1][1]) / yWater;
            double o = (proO[mu].getF()[x][OilMesh.getNY()][1] + proO[mut].getF()[x][OilMesh.getNY()][1]) / yOil;
            double u = (FO[W].getNewField()[x][OilMesh.getNY()][1] * o + FW[W].getNewField()[x][1][1] * w)
                    / (o + w);
            RobinBC robinbc = new RobinBC(0.0, 1.0, u, 1.0);
            FO[W].setNewField(x, OilMesh.getNY() + 1, 1, u);
            FO[W].getBoundary().getBoundaryRegion().setType(x, OilMesh.getNY() + 1, 1, 1, robinbc);
            FW[W].setNewField(x, 0, 1, u);
            FW[W].getBoundary().getBoundaryRegion().setType(x, 0, 1, 1, robinbc);
        }
    }

    /**
     * 更新界面温度和边界条件
     */
    void updateInterfaceT() {
        for (int x = 1; x <= mesh.getNX(); ++x) {
            double yOil = Tool.distance(
                    mesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY() + 1]),
                    mesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY() + 1]),
                    mesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY()]),
                    mesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY()]));
            double yWater = Tool.distance(
                    mesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[1]),
                    mesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[1]),
                    mesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[0]),
                    mesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[0]));
            double w = (-proW[lambda].getF()[x][1][1]) / yWater;
            double o = (-proO[lambda].getF()[x][OilMesh.getNY()][1]) / yOil;
            double t = (FO[T].getNewField()[x][OilMesh.getNY()][1] * o + FW[T].getNewField()[x][1][1] * w)
                    / (o + w);
            RobinBC robinbc = new RobinBC(0.0, 1.0, t, 1.0);
            FO[T].setNewField(x, OilMesh.getNY() + 1, 1, t);
            FO[T].getBoundary().getBoundaryRegion().setType(x, OilMesh.getNY() + 1, 1, 1, robinbc);
            FW[T].setNewField(x, 0, 1, t);
            FW[T].getBoundary().getBoundaryRegion().setType(x, 0, 1, 1, robinbc);
        }
        int xw = 0, xe = WaterMesh.getNX() + 1;
        int yw = 0, ye = WaterMesh.getNY() + 1;
        double dxw = WaterMesh.getDXP()[0], dxe = WaterMesh.getDXP()[WaterMesh.getNX()];
        double dyw = WaterMesh.getDYP()[0], dye = WaterMesh.getDYP()[WaterMesh.getNY()];
        double value;
        //X边界
        //<editor-fold>
        for (int Z = 0; Z <= WaterMesh.getNZ() + 1; Z++) {
            for (int Y = 0; Y <= WaterMesh.getNY() + 1; Y++) {
                //xw
                value = FW[T].getBoundary().getBoundaryRegion().parameter[xw][Y][Z].newValue(dxw, FW[T].getNewField()[xw + 1][Y][Z]);
                FW[T].setNewField(xw, Y, Z, value);
                //xe
                value = FW[T].getBoundary().getBoundaryRegion().parameter[xe][Y][Z].newValue(dxe, FW[T].getNewField()[xe - 1][Y][Z]);
                FW[T].setNewField(xe, Y, Z, value);
            }
        }
        //</editor-fold>
        //Y边界
        //<editor-fold>
        for (int Z = 0; Z <= WaterMesh.getNZ() + 1; Z++) {
            for (int X = 0; X <= WaterMesh.getNX() + 1; X++) {
                //xw
                value = FW[T].getBoundary().getBoundaryRegion().parameter[X][yw][Z].newValue(dyw, FW[T].getNewField()[X][yw + 1][Z]);
                FW[T].setNewField(X, yw, Z, value);
                //xe
                value = FW[T].getBoundary().getBoundaryRegion().parameter[X][ye][Z].newValue(dye, FW[T].getNewField()[X][ye - 1][Z]);
                FW[T].setNewField(X, ye, Z, value);
            }
        }
        //</editor-fold>
        xw = 0;
        xe = OilMesh.getNX() + 1;
        yw = 0;
        ye = OilMesh.getNY() + 1;

        dxw = OilMesh.getDXP()[0];
        dxe = OilMesh.getDXP()[OilMesh.getNX()];
        dyw = OilMesh.getDYP()[0];
        dye = OilMesh.getDYP()[OilMesh.getNY()];
        //X边界
        //<editor-fold>
        for (int Z = 0; Z <= OilMesh.getNZ() + 1; Z++) {
            for (int Y = 0; Y <= OilMesh.getNY() + 1; Y++) {
                //xw
                value = FO[T].getBoundary().getBoundaryRegion().parameter[xw][Y][Z].newValue(dxw, FO[T].getNewField()[xw + 1][Y][Z]);
                FO[T].setNewField(xw, Y, Z, value);
                //xe
                value = FO[T].getBoundary().getBoundaryRegion().parameter[xe][Y][Z].newValue(dxe, FO[T].getNewField()[xe - 1][Y][Z]);
                FO[T].setNewField(xe, Y, Z, value);
            }
        }
        //</editor-fold>
        //Y边界
        //<editor-fold>
        for (int Z = 0; Z <= OilMesh.getNZ() + 1; Z++) {
            for (int X = 0; X <= OilMesh.getNX() + 1; X++) {
                //xw
                value = FO[T].getBoundary().getBoundaryRegion().parameter[X][yw][Z].newValue(dyw, FO[T].getNewField()[X][yw + 1][Z]);
                FO[T].setNewField(X, yw, Z, value);
                //xe
                value = FO[T].getBoundary().getBoundaryRegion().parameter[X][ye][Z].newValue(dye, FO[T].getNewField()[X][ye - 1][Z]);
                FO[T].setNewField(X, ye, Z, value);
            }
        }
        //</editor-fold>
    }

    public TFHeat2() {
    }

    public TFHeat2(double Qso, double Qsw, double Radia) {
        this.Qso = Qso;
        this.Qsw = Qsw;
        this.R = Radia;
    }

    double TOil, TWater;

    public TFHeat2(double Qso, double Qsw, double TOil, double TWater, double Radia) {
        this.Qso = Qso;
        this.Qsw = Qsw;
        this.TOil = TOil;
        this.TWater = TWater;
        this.R = Radia;
    }

    public TFHeat2(double Qso, double Qsw, double TOil, double TWater, RobinBC BC, double Radia) {
        this.Qso = Qso;
        this.Qsw = Qsw;
        this.TOil = TOil;
        this.TWater = TWater;
        try {
            this.bc = BC.clone();

        } catch (CloneNotSupportedException ex) {
            Logger.getLogger(TFHeat2.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        this.R = Radia;
    }

    /**
     * @param Qw Qw水相体积流量
     * @return
     */
    public double F(double Qw) {
        double Qcwater = 0;
        Qcwater = Tool.getTotalVol(FW[W], WaterMesh);
        System.out.println("Qcwate = " + Qcwater);
        return Qcwater - Qw;
    }

    /**
     *
     * @param Qoil 油相体积流量
     * @return
     */
    public double G(double Qoil) {
        double Qcoil = 0;
        Qcoil = Tool.getTotalVol(FO[W], OilMesh);
        System.out.println("Qcoil = " + Qcoil);
        return Qcoil - Qoil;
    }

    /**
     * 输出文件控制
     *
     * @param count 迭代次数
     * @throws FileNotFoundException
     */
    void RESULT(int count) throws FileNotFoundException {
        String dirName = position + this.toString() + (1.0 / (DRT + 1E-30)) + "/";
        Writer.createDir(dirName);
        String fileName = dirName + "RESULT" + count + "." + "DAT";
        Writer.createFile(fileName);
        PrintWriter RES
                = new PrintWriter(new OutputStreamWriter(new FileOutputStream(fileName)));
        RES.println("Title=" + "\"" + "Field" + "\"");
        String var = "Variables=\"X\",\"Y\",\"W\",\"K\",\"Omega\",\"T\",\"S\"";
        RES.println(var);
//        
        RES.println("Zone"
                + " I=" + (WaterMesh.getNX() + 2)
                + " J=" + (WaterMesh.getNY() + 2)
                + " F=POINT");
        int k = 1;
        for (int J = 0; J <= WaterMesh.getNY() + 1; ++J) {
            for (int I = 0; I <= WaterMesh.getNX() + 1; ++I) {
                RES.printf("%16.6E\t", WaterMesh.realX(WaterMesh.gettPx()[I], WaterMesh.gettPy()[J]));
                RES.printf("%16.6E\t", WaterMesh.realY(WaterMesh.gettPx()[I], WaterMesh.gettPy()[J]));
                RES.printf("%16.6E\t", FW[W].getNewField()[I][J][k]);
                RES.printf("%16.6E\t", FW[K].getNewField()[I][J][k]);
                RES.printf("%16.6E\t", FW[Omega].getNewField()[I][J][k]);
                RES.printf("%16.6E\t", FW[T].getNewField()[I][J][k]);
                RES.printf("%16.6E\t", proW[S].getF()[I][J][k]);
                RES.println();
            }
        }
        RES.println("Zone"
                + " I=" + (OilMesh.getNX() + 2)
                + " J=" + (OilMesh.getNY() + 2)
                + " F=POINT");

        for (int J = 0; J <= OilMesh.getNY() + 1; ++J) {
            for (int I = 0; I <= OilMesh.getNX() + 1; ++I) {
                RES.printf("%16.6E\t", OilMesh.realX(OilMesh.gettPx()[I], OilMesh.gettPy()[J]));
                RES.printf("%16.6E\t", OilMesh.realY(OilMesh.gettPx()[I], OilMesh.gettPy()[J]));
                RES.printf("%16.6E\t", FO[W].getNewField()[I][J][k]);
                RES.printf("%16.6E\t", FO[K].getNewField()[I][J][k]);
                RES.printf("%16.6E\t", FO[Omega].getNewField()[I][J][k]);
                RES.printf("%16.6E\t", FO[T].getNewField()[I][J][k]);
                RES.printf("%16.6E\t", proO[S].getF()[I][J][k]);
                RES.println();
            }
        }
        RES.close();
    }

    @Override
    public TFHeat2 clone() throws CloneNotSupportedException {
        TFHeat2 clone = new TFHeat2();
        clone.fvsolution = this.fvsolution.clone();
        clone.controldict = this.controldict.clone();
        clone.sysOil = this.sysOil.clone();
        clone.sysWater = this.sysWater.clone();
        clone.dpdz = this.dpdz;
        clone.hl = this.hl;
        clone.Qso = this.Qso;
        clone.Qsw = this.Qsw;
        clone.TOil = this.TOil;
        clone.TWater = this.TWater;
        clone.R = this.R;
        clone.g = this.g;
        clone.setFluid(testO, testW);
        for (int i = 0; i < proW.length; i++) {
            clone.proW[i] = this.proW[i].clone();
        }
        for (int i = 0; i < proO.length; i++) {
            clone.proO[i] = this.proO[i].clone();
        }
        for (int i = 0; i < FW.length; i++) {
            clone.FW[i] = this.FW[i].clone();
        }
        for (int i = 0; i < FO.length; i++) {
            clone.FO[i] = this.FO[i].clone();
        }
        for (int i = 0; i < coeW.length; i++) {
            clone.coeW[i] = this.coeW[i].clone();
        }
        for (int i = 0; i < coeO.length; i++) {
            clone.coeO[i] = this.coeO[i].clone();
        }
        clone.mesh = clone.bm;
        clone.bm = this.bm.clone();
        clone.OilMesh = this.OilMesh.clone();
        clone.WaterMesh = this.WaterMesh.clone();
        clone.bc = this.bc.clone();
        return clone;
    }

    @Override
    public String toString() {
        return getClass().getName() + Title;
    }
}
