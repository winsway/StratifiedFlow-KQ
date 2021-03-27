/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.KOmega;

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
import com.cup.turbulence.Turbulence;
import static com.cup.util.Flux.faceValue;
import com.cup.util.Tool;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.exp;
import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

/**
 *
 * @see de Sampaio, P.A.B., J.L.H. Faccini and J. Su, Modelling of stratified
 * gas–liquid two-phase flow in horizontal circular pipes. International Journal
 * of Heat and Mass Transfer, 2008. 51(11-12): p. 2752-2761. 液-液
 * <n>1.引入界面摩擦因素</n>
 * <n>2.解决压力对不上的问题</n>
 * <n>3.解决速度场对不上的问题</n>
 * <n>4.采用新的方式更改界面更新方法，避免计算速度降低，同时考虑界面摩擦阻力</n>
 * @author 齐雪宇
 */
public final class komegaTFOWAngeli2Wave1 {

    FvSolution fvsolution = new FvSolution();
    ControlDict controldict = new ControlDict();
    SystemControl sys, sysOil, sysWater;
    double dpdz, hl;
    double Qso, Qsw;
    double R, g = 0;
    ProF[] porW = new ProF[6], porO = new ProF[6];
    int S = 0, rho = 1, mu = 2, mut = 3, lambda = 4, cp = 5;
    Field[] FW = new Field[4], FO = new Field[4];
    Coefficient[] coeW = new Coefficient[4], coeO = new Coefficient[4];
    int W = 0, K = 1, Omega = 2, T = 3;
    Mesh mesh;
    BiSameV2 bm;
    BiSingleV2 OilMesh, WaterMesh;

    public String position, Title;
    /**
     * 物性参数
     */
    double[] testO;
    double[] testW;

    public void setFluid(double[] testO, double[] testW) {
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
            porW[i] = new ProF(WaterMesh.numPx(), WaterMesh.numPy(), WaterMesh.numPz());
            porO[i] = new ProF(OilMesh.numPx(), OilMesh.numPy(), OilMesh.numPz());
        }
        System.out.println("场创建完成!");
    }

    void initField() {
//        油相设定初场
//<editor-fold>
        for (int K = 0; K <= OilMesh.getNZ() + 1; K++) {
            for (int I = 0; I <= OilMesh.getNX() + 1; I++) {
                for (int J = 0; J <= OilMesh.getNY() + 1; J++) {
                    porO[rho].getF()[I][J][K] = testO[0];
                    porO[mu].getF()[I][J][K] = testO[1];
                    porO[lambda].getF()[I][J][K] = testO[2];
                    porO[cp].getF()[I][J][K] = testO[3];
                }
            }
        }
//</editor-fold>
//        水相设定初场
//<editor-fold>
        for (int K = 0; K <= WaterMesh.getNZ() + 1; K++) {
            for (int I = 0; I <= WaterMesh.getNX() + 1; I++) {
                for (int J = 0; J <= WaterMesh.getNY() + 1; J++) {
                    porW[rho].getF()[I][J][K] = testW[0];
                    porW[mu].getF()[I][J][K] = testW[1];
                    porW[lambda].getF()[I][J][K] = testW[2];
                    porW[cp].getF()[I][J][K] = testW[3];
                }
            }
        }
//</editor-fold>
    }

    void setBoundary() {
        RobinBC[] temp = new RobinBC[5];
        int xw = 0, xe = WaterMesh.getNX() + 1;
        int yw = 0, ye = WaterMesh.getNY() + 1;
        int zw = 0, ze = WaterMesh.getNZ() + 1;
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
                temp[T] = new RobinBC(1.0, 1.0, 0.0, 1.0);
                FW[T].getBoundary().getBoundaryRegion().setType(xw, Y, Z, 3, temp[T]);
                //xe
                temp[W] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[W].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 1, temp[W]);
                temp[K] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[K].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 1, temp[K]);
                temp[Omega] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[Omega].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 1, temp[Omega]);
                temp[T] = new RobinBC(1.0, 1.0, 0.0, 1.0);
                FW[T].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 3, temp[T]);
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
                FW[T].getBoundary().getBoundaryRegion().setType(X, yw, Z, 1, temp[T]);
                //ye
                temp[W] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[W].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, temp[W]);
                temp[K] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[K].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, temp[K]);
                temp[Omega] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FW[Omega].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, temp[Omega]);
                temp[T] = new RobinBC(1.0, 1.0, 0.0, 1.0);
                FW[T].getBoundary().getBoundaryRegion().setType(X, ye, Z, 3, temp[T]);
            }
        }
        //</editor-fold>
        //Z边界
        //<editor-fold>
        for (int Y = 0; Y <= WaterMesh.getNY() + 1; Y++) {
            for (int X = 0; X <= WaterMesh.getNX() + 1; X++) {
                //zw
                temp[W] = new RobinBC(0, 1.0, 0, 1.0);
                FW[W].getBoundary().getBoundaryRegion().setType(X, Y, zw, 7, temp[W]);
                temp[K] = new RobinBC(0, 1.0, 0, 1.0);
                FW[K].getBoundary().getBoundaryRegion().setType(X, Y, zw, 7, temp[K]);
                temp[Omega] = new RobinBC(0, 1.0, 0, 1.0);
                FW[Omega].getBoundary().getBoundaryRegion().setType(X, Y, zw, 7, temp[Omega]);
                temp[T] = new RobinBC(0, 1.0, 0, 1.0);
                FW[T].getBoundary().getBoundaryRegion().setType(X, Y, zw, 7, temp[T]);
                //ze
                temp[W] = new RobinBC(0, 1.0, 0, 1.0);
                FW[W].getBoundary().getBoundaryRegion().setType(X, Y, ze, 7, temp[W]);
                temp[K] = new RobinBC(0, 1.0, 0, 1.0);
                FW[K].getBoundary().getBoundaryRegion().setType(X, Y, ze, 7, temp[K]);
                temp[Omega] = new RobinBC(0, 1.0, 0, 1.0);
                FW[Omega].getBoundary().getBoundaryRegion().setType(X, Y, ze, 7, temp[Omega]);
                temp[T] = new RobinBC(0, 1.0, 0, 1.0);
                FW[T].getBoundary().getBoundaryRegion().setType(X, Y, ze, 7, temp[T]);
            }
        }
        //</editor-fold>
        xw = 0;
        xe = OilMesh.getNX() + 1;
        yw = 0;
        ye = OilMesh.getNY() + 1;
        zw = 0;
        ze = OilMesh.getNZ() + 1;
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
                temp[T] = new RobinBC(1.0, 1.0, 0.0, 1.0);
                FO[T].getBoundary().getBoundaryRegion().setType(xw, Y, Z, 3, temp[T]);
                //xe
                temp[W] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[W].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 1, temp[W]);
                temp[K] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[K].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 1, temp[K]);
                temp[Omega] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[Omega].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 1, temp[Omega]);
                temp[T] = new RobinBC(1.0, 1.0, 0.0, 1.0);
                FO[T].getBoundary().getBoundaryRegion().setType(xe, Y, Z, 3, temp[T]);
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
                temp[T] = new RobinBC(1.0, 1.0, 0.0, 1.0);
                FO[T].getBoundary().getBoundaryRegion().setType(X, yw, Z, 3, temp[T]);
                //ye
                temp[W] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[W].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, temp[W]);
                temp[K] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[K].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, temp[K]);
                temp[Omega] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[Omega].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, temp[Omega]);
                temp[T] = new RobinBC(0.0, 1.0, 0.0, 1.0);
                FO[T].getBoundary().getBoundaryRegion().setType(X, ye, Z, 1, temp[T]);
            }
        }
        //</editor-fold>
        //Z边界
        //<editor-fold>
        for (int Y = 0; Y <= OilMesh.getNY() + 1; Y++) {
            for (int X = 0; X <= OilMesh.getNX() + 1; X++) {
                //zw
                temp[W] = new RobinBC(0, 1.0, 0, 1.0);
                FO[W].getBoundary().getBoundaryRegion().setType(X, Y, zw, 7, temp[W]);
                temp[K] = new RobinBC(0, 1.0, 0, 1.0);
                FO[K].getBoundary().getBoundaryRegion().setType(X, Y, zw, 7, temp[K]);
                temp[Omega] = new RobinBC(0, 1.0, 0, 1.0);
                FO[Omega].getBoundary().getBoundaryRegion().setType(X, Y, zw, 7, temp[Omega]);
                temp[T] = new RobinBC(0, 1.0, 0, 1.0);
                FO[T].getBoundary().getBoundaryRegion().setType(X, Y, zw, 7, temp[T]);
                //ze
                temp[W] = new RobinBC(0, 1.0, 0, 1.0);
                FO[W].getBoundary().getBoundaryRegion().setType(X, Y, ze, 7, temp[W]);
                temp[K] = new RobinBC(0, 1.0, 0, 1.0);
                FO[K].getBoundary().getBoundaryRegion().setType(X, Y, ze, 7, temp[K]);
                temp[Omega] = new RobinBC(0, 1.0, 0, 1.0);
                FO[Omega].getBoundary().getBoundaryRegion().setType(X, Y, ze, 7, temp[Omega]);
                temp[T] = new RobinBC(0, 1.0, 0, 1.0);
                FO[T].getBoundary().getBoundaryRegion().setType(X, Y, ze, 7, temp[T]);
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
        System.out.println("begin LiquidMesh");
        WaterMesh = new BiSingleV2(12, bm.theta3 - bm.theta2,
                bm.getNX(), bm.getNY() - bm.interhL);
        WaterMesh.theta1 = bm.theta2;
        WaterMesh.theta2 = bm.theta3;
        WaterMesh.a = bm.a;
        WaterMesh.radia = bm.radia;
        WaterMesh.blockMesh();
        return true;
    }

    Boolean Log5 = Boolean.FALSE;

    public boolean Solver(int i) throws FileNotFoundException {
        if (i == 3) {
            scheme3();
        }
        if (i == 4) {
            scheme4();
        }
        if (i == 5) {
            scheme5();
        }
        return true;
    }

    /**
     * 单相双流体层流方程求解
     */
    void scheme3() throws FileNotFoundException {
        double resw = 0, reso = 0;
        int count = 0;
        double Qw, Qo;
        RESULT(count);
//        进入大循环
        do {
            count++;
            System.out.println("计算速度场W");
            cacW(dpdz, WaterMesh, FW[W], porW[mu].getF(),
                    porW[rho].getF(), porW[mut].getF(), sysWater);
            cacW(dpdz, OilMesh, FO[W], porO[mu].getF(),
                    porO[rho].getF(), porO[mut].getF(), sysOil);
            updateInterfaceW();
            resw = getRes(FW[W].getNewField(), FW[W].getOldField());
            reso = getRes(FO[W].getNewField(), FO[W].getOldField());
            FW[W].copyNew2Old();
            FO[W].copyNew2Old();
            System.out.println("resw=" + resw + " reso=" + reso);
            Qw = Tool.getTotalVol(FW[W], WaterMesh);
            Qo = Tool.getTotalVol(FO[W], OilMesh);
            System.out.println("计算总流量=Qw" + Qw);
            System.out.println("计算总流量=Qo" + Qo);
            System.out.println("总流量=" + (Qw + Qo));
        } while (max(resw, reso) > 1e-3 || count < 500);
        RESULT(count);
    }

    /**
     * 双流体方程+湍流模型求解
     */
    void scheme4() throws FileNotFoundException {
        double resw = 0, reso = 0;
        double reswt = 0, resot = 0;
        int count = 0;
        RESULT(count);
//      假定速度场
        cacW(dpdz, WaterMesh, FW[W], porW[mu].getF(),
                porW[rho].getF(), porW[mut].getF(), sysWater);
        cacW(dpdz, OilMesh, FO[W], porO[mu].getF(),
                porO[rho].getF(), porO[mut].getF(), sysOil);
        FW[W].copyNew2Old();
        FO[W].copyNew2Old();
        System.out.println("初始化湍流动能！");
        setK(FO[K], FO[W]);
        FO[K].copyNew2Old();
        setK(FW[K], FW[W]);
        FW[K].copyNew2Old();
        System.out.println("初始化湍流能量耗散率！");
        setOmega(FO[Omega], FO[K], OilMesh);
        FO[Omega].copyNew2Old();
        setOmega(FW[Omega], FW[K], WaterMesh);
        FW[Omega].copyNew2Old();
//        进入大循环
        do {
            count++;
            Turbulence.cacScenter(FO[W], OilMesh, porO[S].getF());
            Turbulence.cacScenter(FW[W], WaterMesh, porW[S].getF());
            //<editor-fold>湍流粘度
            System.out.println("计算湍流粘度！");
            cacMut(porO[rho].getF(), porO[mut].getF(), FO[K], FO[Omega]);
            cacMut(porW[rho].getF(), porW[mut].getF(), FW[K], FW[Omega]);
            //</editor-fold>
            //<editor-fold>湍流动能
            System.out.println("计算湍流动能！");
            cacKOil();
            cacKWater();
            FO[K].copyNew2Old();
            FW[K].copyNew2Old();
            //</editor-fold>
            //<editor-fold>湍流Omega
            System.out.println("计算Omega！");
            cacOmegaOil();
            cacOmegaWater();
            resot = getRes(FO[Omega].getNewField(), FO[Omega].getOldField());
            reswt = getRes(FW[Omega].getNewField(), FW[Omega].getOldField());
            System.out.println("湍流残差resot = " + resot + "湍流残差reswt = " + reswt);
            FO[Omega].copyNew2Old();
            FW[Omega].copyNew2Old();
            //</editor-fold>
            //<editor-fold>计算速度
            System.out.println("计算速度场W");
            cacW(dpdz, WaterMesh, FW[W], porW[mu].getF(),
                    porW[rho].getF(), porW[mut].getF(), sysWater);
            cacW(dpdz, OilMesh, FO[W], porO[mu].getF(),
                    porO[rho].getF(), porO[mut].getF(), sysOil);
            updateInterfaceW();
            resw = getRes(FW[W].getNewField(), FW[W].getOldField());
            reso = getRes(FO[W].getNewField(), FO[W].getOldField());
            FW[W].copyNew2Old();
            FO[W].copyNew2Old();
            System.out.println("resw=" + resw + " reso=" + reso);
            //</editor-fold>
        } while (max(resw, reso) > 1e-3 || count < 500);
        RESULT(count);
    }

    /**
     * 双流体方程+湍流模型求解+界面摩擦
     */
    void scheme5() throws FileNotFoundException {
        double resw = 0, reso = 0;
        double reswt = 0, resot = 0;
        int count = 0;
        RESULT(count);
//      假定速度场
        cacW(dpdz, WaterMesh, FW[W], porW[mu].getF(),
                porW[rho].getF(), porW[mut].getF(), sysWater);
        cacW(dpdz, OilMesh, FO[W], porO[mu].getF(),
                porO[rho].getF(), porO[mut].getF(), sysOil);
        FW[W].copyNew2Old();
        FO[W].copyNew2Old();
        System.out.println("初始化湍流动能！");
        setK(FO[K], FO[W]);
        FO[K].copyNew2Old();
        setK(FW[K], FW[W]);
        FW[K].copyNew2Old();
        System.out.println("初始化湍流能量耗散率！");
        setOmega(FO[Omega], FO[K], OilMesh);
        FO[Omega].copyNew2Old();
        setOmega(FW[Omega], FW[K], WaterMesh);
        FW[Omega].copyNew2Old();
//        进入大循环
        do {
            count++;
            Turbulence.cacScenter1(FO[W], OilMesh, porO[S].getF());
            Turbulence.cacScenter1(FW[W], WaterMesh, porW[S].getF());
            //<editor-fold>湍流粘度
            System.out.println("计算湍流粘度！");
            cacMut(porO[rho].getF(), porO[mut].getF(), FO[K], FO[Omega]);
            cacMut(porW[rho].getF(), porW[mut].getF(), FW[K], FW[Omega]);
            //</editor-fold>
            cacVelocity();
            cacTauInter();
            //<editor-fold>湍流动能
            System.out.println("计算湍流动能！");
            cacKOil();
            cacKWater();
            FO[K].copyNew2Old();
            FW[K].copyNew2Old();
            //</editor-fold>
            //<editor-fold>湍流Omega
            System.out.println("计算Omega！");
            cacOmegaOil();
            cacOmegaWater();
            resot = getRes(FO[Omega].getNewField(), FO[Omega].getOldField());
            reswt = getRes(FW[Omega].getNewField(), FW[Omega].getOldField());
            System.out.println("湍流残差resot = " + resot + " 湍流残差reswt = " + reswt);
            FO[Omega].copyNew2Old();
            FW[Omega].copyNew2Old();
            //</editor-fold>
            //<editor-fold>计算速度
            System.out.println("计算速度场W");
            cacW(dpdz, WaterMesh, FW[W], porW[mu].getF(),
                    porW[rho].getF(), porW[mut].getF(), sysWater);
            cacW(dpdz, OilMesh, FO[W], porO[mu].getF(),
                    porO[rho].getF(), porO[mut].getF(), sysOil);
            updateInterfaceW();
            resw = getRes(FW[W].getNewField(), FW[W].getOldField());
            reso = getRes(FO[W].getNewField(), FO[W].getOldField());
            FW[W].copyNew2Old();
            FO[W].copyNew2Old();
            System.out.println("resw=" + resw + " reso=" + reso);
            //</editor-fold>
        } while (max(resw, reso) > 1e-3 && count < 10000);
        RESULT(count);
    }

    /**
     * 计算油相湍流动能场
     */
    void cacKOil() {
        double URF = 0.9;
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;
        double simga2 = 0.5;
        double betastar = 0.09;
        double GEN;
        double CTRANS = 11.63;
        double CAPPA = 0.41;
        Label flagK = new Label("KOil");
        for (int Z = 1; Z <= OilMesh.getNZ(); ++Z) {
            for (int Y = 1; Y <= OilMesh.getNY(); ++Y) {
                for (int X = 1; X <= OilMesh.getNX(); ++X) {
                    //<editor-fold>
                    flagK.setFlag(FO[K].getBoundary(), X, Y);
                    double vol
                            = OilMesh.getDXU()[X] * OilMesh.getDYV()[Y]
                            * OilMesh.J(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]);
                    Fw = 0;
                    Fe = 0;
                    Fs = 0;
                    Fn = 0;
//
                    mup = (porO[mu].getF()[X][Y][Z] + porO[mut].getF()[X][Y][Z] * simga2);
                    muw = (porO[mu].getF()[X - 1][Y][Z] + porO[mut].getF()[X - 1][Y][Z] * simga2);
                    mue = (porO[mu].getF()[X + 1][Y][Z] + porO[mut].getF()[X + 1][Y][Z] * simga2);
                    mus = (porO[mu].getF()[X][Y - 1][Z] + porO[mut].getF()[X][Y - 1][Z] * simga2);
                    mun = (porO[mu].getF()[X][Y + 1][Z] + porO[mut].getF()[X][Y + 1][Z] * simga2);
//得到gamma
                    gamw = faceValue(muw, mup - muw, OilMesh.gettPx()[X - 1], OilMesh.getTx()[X - 1], OilMesh.gettPx()[X]);
                    game = faceValue(mup, mue - mup, OilMesh.gettPx()[X], OilMesh.getTx()[X], OilMesh.gettPx()[X + 1]);
                    gams = faceValue(mus, mup - mus, OilMesh.gettPy()[Y - 1], OilMesh.getTy()[Y - 1], OilMesh.gettPy()[Y]);
                    gamn = faceValue(mup, mun - mup, OilMesh.gettPy()[Y], OilMesh.getTy()[Y], OilMesh.gettPy()[Y + 1]);
//GET d
                    Dw = gamw * OilMesh.getDYV()[Y] / OilMesh.getDXP()[X - 1];
                    De = game * OilMesh.getDYV()[Y] / OilMesh.getDXP()[X];
                    Ds = gams * OilMesh.getDXU()[X] / OilMesh.getDYP()[Y - 1];
                    Dn = gamn * OilMesh.getDXU()[X] / OilMesh.getDYP()[Y];
//系数
                    coeO[K].aw[X][Y][Z] = (1 - flagK.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeO[K].ae[X][Y][Z] = (1 - flagK.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeO[K].as[X][Y][Z] = (1 - flagK.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeO[K].an[X][Y][Z] = (1 - flagK.n) * (Math.max(-Fn, 0.0) + Dn);//n
//设定附加源项
                    Sp = -betastar * porO[rho].getF()[X][Y][Z] * FO[Omega].getNewField()[X][Y][Z];
                    Spad
                            = Spadw(FO[K].getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, OilMesh.getDXP()[X - 1], vol,
                                    FO[K].getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Spade(FO[K].getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, OilMesh.getDXP()[X], vol,
                                    FO[K].getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Spadw(FO[K].getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, OilMesh.getDYP()[Y - 1], vol,
                                    FO[K].getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Spade(FO[K].getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, OilMesh.getDYP()[Y], vol,
                                    FO[K].getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
//                    得到ap系数
                    coeO[K].ap[X][Y][Z]
                            = coeO[K].aw[X][Y][Z] + coeO[K].ae[X][Y][Z]
                            + coeO[K].as[X][Y][Z] + coeO[K].an[X][Y][Z]
                            + (Fe - Fw) + (Fn - Fs) - (Sp + Spad) * vol;

//                    得到源项
                    //<editor-fold>
                    if ((flagK.w == 1) || (flagK.e == 1)
                            || (flagK.s == 1) || (flagK.n == 1)) {
                        double VISS = porO[mu].getF()[X][Y][Z];
                        double YPL, DN = 0;
                        double TAU = 0;
                        double CMU = 0.09;
                        double CMU25 = sqrt(sqrt(CMU));
                        double ELOG = 8.342;
                        double CK = CMU25 * sqrt(max(0, FO[K].getOldField()[X][Y][Z]));
                        double VISCW, VISW;
                        double YPL1, epl = 0, es, Uint;
                        if (flagK.w == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realY(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realX(OilMesh.gettPx()[X - 1], OilMesh.gettPy()[Y]),
                                    OilMesh.realY(OilMesh.gettPx()[X - 1], OilMesh.gettPy()[Y]));
                            YPL = porO[rho].getF()[X][Y][Z] * CK * DN / porO[mu].getF()[X][Y][Z];
                            VISCW = YPL * porO[mu].getF()[X][Y][Z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(porO[mu].getF()[X][Y][Z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (FO[W].getNewField()[X][Y][Z] - FO[W].getNewField()[X - 1][Y][Z]) / DN;
                        }
                        if (flagK.e == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realY(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realX(OilMesh.gettPx()[X + 1], OilMesh.gettPy()[Y]),
                                    OilMesh.realY(OilMesh.gettPx()[X + 1], OilMesh.gettPy()[Y]));
                            YPL = porO[rho].getF()[X][Y][Z] * CK * DN / porO[mu].getF()[X][Y][Z];
                            VISCW = YPL * porO[mu].getF()[X][Y][Z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(porO[mu].getF()[X][Y][Z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (FO[W].getNewField()[X][Y][Z] - FO[W].getNewField()[X + 1][Y][Z]) / DN;
                        }
                        if (flagK.s == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realY(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realX(OilMesh.gettPx()[X], OilMesh.gettPy()[Y - 1]),
                                    OilMesh.realY(OilMesh.gettPx()[X], OilMesh.gettPy()[Y - 1]));
                            YPL = porO[rho].getF()[X][Y][Z] * CK * DN / porO[mu].getF()[X][Y][Z];
                            VISCW = YPL * porO[mu].getF()[X][Y][Z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(porO[mu].getF()[X][Y][Z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (FO[W].getNewField()[X][Y][Z] - FO[W].getNewField()[X][Y - 1][Z]) / DN;
                        }
                        if (flagK.n == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realY(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realX(OilMesh.gettPx()[X], OilMesh.gettPy()[Y + 1]),
                                    OilMesh.realY(OilMesh.gettPx()[X], OilMesh.gettPy()[Y + 1]));
                            YPL = porO[rho].getF()[X][Y][Z] * CK * DN / porO[mu].getF()[X][Y][Z];
                            VISCW = YPL * porO[mu].getF()[X][Y][Z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(porO[mu].getF()[X][Y][Z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (FO[W].getNewField()[X][Y][Z] - FO[W].getNewField()[X][Y + 1][Z]) / DN;
                        }
                        GEN
                                = abs(TAU) * CMU25 * sqrt(max(0, FO[K].getOldField()[X][Y][Z]))
                                / (DN * CAPPA);
                    } else {
                        GEN = (porO[mut].getF()[X][Y][Z] * porO[S].getF()[X][Y][Z]);
                    }
                    //</editor-fold>
                    Sc = GEN;
                    Scad
                            = Scadw(FO[K].getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, OilMesh.getDXP()[X - 1], vol,
                                    FO[K].getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Scade(FO[K].getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, OilMesh.getDXP()[X], vol,
                                    FO[K].getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Scadw(FO[K].getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, OilMesh.getDYP()[Y - 1], vol,
                                    FO[K].getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Scade(FO[K].getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, OilMesh.getDYP()[Y], vol,
                                    FO[K].getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
                    coeO[K].b[X][Y][Z] = (Sc + Scad) * vol;
                    //</editor-fold>
//Under-relaxizaion
                    coeO[K].b[X][Y][Z]
                            = coeO[K].b[X][Y][Z] + (URFFI - 1.0) * coeO[K].ap[X][Y][Z]
                            * FO[K].getOldField()[X][Y][Z];
                    coeO[K].ap[X][Y][Z] = (coeO[K].ap[X][Y][Z]) * URFFI;
                }
            }
        }
        FO[K].getLog().rl = MatrixSolver.gaussSeidel(FO[K].getNewField(), coeO[K], sysOil);
        FO[K].getLog().log("Koil");
    }

    /**
     * 计算水相湍流动能场
     */
    void cacKWater() {
        double URF = 0.9;
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;
        double simga2 = 0.5;
        double betastar = 0.09;
        double GEN;
        double CTRANS = 11.63;
        double CAPPA = 0.41;
        Label flagK = new Label("KOil");
        for (int Z = 1; Z <= WaterMesh.getNZ(); ++Z) {
            for (int Y = 1; Y <= WaterMesh.getNY(); ++Y) {
                for (int X = 1; X <= WaterMesh.getNX(); ++X) {
                    //<editor-fold>
                    flagK.setFlag(FW[K].getBoundary(), X, Y);
                    double vol
                            = WaterMesh.getDXU()[X] * WaterMesh.getDYV()[Y]
                            * WaterMesh.J(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]);
                    Fw = 0;
                    Fe = 0;
                    Fs = 0;
                    Fn = 0;
//
                    mup = (porW[mu].getF()[X][Y][Z] + porW[mut].getF()[X][Y][Z] * simga2);
                    muw = (porW[mu].getF()[X - 1][Y][Z] + porW[mut].getF()[X - 1][Y][Z] * simga2);
                    mue = (porW[mu].getF()[X + 1][Y][Z] + porW[mut].getF()[X + 1][Y][Z] * simga2);
                    mus = (porW[mu].getF()[X][Y - 1][Z] + porW[mut].getF()[X][Y - 1][Z] * simga2);
                    mun = (porW[mu].getF()[X][Y + 1][Z] + porW[mut].getF()[X][Y + 1][Z] * simga2);
//得到gamma
                    gamw = faceValue(muw, mup - muw, WaterMesh.gettPx()[X - 1], WaterMesh.getTx()[X - 1], WaterMesh.gettPx()[X]);
                    game = faceValue(mup, mue - mup, WaterMesh.gettPx()[X], WaterMesh.getTx()[X], WaterMesh.gettPx()[X + 1]);
                    gams = faceValue(mus, mup - mus, WaterMesh.gettPy()[Y - 1], WaterMesh.getTy()[Y - 1], WaterMesh.gettPy()[Y]);
                    gamn = faceValue(mup, mun - mup, WaterMesh.gettPy()[Y], WaterMesh.getTy()[Y], WaterMesh.gettPy()[Y + 1]);
//GET d
                    Dw = gamw * WaterMesh.getDYV()[Y] / WaterMesh.getDXP()[X - 1];
                    De = game * WaterMesh.getDYV()[Y] / WaterMesh.getDXP()[X];
                    Ds = gams * WaterMesh.getDXU()[X] / WaterMesh.getDYP()[Y - 1];
                    Dn = gamn * WaterMesh.getDXU()[X] / WaterMesh.getDYP()[Y];
//系数
                    coeW[K].aw[X][Y][Z] = (1 - flagK.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeW[K].ae[X][Y][Z] = (1 - flagK.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeW[K].as[X][Y][Z] = (1 - flagK.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeW[K].an[X][Y][Z] = (1 - flagK.n) * (Math.max(-Fn, 0.0) + Dn);//n
//设定附加源项
                    Sp = -betastar * porW[rho].getF()[X][Y][Z] * FW[Omega].getNewField()[X][Y][Z];
                    Spad
                            = Spadw(FW[K].getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, WaterMesh.getDXP()[X - 1], vol,
                                    FW[K].getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Spade(FW[K].getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, WaterMesh.getDXP()[X], vol,
                                    FW[K].getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Spadw(FW[K].getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, WaterMesh.getDYP()[Y - 1], vol,
                                    FW[K].getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Spade(FW[K].getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, WaterMesh.getDYP()[Y], vol,
                                    FW[K].getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
//                    得到ap系数
                    coeW[K].ap[X][Y][Z]
                            = coeW[K].aw[X][Y][Z] + coeW[K].ae[X][Y][Z]
                            + coeW[K].as[X][Y][Z] + coeW[K].an[X][Y][Z]
                            + (Fe - Fw) + (Fn - Fs) - (Sp + Spad) * vol;
//                    得到源项
//<editor-fold>
                    if ((flagK.w == 1) || (flagK.e == 1)
                            || (flagK.s == 1) || (flagK.n == 1)) {
                        double VISS = porW[mu].getF()[X][Y][Z];
                        double YPL, DN = 0;
                        double TAU = 0;
                        double CMU = 0.09;
                        double CMU25 = sqrt(sqrt(CMU));
                        double ELOG = 8.342;
                        double CK = CMU25 * sqrt(max(0, FW[K].getOldField()[X][Y][Z]));
                        double VISCW, VISW;
                        double YPL1, epl = 0, es, Uint;
                        if (flagK.w == 1) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[X - 1], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X - 1], WaterMesh.gettPy()[Y]));
                            YPL = porW[rho].getF()[X][Y][Z] * CK * DN / porW[mu].getF()[X][Y][Z];
                            VISCW = YPL * porW[mu].getF()[X][Y][Z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(porW[mu].getF()[X][Y][Z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (FW[W].getNewField()[X][Y][Z] - FW[W].getNewField()[X - 1][Y][Z]) / DN;
                        }
                        if (flagK.e == 1) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[X + 1], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X + 1], WaterMesh.gettPy()[Y]));
                            YPL = porW[rho].getF()[X][Y][Z] * CK * DN / porW[mu].getF()[X][Y][Z];
                            VISCW = YPL * porW[mu].getF()[X][Y][Z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(porW[mu].getF()[X][Y][Z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (FW[W].getNewField()[X][Y][Z] - FW[W].getNewField()[X + 1][Y][Z]) / DN;
                        }
                        if (flagK.s == 1) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y - 1]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y - 1]));
                            YPL = porW[rho].getF()[X][Y][Z] * CK * DN / porW[mu].getF()[X][Y][Z];
//
                            VISCW = YPL * porW[mu].getF()[X][Y][Z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(porW[mu].getF()[X][Y][Z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (FW[W].getNewField()[X][Y][Z] - FW[W].getNewField()[X][Y - 1][Z]) / DN;
                        }
                        if (flagK.n == 1) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y + 1]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y + 1]));
                            YPL = porW[rho].getF()[X][Y][Z] * CK * DN / porW[mu].getF()[X][Y][Z];
                            VISCW = YPL * porW[mu].getF()[X][Y][Z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(porW[mu].getF()[X][Y][Z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (FW[W].getNewField()[X][Y][Z] - FW[W].getNewField()[X][Y + 1][Z]) / DN;
                        }
                        GEN
                                = abs(TAU) * CMU25 * sqrt(max(0, FW[K].getOldField()[X][Y][Z]))
                                / (DN * CAPPA);
                    } else {
                        GEN = (porW[mut].getF()[X][Y][Z] * porW[S].getF()[X][Y][Z]);
                    }
//</editor-fold>
                    Sc = GEN;
                    Scad
                            = Scadw(FW[K].getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, WaterMesh.getDXP()[X - 1], vol,
                                    FW[K].getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Scade(FW[K].getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, WaterMesh.getDXP()[X], vol,
                                    FW[K].getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Scadw(FW[K].getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, WaterMesh.getDYP()[Y - 1], vol,
                                    FW[K].getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Scade(FW[K].getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, WaterMesh.getDYP()[Y], vol,
                                    FW[K].getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
                    coeW[K].b[X][Y][Z] = (Sc + Scad) * vol;
                    //</editor-fold>
//Under-relaxizaion
                    coeW[K].b[X][Y][Z]
                            = coeW[K].b[X][Y][Z] + (URFFI - 1.0) * coeW[K].ap[X][Y][Z]
                            * FW[K].getOldField()[X][Y][Z];
                    coeW[K].ap[X][Y][Z] = coeW[K].ap[X][Y][Z] * URFFI;
                }
            }
        }
        FW[K].getLog().rl = MatrixSolver.gaussSeidel(FW[K].getNewField(), coeW[K], sysWater);
        FW[K].getLog().log("Kwater");
    }

    /**
     * 计算油相Omega
     */
    void cacOmegaOil() {
        double URF = 0.9;
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;
        double sigma1 = 0.5;
        double alpha = 0.555;
        double beta = 0.075;
        Label flagW = new Label("OmegaOil");
        for (int Z = 1; Z <= OilMesh.getNZ(); ++Z) {
            for (int Y = 1; Y <= OilMesh.getNY(); ++Y) {
                for (int X = 1; X <= OilMesh.getNX(); ++X) {
                    //<editor-fold>
                    flagW.setFlag(FO[Omega].getBoundary(), X, Y);
                    double vol
                            = OilMesh.getDXU()[X] * OilMesh.getDYV()[Y]
                            * OilMesh.J(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]);
                    Fw = 0;
                    Fe = 0;
                    Fs = 0;
                    Fn = 0;
//
                    mup = (porO[mu].getF()[X][Y][Z] + porO[mut].getF()[X][Y][Z] * sigma1);
                    muw = (porO[mu].getF()[X - 1][Y][Z] + porO[mut].getF()[X - 1][Y][Z] * sigma1);
                    mue = (porO[mu].getF()[X + 1][Y][Z] + porO[mut].getF()[X + 1][Y][Z] * sigma1);
                    mus = (porO[mu].getF()[X][Y - 1][Z] + porO[mut].getF()[X][Y - 1][Z] * sigma1);
                    mun = (porO[mu].getF()[X][Y + 1][Z] + porO[mut].getF()[X][Y + 1][Z] * sigma1);
//得到gamma
                    gamw = faceValue(muw, mup - muw, OilMesh.gettPx()[X - 1], OilMesh.getTx()[X - 1], OilMesh.gettPx()[X]);
                    game = faceValue(mup, mue - mup, OilMesh.gettPx()[X], OilMesh.getTx()[X], OilMesh.gettPx()[X + 1]);
                    gams = faceValue(mus, mup - mus, OilMesh.gettPy()[Y - 1], OilMesh.getTy()[Y - 1], OilMesh.gettPy()[Y]);
                    gamn = faceValue(mup, mun - mup, OilMesh.gettPy()[Y], OilMesh.getTy()[Y], OilMesh.gettPy()[Y + 1]);
//GET d
                    Dw = gamw * OilMesh.getDYV()[Y] / OilMesh.getDXP()[X - 1];
                    De = game * OilMesh.getDYV()[Y] / OilMesh.getDXP()[X];
                    Ds = gams * OilMesh.getDXU()[X] / OilMesh.getDYP()[Y - 1];
                    Dn = gamn * OilMesh.getDXU()[X] / OilMesh.getDYP()[Y];
//系数
                    coeO[Omega].aw[X][Y][Z] = (1 - flagW.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeO[Omega].ae[X][Y][Z] = (1 - flagW.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeO[Omega].as[X][Y][Z] = (1 - flagW.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeO[Omega].an[X][Y][Z] = (1 - flagW.n) * (Math.max(-Fn, 0.0) + Dn);//n
//设定附加源项
                    Sp = -FO[Omega].getOldField()[X][Y][Z] * porO[rho].getF()[X][Y][Z] * beta;
                    Spad
                            = Spadw(FO[Omega].getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, OilMesh.getDXP()[X - 1], vol,
                                    FO[Omega].getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Spade(FO[Omega].getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, OilMesh.getDXP()[X], vol,
                                    FO[Omega].getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Spadw(FO[Omega].getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, OilMesh.getDYP()[Y - 1], vol,
                                    FO[Omega].getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Spade(FO[Omega].getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, OilMesh.getDYP()[Y], vol,
                                    FO[Omega].getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
//得到ap系数
                    coeO[Omega].ap[X][Y][Z]
                            = coeO[Omega].aw[X][Y][Z] + coeO[Omega].ae[X][Y][Z]
                            + coeO[Omega].as[X][Y][Z] + coeO[Omega].an[X][Y][Z]
                            + (Fe - Fw) + (Fn - Fs) - (Sp + Spad) * vol;
//得到源项
                    Sc
                            = alpha * (porO[mut].getF()[X][Y][Z] * porO[S].getF()[X][Y][Z])
                            * FO[Omega].getOldField()[X][Y][Z] / FO[K].getOldField()[X][Y][Z];
                    Scad
                            = Scadw(FO[Omega].getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, OilMesh.getDXP()[X - 1], vol,
                                    FO[Omega].getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Scade(FO[Omega].getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, OilMesh.getDXP()[X], vol,
                                    FO[Omega].getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Scadw(FO[Omega].getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, OilMesh.getDYP()[Y - 1], vol,
                                    FO[Omega].getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Scade(FO[Omega].getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, OilMesh.getDYP()[Y], vol,
                                    FO[Omega].getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
                    coeO[Omega].b[X][Y][Z] = (Sc + Scad) * vol;
                    //</editor-fold>
                }
            }
        }
        double DN = 1E-30;
        double YPL1, YP, epl = 0, es, Uint = 1.0;
        for (int Z = 1; Z <= OilMesh.getNZ(); ++Z) {
            for (int Y = 1; Y <= OilMesh.getNY(); ++Y) {
                for (int X = 1; X <= OilMesh.getNX(); ++X) {
                    //<editor-fold>
                    flagW.setFlag(FO[Omega].getBoundary(), X, Y);
                    if ((flagW.w == 1) || (flagW.e == 1)
                            || (flagW.s == 1) || (flagW.n == 1)) {
                        if (flagW.w == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realY(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realX(OilMesh.gettPx()[X - 1], OilMesh.gettPy()[Y]),
                                    OilMesh.realY(OilMesh.gettPx()[X - 1], OilMesh.gettPy()[Y]));
                        }
                        if (flagW.e == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realY(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realX(OilMesh.gettPx()[X + 1], OilMesh.gettPy()[Y]),
                                    OilMesh.realY(OilMesh.gettPx()[X + 1], OilMesh.gettPy()[Y]));
                        }
                        if (flagW.s == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realY(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realX(OilMesh.gettPx()[X], OilMesh.gettPy()[Y - 1]),
                                    OilMesh.realY(OilMesh.gettPx()[X], OilMesh.gettPy()[Y - 1]));

                        }
                        if (flagW.n == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realY(OilMesh.gettPx()[X], OilMesh.gettPy()[Y]),
                                    OilMesh.realX(OilMesh.gettPx()[X], OilMesh.gettPy()[Y + 1]),
                                    OilMesh.realY(OilMesh.gettPx()[X], OilMesh.gettPy()[Y + 1]));

                        }
                        coeO[Omega].ap[X][Y][Z] = 1.0;
                        coeO[Omega].aw[X][Y][Z] = 0;//w
                        coeO[Omega].ae[X][Y][Z] = 0;//e
                        coeO[Omega].as[X][Y][Z] = 0;//s
                        coeO[Omega].an[X][Y][Z] = 0;//n
                        coeO[Omega].b[X][Y][Z]
                                = 6.0 * porO[mu].getF()[X][Y][Z]
                                / (porO[rho].getF()[X][Y][Z] * beta * DN * DN);
                    }
                    //</editor-fold>
                    //under-relaxization
                    coeO[Omega].b[X][Y][Z]
                            = coeO[Omega].b[X][Y][Z] + (URFFI - 1.0) * coeO[Omega].ap[X][Y][Z]
                            * FO[Omega].getOldField()[X][Y][Z];
                    coeO[Omega].ap[X][Y][Z] = (coeO[Omega].ap[X][Y][Z]) * URFFI;
                }
            }
        }
        FO[Omega].getLog().rl = MatrixSolver.gaussSeidel(FO[Omega].getNewField(), coeO[Omega], sysOil);
        FO[Omega].getLog().log("OmegaOil");
    }

    /**
     * 计算水相Omega
     */
    void cacOmegaWater() {
        double URF = 0.9;
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;
        double sigma1 = 0.5;
        double alpha = 0.555;
        double beta = 0.075;
        Label flagW = new Label("OmegaWater");
        for (int Z = 1; Z < FW[Omega].getNewField()[0][0].length - 1; ++Z) {
            for (int Y = 1; Y < FW[Omega].getNewField()[0].length - 1; ++Y) {
                for (int X = 1; X < FW[Omega].getNewField().length - 1; ++X) {
                    //<editor-fold>
                    flagW.setFlag(FW[Omega].getBoundary(), X, Y);
                    double vol
                            = WaterMesh.getDXU()[X] * WaterMesh.getDYV()[Y]
                            * WaterMesh.J(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]);
                    Fw = 0;
                    Fe = 0;
                    Fs = 0;
                    Fn = 0;
//
                    mup = (porW[mu].getF()[X][Y][Z] + porW[mut].getF()[X][Y][Z] * sigma1);
                    muw = (porW[mu].getF()[X - 1][Y][Z] + porW[mut].getF()[X - 1][Y][Z] * sigma1);
                    mue = (porW[mu].getF()[X + 1][Y][Z] + porW[mut].getF()[X + 1][Y][Z] * sigma1);
                    mus = (porW[mu].getF()[X][Y - 1][Z] + porW[mut].getF()[X][Y - 1][Z] * sigma1);
                    mun = (porW[mu].getF()[X][Y + 1][Z] + porW[mut].getF()[X][Y + 1][Z] * sigma1);
//得到gamma
                    gamw = faceValue(muw, mup - muw, WaterMesh.gettPx()[X - 1], WaterMesh.getTx()[X - 1], WaterMesh.gettPx()[X]);
                    game = faceValue(mup, mue - mup, WaterMesh.gettPx()[X], WaterMesh.getTx()[X], WaterMesh.gettPx()[X + 1]);
                    gams = faceValue(mus, mup - mus, WaterMesh.gettPy()[Y - 1], WaterMesh.getTy()[Y - 1], WaterMesh.gettPy()[Y]);
                    gamn = faceValue(mup, mun - mup, WaterMesh.gettPy()[Y], WaterMesh.getTy()[Y], WaterMesh.gettPy()[Y + 1]);
//GET d
                    Dw = gamw * WaterMesh.getDYV()[Y] / WaterMesh.getDXP()[X - 1];
                    De = game * WaterMesh.getDYV()[Y] / WaterMesh.getDXP()[X];
                    Ds = gams * WaterMesh.getDXU()[X] / WaterMesh.getDYP()[Y - 1];
                    Dn = gamn * WaterMesh.getDXU()[X] / WaterMesh.getDYP()[Y];
//系数
                    coeW[Omega].aw[X][Y][Z] = (1 - flagW.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeW[Omega].ae[X][Y][Z] = (1 - flagW.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeW[Omega].as[X][Y][Z] = (1 - flagW.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeW[Omega].an[X][Y][Z] = (1 - flagW.n) * (Math.max(-Fn, 0.0) + Dn);//n
//设定附加源项
                    Sp = -FW[Omega].getOldField()[X][Y][Z] * porW[rho].getF()[X][Y][Z] * beta;
                    Spad
                            = Spadw(FW[Omega].getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, WaterMesh.getDXP()[X - 1], vol,
                                    FW[Omega].getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Spade(FW[Omega].getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, WaterMesh.getDXP()[X], vol,
                                    FW[Omega].getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Spadw(FW[Omega].getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, WaterMesh.getDYP()[Y - 1], vol,
                                    FW[Omega].getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Spade(FW[Omega].getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, WaterMesh.getDYP()[Y], vol,
                                    FW[Omega].getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
//得到ap系数
                    coeW[Omega].ap[X][Y][Z]
                            = coeW[Omega].aw[X][Y][Z] + coeW[Omega].ae[X][Y][Z]
                            + coeW[Omega].as[X][Y][Z] + coeW[Omega].an[X][Y][Z]
                            + (Fe - Fw) + (Fn - Fs) - (Sp + Spad) * vol;
//得到源项
                    Sc
                            = alpha * (porW[mut].getF()[X][Y][Z] * porW[S].getF()[X][Y][Z])
                            * FW[Omega].getOldField()[X][Y][Z] / FW[K].getOldField()[X][Y][Z];
                    Scad
                            = Scadw(FW[Omega].getBoundary().getType(X - 1, Y, Z),
                                    Dw, Fw, WaterMesh.getDXP()[X - 1], vol,
                                    FW[Omega].getBoundary().getBoundaryRegion().parameter[X - 1][Y][Z])
                            + Scade(FW[Omega].getBoundary().getType(X + 1, Y, Z),
                                    De, Fe, WaterMesh.getDXP()[X], vol,
                                    FW[Omega].getBoundary().getBoundaryRegion().parameter[X + 1][Y][Z])
                            + Scadw(FW[Omega].getBoundary().getType(X, Y - 1, Z),
                                    Ds, Fs, WaterMesh.getDYP()[Y - 1], vol,
                                    FW[Omega].getBoundary().getBoundaryRegion().parameter[X][Y - 1][Z])
                            + Scade(FW[Omega].getBoundary().getType(X, Y + 1, Z),
                                    Dn, Fn, WaterMesh.getDYP()[Y], vol,
                                    FW[Omega].getBoundary().getBoundaryRegion().parameter[X][Y + 1][Z]);
                    coeW[Omega].b[X][Y][Z] = (Sc + Scad) * vol;
                    //</editor-fold>
                }
            }
        }
        double DN = 1E-30;
        double Uint = 1.0, es, epl = 0, YPL1, YP;
        for (int Z = 1; Z <= WaterMesh.getNZ(); ++Z) {
            for (int Y = 1; Y <= WaterMesh.getNY(); ++Y) {
                for (int X = 1; X <= WaterMesh.getNX(); ++X) {
                    //<editor-fold>
                    flagW.setFlag(FW[Omega].getBoundary(), X, Y);
                    if ((flagW.w == 1) || (flagW.e == 1)
                            || (flagW.s == 1) || (flagW.n == 1)) {
                        if (flagW.w == 1) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[X - 1], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X - 1], WaterMesh.gettPy()[Y]));
                        }
                        if (flagW.e == 1) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[X + 1], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X + 1], WaterMesh.gettPy()[Y]));
                        }
                        if (flagW.s == 1) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y - 1]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y - 1]));
                        }
                        if (flagW.n == 1) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y + 1]),
                                    WaterMesh.realY(WaterMesh.gettPx()[X], WaterMesh.gettPy()[Y + 1]));
                        }
                        coeW[Omega].ap[X][Y][Z] = 1.0;
                        coeW[Omega].aw[X][Y][Z] = 0;//w
                        coeW[Omega].ae[X][Y][Z] = 0;//e
                        coeW[Omega].as[X][Y][Z] = 0;//s
                        coeW[Omega].an[X][Y][Z] = 0;//n
                        coeW[Omega].b[X][Y][Z]
                                = 6.0 * porW[mu].getF()[X][Y][Z]
                                / (porW[rho].getF()[X][Y][Z] * beta * DN * DN);
                    }
                    //</editor-fold>
                    //under-relaxization
                    coeW[Omega].b[X][Y][Z]
                            = coeW[Omega].b[X][Y][Z] + (URFFI - 1.0) * coeW[Omega].ap[X][Y][Z]
                            * FW[Omega].getOldField()[X][Y][Z];
                    coeW[Omega].ap[X][Y][Z] = coeW[Omega].ap[X][Y][Z] * URFFI;
                }
            }
        }
        FW[Omega].getLog().rl = MatrixSolver.gaussSeidel(FW[Omega].getNewField(), coeW[Omega], sysWater);
        FW[Omega].getLog().log("OmegaWater");
    }

    /**
     * 计算速度场
     *
     * @param dpdz
     */
    void cacW(double dpdz, Mesh mesh, Field W,
            double[][][] mu, double[][][] rho,
            double[][][] Mut, SystemControl sys
    ) {
        double URF = 0.9;
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
     * 引入温度因素改变粘度
     */
    void cacMu() {

    }

    /**
     * 更新界面速度和边界条件
     */
    void updateInterfaceW() {
        double tau = abs(Tauo) < abs(Tauw) ? Tauw : Tauo;
        double tauo, tauw;
        if (wso > wsw) {
            tauo = tau;
            tauw = -tau;
        } else {
            tauo = -tau;
            tauw = tau;
        }
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
            double w = (porW[mu].getF()[x][1][1]) / yWater;
            double o = (porO[mu].getF()[x][OilMesh.getNY()][1]) / yOil;
            double A = (tauw) / w + FW[W].getNewField()[x][1][1];
            double B = (tauo) / o + FO[W].getNewField()[x][OilMesh.getNY()][1];
            double u = (A + B) / 2.0;
            u = 1.1 * cacWave(x);
            RobinBC robinbc = new RobinBC(0.0, 1.0, u, 1.0);
            FO[W].setNewField(x, OilMesh.getNY() + 1, 1, u);
            FO[W].getBoundary().getBoundaryRegion().setType(x, OilMesh.getNY() + 1, 1, 1, robinbc);
            FW[W].setNewField(x, 0, 1, u);
            FW[W].getBoundary().getBoundaryRegion().setType(x, 0, 1, 1, robinbc);
        }
    }

    double cacWave(int x) {
        double temp;
        double distance;
        distance = Tool.distance(
                0, 0,
                OilMesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY() + 1]),
                OilMesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY() + 1]));
        double alpha = 1.0;
        temp = max(0, (1.0 - distance / bm.a));
        temp = temp - 1.0;
        temp = alpha * (-temp * temp + 1.0);
        return temp;
    }

    double Tauo = 0, Tauw = 0;
    double eso = 0, esw = 0;
    double wso, wsw;

    /**
     * 求解界面的摩擦
     */
    void cacTauInter() {
        double fint;
        double leta;
        //water
        int y = 1, z = 1;
        for (int x = 1; x <= WaterMesh.getNX(); ++x) {
            leta = sqrt(WaterMesh.alpha(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]));
            fint = cacFintES(x, y, z);
            Tauw
                    = Tauw
                    + WaterMesh.getDXU()[x] * leta / (2.0 * WaterMesh.a)
                    * fint * porW[rho].getF()[x][y][z]
                    * (wso - wsw) * abs(wso - wsw)
                    / 2.0;
        }
        //oil
        y = OilMesh.getNY();
        for (int x = 1; x <= OilMesh.getNX(); ++x) {
            leta = sqrt(OilMesh.alpha(OilMesh.gettPx()[x], OilMesh.gettPy()[y]));
            fint = cacFintES(x, y, z);
            Tauo
                    = Tauo
                    + OilMesh.getDXU()[x] * leta / (2.0 * OilMesh.a)
                    * fint * porO[rho].getF()[x][y][z]
                    * (wso - wsw) * abs(wso - wsw)
                    / 2.0;
        }
    }

    /**
     * 计算表观速度
     */
    void cacVelocity() {
        wso = Tool.getTotalVol(FO[W], OilMesh) / (PI * this.R * this.R);
        wsw = Tool.getTotalVol(FW[W], WaterMesh) / (PI * this.R * this.R);
    }

    /**
     * 计算粗糙度
     *
     * @param Tau
     * @param rho
     * @param x
     * @param y
     * @param z
     * @return
     */
    double cacFintES(int x, int y, int z) {
        double hw = this.hl;
        double H = hw / this.R - 1.0;
        double Si = 2 * this.R * sqrt(1.0 - H * H);
        double So = 2 * this.R * Math.acos(H);
        double Sw = 2 * this.R - So;
        double Ao = this.R * 2.0 / 4.0 * (So - Si * H);
        double Aw = PI * this.R * this.R - Ao;
        double Dint, Reint;
//     判断快慢相   
//       <editor-fold>
        if (wso > wsw) {
            Dint = 4.0 * Ao / (So + Si);
            Reint = porO[rho].getF()[x][OilMesh.getNY()][z] * wso * Dint
                    / porO[mu].getF()[x][OilMesh.getNY()][z];
        } else if (wsw > wso) {
            Dint = 4.0 * Aw / (Sw + Si);
            Reint = porW[rho].getF()[x][1][z] * wsw * Dint
                    / porW[mu].getF()[x][1][z];
        } else {
            Dint = 4.0 * Aw / Sw;
            Reint = porW[rho].getF()[x][1][z] * wsw * Dint
                    / porW[mu].getF()[x][1][z];
        }
//       </editor-fold>
        double e = cacEs(x, y, z);
        double fint = cacFint(Reint, e, Dint);
        return fint;
    }

    double cacEs(int x, int y, int z) {
        return 1.0E-3;
    }

    double cacFint(double Reo, double e, double Do) {
        double temp;
        temp = pow(-3.6 * Math.log10(6.9 / Reo + pow(e / 3.7 / Do, 1.1)), -2.0);
        return temp;
    }

    public komegaTFOWAngeli2Wave1(double Qso, double Qsw, double Radia) {
        this.Qso = Qso;
        this.Qsw = Qsw;
        this.R = Radia;
    }

    /**
     * 计算湍流粘度
     *
     * @param rho
     * @param Mut
     * @param K
     * @param omega
     */
    public void cacMut(double[][][] rho, double[][][] Mut, Field K, Field omega) {
        double alpha = 1.0;
        for (int z = 1; z < K.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < K.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < K.getNewField().length - 1; ++x) {
                    Mut[x][y][z]
                            = alpha * rho[x][y][z] * K.getNewField()[x][y][z]
                            / (omega.getNewField()[x][y][z] + 1e-30);
                }
            }
        }
    }

    /**
     * 设定湍流动能
     *
     * @param K
     * @param W
     */
    public void setK(Field K, Field W) {
        for (int z = 1; z < K.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < K.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < K.getNewField().length - 1; ++x) {
                    double value = 1e-4 * pow(W.getNewField()[x][y][z], 2.0);
                    K.setNewField(x, y, z, value);
                }
            }
        }
    }

    /**
     * 设定能量耗散率
     *
     * @param omega
     * @param K
     * @param bm
     */
    public void setOmega(Field omega, Field K, BiSingleV2 bm) {
        for (int z = 1; z < omega.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < omega.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < omega.getNewField().length - 1; ++x) {
                    double value
                            = pow(0.09, 0.75)
                            * pow(K.getNewField()[x][y][z], 3.0 / 2.0)
                            / Turbulence.turbulenceLengthscale(bm.radia * 2);
                    value = value / K.getNewField()[x][y][z];
                    omega.setNewField(x, y, z, value);
                }
            }
        }
    }

    /**
     * @param Qw Qw
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
     * @param Qoil Qoil
     * @return
     */
    public double G(double Qoil) {
        double Qcoil = 0;
        Qcoil = Tool.getTotalVol(FO[W], OilMesh);
        System.out.println("Qcoil = " + Qcoil);
        return Qcoil - Qoil;
    }

    void RESULT(int count) throws FileNotFoundException {
        String dirName = position + this.toString() + "/";
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
                RES.printf("%16.6E\t", porW[S].getF()[I][J][k]);
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
                RES.printf("%16.6E\t", porO[S].getF()[I][J][k]);
                RES.println();
            }
        }
        RES.close();
    }

    @Override
    public String toString() {
        return getClass().getName() + Title;
    }
}
