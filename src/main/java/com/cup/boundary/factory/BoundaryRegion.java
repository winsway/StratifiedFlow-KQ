/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.boundary.factory;

/**
 * 根据行压缩稀疏矩阵定义边界
 *
 * @author winsway
 */
public class BoundaryRegion extends AbstractBoundary {

    /**
     * 边界条件设定:内场(0),第一类边界（1），第二类边界（2），第三类边界（3）
     */
    public int[][][] ifBound;

    /**
     * 边界设定,通用边界
     */
    public RobinBC[][][] parameter;

    public BoundaryRegion(int X, int Y, int Z) {
        super(X, Y, Z);
        ifBound = new int[X][Y][Z];
        parameter = new RobinBC[X][Y][Z];
    }

    @Override
    public String getBoundaryTypes() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setType(int x, int y, int z, int value) {
        ifBound[x][y][z] = value;
    }

    public void setType(int x, int y, int z, int value, RobinBC robinbc) {
        ifBound[x][y][z] = value;
        parameter[x][y][z] = new RobinBC(robinbc);
    }

    @Override
    public int judgeBoundary(int m) {
        return (m == 0) ? 0 : 1;
    }

    @Override
    public int judgeBoundary(int x, int y, int z) {
        return judgeBoundary(ifBound[x][y][z]);
    }

    @Override
    public int getType(int x, int y, int z) {
        return ifBound[x][y][z];
    }

    @Override
    public BoundaryRegion getBoundaryRegion() {
        return this;
    }

    @Override
    public BoundaryRegion clone() {
        BoundaryRegion clone = new BoundaryRegion(numX, numY, numZ);
        for (int z = 0; z < this.numZ; ++z) {
            for (int j = 0; j < this.numY; ++j) {
                for (int i = 0; i < this.numX; ++i) {
                    clone.ifBound[i][j][z] = this.ifBound[i][j][z];
                    clone.parameter[i][j][z] = new RobinBC(this.parameter[i][j][z]);
                }
            }
        }
        return clone; //To change body of generated methods, choose Tools | Templates.
    }

}
