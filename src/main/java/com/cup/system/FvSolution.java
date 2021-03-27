package com.cup.system;

/**
 * <n>这里需要改进，在选择求解器，设置精度等方面需要新建类。</n>
 * <n>最好定义一个容器进行存放</n>
 *
 * @author winsway
 */
public class FvSolution {

    /**
     * 时间精度
     */
    double ttol;
    /**
     * 初始残差
     */
    double initR;
    /**
     * 绝对残差
     */
    double atol;
    /**
     * 相对残差
     */
    double rtol;
    /**
     * 最大迭代次数
     */
    int maxIter;

    /**
     * 构造默认的求解精度
     */
    public FvSolution() {
        this.maxIter = 100000;
        this.rtol = 1e-6;
        this.atol = 1e-8;
        this.ttol = 1e-3;
    }

    /**
     * 设定最大迭代次数
     *
     * @param maxIter 最大迭代次数
     */
    public void setMaxIterations(int maxIter) {
        this.maxIter = maxIter;
    }

    /**
     * 设定相对收敛容差
     *
     * @param rtol Relative convergence tolerance (to initial residual)
     */
    public void setRelativeTolerance(double rtol) {
        this.rtol = rtol;
    }

    /**
     * 设定绝对容差
     *
     * @param atol Absolute convergence tolerance
     */
    public void setAbsoluteTolerance(double atol) {
        this.atol = atol;
    }

    /**
     * 设定时间收敛标准
     *
     * @param ttol 时间容差
     */
    public void setTimeTolerance(double ttol) {
        this.ttol = ttol;
    }

    /**
     *
     * @return 时间精度
     */
    public double getTtol() {
        return ttol;
    }

    /**
     *
     * @return 初始残差
     */
    public double getInitR() {
        return initR;
    }

    /**
     *
     * @return 返回绝对残差
     */
    public double getAtol() {
        return atol;
    }

    /**
     *
     * @return 相对残差
     */
    public double getRtol() {
        return rtol;
    }

    /**
     *
     * @return 最大迭代次数
     */
    public int getMaxIter() {
        return maxIter;
    }

    @Override
    public FvSolution clone() throws CloneNotSupportedException {
        FvSolution clone = new FvSolution();
        clone.atol = this.atol;
        clone.initR = this.initR;
        clone.maxIter = this.maxIter;
        clone.rtol = this.rtol;
        clone.ttol = this.ttol;
        return clone;
    }

}
