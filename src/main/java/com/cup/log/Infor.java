package com.cup.log;

/**
 * <li>日志的输出控制</li>
 * <li>日志需要输出到一个文本文件中</li>
 * @author winsway
 */
public class Infor {

    public double CourantNumbermean, CourantNumbermax;
    public ResidualLog rl = new ResidualLog();
    public double timeResidual;

    public void timeR(double fieldN, double fieldO) {
        if (Math.abs(fieldN - fieldO) > this.timeResidual) {
            this.timeResidual = Math.abs(fieldN - fieldO);
        }
    }

    public void log(String var) {
        System.out.printf("%s   timeResidual=%e  Final_residual=%e  iteration=%d\n",
                var,
                this.timeResidual,
                rl.Final_residual,
                rl.No_Iteration
        );
        this.setZero();
    }

    public void setZero() {
        this.timeResidual = 0;
        rl.setZero();
    }
}
