package com.cup.system;

import com.cup.stratifiedflow.*;
import com.cup.system.*;

/**
 * <li>1.输出控制</li>
 * <li>2.需要加入场的读入</li>
 * <li>3.控制瞬态和稳态</li>
 *
 * @author winsway
 */
public class ControlDict {

    public double endTime = 20;
    public double deltaT = 0.01;
    public int writeInterval = 100;
    public String position;
    public int count;
    public int max_it;

    public ControlDict() {

    }

    public ControlDict(double endTime, int writeInterval) {
        this.endTime = endTime;
        this.writeInterval = writeInterval;
    }

    public boolean writeDate() {
        if (count % writeInterval == 0 || count == 0) {
            return true;
        } else {
            return false;
        }
    }

    public double getTime() {
        return count * deltaT;
    }

    public boolean judgeStop() {
        if (this.getTime() < this.endTime) {
            return true;
        } else {
            return false;
        }
    }

    public int getMax_it() {
        max_it = (int) (endTime / deltaT);
        return max_it;
    }

    @Override
    public ControlDict clone() throws CloneNotSupportedException {
        ControlDict clone = new ControlDict();
        clone.count = this.count;
        clone.deltaT = this.deltaT;
        clone.endTime = this.endTime;
        clone.position = this.position;
        clone.writeInterval = this.writeInterval;
        return clone;
    }

}
