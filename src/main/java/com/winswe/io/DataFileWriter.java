package com.winswe.io;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

public class DataFileWriter implements AutoCloseable {

    public final PrintWriter fileWriter;
    public final String commentStr;

    public DataFileWriter(File file) throws FileNotFoundException, IOException {
        this(file, "%");
    }

    /**
     *
     * @param file 读取的数据文件
     * @param commentStr 注释标记符
     * @throws FileNotFoundException
     */
    public DataFileWriter(File file, String commentStr) throws FileNotFoundException, IOException {
        if (!file.getParentFile().exists()) {
            file.getParentFile().mkdirs();
        }
        if (!file.exists()) {
            file.createNewFile();
        }
        fileWriter = new PrintWriter(file);
        this.commentStr = commentStr;
    }

    @Override
    public void close() {
        fileWriter.close();
        System.out.println("文件创建完成！");
    }

    public void writerDouble(double var) {
        fileWriter.format("%e\t", var);
    }

    public void writerIndex(int var) {
        fileWriter.format("%6d", var);
    }

    public void writerDataLine(double[] var) {
        for (int i = 0; i < var.length; i++) {
            writerDouble(var[i]);
        }
        fileWriter.println();
    }

    public void println() {
        fileWriter.println();
    }

}
