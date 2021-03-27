package com.winswe.io;

/**
 * writer和input需要进一步改进
 *
 * @author 齐雪宇
 */
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;

public class Writer {

    /**
     *
     * @param destFileName
     * @return
     */
    public static boolean createFile(String destFileName) {
        File file = new File(destFileName);
        if (file.exists()) {
            System.out.println("Create a file" + destFileName + "Failure，the file has exist！");
            return false;
        }
        if (destFileName.endsWith(File.separator)) {
            System.out.println("Create a file" + destFileName + "Failure，The target file can not be directory!");
            return false;
        }
        //判断目标文件所在的目录是否存在  
        if (!file.getParentFile().exists()) {
            //如果目标文件所在的目录不存在，则创建父目录  
            System.out.println("The destination file directory does not exist, prepare to create it!");
            if (!file.getParentFile().mkdirs()) {
                System.out.println("Fail to create destinatory directory！");
                return false;
            }
        }
        //创建目标文件  
        try {
            if (file.createNewFile()) {
                System.out.println("Create a file" + destFileName + "Successful!");
                return true;
            } else {
                System.out.println("Create a file" + destFileName + "Failure!");
                return false;
            }
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("Create a file" + destFileName + "Failure！" + e.getMessage());
            return false;
        }
    }

    /**
     *
     * @param destDirName
     * @return
     */
    public static boolean createDir(String destDirName) {
        File dir = new File(destDirName);
        if (dir.exists()) {
            System.out.println("Create direction" + destDirName + "Failure，the direction has exist!");
            return false;
        }
        if (!destDirName.endsWith(File.separator)) {
            destDirName = destDirName + File.separator;
        }
        //创建目录  
        if (dir.mkdirs()) {
            System.out.println("Create direction" + destDirName + "Successful！");
            return true;
        } else {
            System.out.println("Create direction" + destDirName + "Failure！");
            return false;
        }
    }

    /**
     *
     * @param fileName
     * @param titlename
     * @param phi
     * @param x
     * @param y
     * @throws FileNotFoundException
     */
    static private void inputField(String fileName, String titlename, double[] phi,
            double dx, int count) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(fileName)));
        //文件抬头
        pw.println("Title=" + "\"" + titlename + "\"");
        String var = "Variables=\"X\"";
        pw.println(var + "\"" + titlename + "\"");
        pw.println("Zone" + " I=" + phi.length + " F=POINT");
        for (int i = 0; i < phi.length; i++) {
            pw.printf("%16.6f\t", dx * i);
            pw.printf("%16.6f\t", phi[i]);
            pw.println();
        }
        pw.close();
    }

    static private void inputField(String name, String titlename, double[][] phi, double[] x, double[] y) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(name)));
        //文件抬头
        pw.println("Title=" + "\"" + titlename + "\"");
        String var = "Variables=\"X\",\"Y\"";
        pw.println(var + "\"" + titlename + "\"");
        pw.println("Zone" + " I=" + phi.length + " J=" + phi[0].length + " F=POINT");
        for (int j = 0; j < phi[0].length; j++) {
            for (int i = 0; i < phi.length; i++) {
                pw.printf("%16.6f\t", x[i]);
                pw.printf("%16.6f\t", y[j]);
                pw.printf("%16.6f\t", phi[i][j]);
                pw.println();
            }
        }
        pw.close();
    }

    /**
     *
     * @param name 文件名字
     * @param titlename 标题
     * @param phi 场
     * @param U U场
     * @param V V场
     * @param x x坐标
     * @param y y坐标
     * @throws FileNotFoundException
     */
    static private void inputField(String name, String titlename, double[][] phi, double[][] U, double[][] V, double[] x, double[] y) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(name)));
        //文件抬头
        pw.println("Title=" + "\"" + titlename + "\"");
        String var = "Variables=\"X\",\"Y\"";
        pw.println(var
                + "\"" + titlename + "\""
                + "\"" + "U" + "\""
                + "\"" + "V" + "\"");
        pw.println("Zone" + " I=" + phi.length + " J=" + phi[0].length + " F=POINT");
        for (int j = 0; j < phi[0].length; j++) {
            for (int i = 0; i < phi.length; i++) {
                pw.printf("%16.6f\t", x[i]);
                pw.printf("%16.6f\t", y[j]);
                pw.printf("%16.6f\t", phi[i][j]);
                pw.printf("%16.6f\t", U[i][j]);
                pw.printf("%16.6f\t", V[i][j]);
                pw.println();
            }
        }
        pw.close();
    }

    /**
     *
     * @param dirName 路径名
     * @param formate 格式“dat”
     * @param titlename 文件名
     * @param phi 场
     * @param dx
     * @param count
     * @throws FileNotFoundException
     */
    static public void writerArrayField(String dirName, String formate, String titlename, double[] phi, double dx, int count) throws FileNotFoundException {
        //创建目录  
        Writer.createDir(dirName);
        //创建文件  
        String str = String.format("%s", titlename);
        String fileName = dirName + str + "." + formate;
        Writer.createFile(fileName);
        Writer.inputField(fileName, titlename, phi, dx, count);
    }

    /**
     *
     * @param dirName 路径名
     * @param num 第几号文件
     * @param formate 格式“dat”
     * @param titlename 文件名
     * @param phi 场
     * @param x x坐标
     * @param y y坐标
     * @throws FileNotFoundException
     */
    static public void writerArrayField(String dirName, double num, String formate, String titlename, double[][] phi, double[] x, double[] y) throws FileNotFoundException {
        //创建目录  
        Writer.createDir(dirName);
        //创建文件  
        String str = String.format("%s_%f", titlename, num);
        String fileName = dirName + str + "." + formate;
        Writer.createFile(fileName);
        Writer.inputField(fileName, titlename, phi, x, y);
    }

    /**
     *
     * @param dirName 路径名
     * @param num 第几号文件
     * @param formate 格式“dat”
     * @param titlename 文件名
     * @param phi 场
     * @param U U速度
     * @param V V速度
     * @param x x坐标
     * @param y y坐标
     * @throws FileNotFoundException
     */
    static public void writerArrayField(String dirName, double num, String formate, String titlename,
            double[][] phi, double[][] U, double[][] V, double[] x, double[] y) throws FileNotFoundException {
        //创建目录  
        Writer.createDir(dirName);
        //创建文件  
        String str = String.format("%s_%f", titlename, num);
        String fileName = dirName + str + "." + formate;
        Writer.createFile(fileName);
        Writer.inputField(fileName, titlename, phi, U, V, x, y);
    }

    static public void writerMesh(String dirName, double num, String formate, String titlename,
            double[][][] phi, double[] x, double[] y, double[] z) throws FileNotFoundException {
        //创建目录  
        Writer.createDir(dirName);
        //创建文件  
        String str = String.format("%s_%f", titlename, num);
        String fileName = dirName + str + "." + formate;
        Writer.createFile(fileName);
        Writer.inputField(fileName, titlename, phi, x, y, z);
    }

    private static void inputField(String fileName, String titlename,
            double[][][] phi, double[] x, double[] y, double[] z) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(fileName)));
        //文件抬头
        pw.println("Title=" + "\"" + titlename + "\"");
        String var = "Variables=\"X\",\"Y\",\"Z\"";
        pw.println(var + "\"" + titlename + "\"");
        pw.println("Zone" + " I=" + phi.length + " J=" + phi[0].length + " K=" + phi[0][0].length + " F=POINT");
        for (int k = 0; k < z.length; ++k) {
            for (int j = 0; j < y.length; ++j) {
                for (int i = 0; i < x.length; ++i) {
                    pw.printf("%16.6f\t", x[i]);
                    pw.printf("%16.6f\t", y[j]);
                    pw.printf("%16.6f\t", z[k]);
                    pw.printf("%16.6f\t", phi[i][j][k]);
                    pw.println();
                }
            }
        }
        pw.close();
    }

    static public void writerArrayField(String dirName, double num, String formate, String titlename,
            double[][][] phi, double[][][] U, double[][][] V, double[] x, double[] y, double[] z) throws FileNotFoundException {
        //创建目录  
        Writer.createDir(dirName);
        //创建文件  
        String str = String.format("%s_%f", titlename, num);
        String fileName = dirName + str + "." + formate;
        Writer.createFile(fileName);
        Writer.inputField(fileName, titlename, phi, U, V, x, y, z);
    }

    private static void inputField(String fileName, String titlename,
            double[][][] phi, double[][][] U, double[][][] V, double[] x, double[] y, double[] z) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(fileName)));
        //文件抬头
        pw.println("Title=" + "\"" + titlename + "\"");
        String var = "Variables=\"X\",\"Y\"";
        pw.println(var
                + "\"" + titlename + "\""
                + "\"" + "U" + "\""
                + "\"" + "V" + "\"");
//        pw.println("Zone" + " I=" + phi.length + " J=" + phi[0].length + " K=" + phi[0][0].length + " F=POINT");
        pw.println("Zone" + " I=" + phi.length + " J=" + phi[0].length + " F=POINT");
        if (z.length == 3) {
            for (int k = 1; k < z.length - 1; ++k) {
                for (int j = 0; j < y.length; ++j) {
                    for (int i = 0; i < x.length; ++i) {
                        pw.printf("%16.6f\t", x[i]);
                        pw.printf("%16.6f\t", y[j]);
//                        pw.printf("%16.6f\t", z[k]);
                        pw.printf("%16.6f\t", phi[i][j][k]);
                        pw.printf("%16.6f\t", U[i][j][k]);
                        pw.printf("%16.6f\t", V[i][j][k]);
                        pw.println();
                    }
                }
            }
        }
        pw.close();
    }

    /**
     * **************Debug program******************
     */
    /**
     * 2D 数据
     *
     * @param name 第几号文件
     * @param titlename 文件名
     * @param phi 数据
     * @throws FileNotFoundException
     */
    static private void inputFieldDebug(String name, String titlename, double[][] phi) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(name)));
        //文件抬头
        pw.println("Title=" + "\"" + titlename + "\"");
        pw.println("Zone" + " I=" + phi.length + " J=" + phi[0].length + " F=POINT");
        for (int j = phi[0].length - 1; j >= 0; j--) {
            for (int i = 0; i < phi.length; i++) {
                pw.printf("(%d,%d)%.4f\t", i, j, phi[i][j]);
            }
            pw.println();
        }
        pw.close();
    }

    /**
     * 3D 数据
     *
     * @param name 第几号文件
     * @param titlename 文件名
     * @param phi 数据
     * @throws FileNotFoundException
     */
    static private void inputFieldDebug(String name, String titlename, double[][][] phi) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(name)));
        //文件抬头
        pw.println("Title=" + "\"" + titlename + "\"");
        pw.println("Zone"
                + " I=" + phi.length
                + " J=" + phi[0].length
                + " K=" + phi[0][0].length + " F=POINT");
        for (int k = 0; k < phi[0][0].length; k++) {
            pw.println(titlename + "[" + k + "]");
            for (int j = phi[0].length - 1; j >= 0; j--) {
                for (int i = 0; i < phi.length; i++) {
                    pw.printf("(%d,%d)%.4f\t", i, j, phi[i][j][k]);
                }
                pw.println();
            }
        }
        pw.close();
    }

    /**
     * 3D 数据
     *
     * @param name 第几号文件
     * @param titlename 文件名
     * @param phi 数据
     * @throws FileNotFoundException
     */
    static private void inputFieldDebug(String name, String titlename, int[][][] phi) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(name)));
        //文件抬头
        pw.println("Title=" + "\"" + titlename + "\"");
        pw.println("Zone"
                + " I=" + phi.length
                + " J=" + phi[0].length
                + " K=" + phi[0][0].length + " F=POINT");
        for (int k = 0; k < phi[0][0].length; k++) {
            pw.println(titlename + "[" + k + "]");
            for (int j = phi[0].length - 1; j >= 0; j--) {
                for (int i = 0; i < phi.length; i++) {
                    pw.printf("(%d,%d)%d\t", i, j, phi[i][j][k]);
                }
                pw.println();
            }
        }
        pw.close();
    }

    /**
     * 调试使用2D数据
     *
     * @param dirName 路径名
     * @param num 第几号文件
     * @param formate 格式“dat”
     * @param titlename 文件名
     * @param phi 数据
     * @throws FileNotFoundException
     */
    static public void writer2DArrayDebug(String dirName, int num, String formate, String titlename, double[][] phi) throws FileNotFoundException {
        //创建目录  
        Writer.createDir(dirName);
        //创建文件  
        String str = String.format("%d", num);
        String fileName = dirName + str + "." + formate;
        Writer.createFile(fileName);
        Writer.inputFieldDebug(fileName, titlename, phi);
    }

    /**
     * 调试使用3D数据
     *
     * @param dirName 路径名
     * @param num 第几号文件
     * @param formate 格式“dat”
     * @param titlename 文件名
     * @param phi 数据
     * @throws FileNotFoundException
     */
    static public void writer3DArrayDebug(String dirName, int num, String formate, String titlename, double[][][] phi) throws FileNotFoundException {
        //创建目录  
        Writer.createDir(dirName);
        //创建文件  
        String str = String.format("%d", num);
        String fileName = dirName + str + "." + formate;
        Writer.createFile(fileName);
        Writer.inputFieldDebug(fileName, titlename, phi);
    }

    /**
     * 调试使用3D数据
     *
     * @param dirName 路径名
     * @param num 第几号文件
     * @param formate 格式“dat”
     * @param titlename 文件名
     * @param phi 数据
     * @throws FileNotFoundException
     */
    static public void writer3DArrayDebug(String dirName, int num, String formate, String titlename, int[][][] phi) throws FileNotFoundException {
        //创建目录  
        Writer.createDir(dirName);
        //创建文件  
        String str = String.format("%d", num);
        String fileName = dirName + str + "." + formate;
        Writer.createFile(fileName);
        Writer.inputFieldDebug(fileName, titlename, phi);
    }

}
