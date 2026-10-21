package cn.org.higis.radar;

public class Radar {

    static {
        System.loadLibrary("radarjni");
    }

    Radar(String name, int N, float d, float f, float G, float Pt, float mintheta, float maxtheta, float minphi,
            float maxphi, float maxRange) {
        init(name, N, d, f, G, Pt, mintheta, maxtheta, minphi, maxphi, maxRange);
    }

    // Declare a native method sayHello() that receives no arguments and returns
    // void
    private native void init(String name, int N, float d, float f, float G, float Pt, float mintheta, float maxtheta,
            float minphi, float maxphi, float maxRange);

    public native void move(Position pos, float basePhi, float baseTheta, float timestamp);

    public native void setGridParams(float range, int level);

    public native void capablePowerDensityGrid();

    public native void capableRegion(float rcs);

    public native void dispose();

    static native void begin_mpi();

    static native void end_mpi();

    private long radarImpl = 0; // c++ radar object
    private int level = 10;
    private float range = 320000;
    public String dbHost = "127.0.0.1";
    public int dbPort = 9000;
    public int id = 0;
}