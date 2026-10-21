package cn.org.higis.radar;

public class Radar {

    static {
        System.loadLibrary("radarjni");
    }

    Radar(String name) {
        switch(name)
        {
            case "AN/APG-66":
                initAnApg66();
                break;
            case "AN/APS-145":
                initAnAps145();
                break;
            case "AN/APS-115":
                initAnAps115();
                break;
            case "AN/APY-1":
                initAnApy1();
                break;
            case "AN/APY-2":
                initAnApy2();
                break;
            case "AN/SPS-49":
                initAnSps49();
                break;
            case "SPY-1D(V)":
                initSpy1d();
                break;
            case "AN/MPQ-65":
                initAnMpq65();
                break;
            default:
                System.out.println("no such radar");
        }
    }

    private native void initAnApg66();
    private native void initAnAps145();
    private native void initAnAps115();
    private native void initAnApy1();
    private native void initAnApy2();
    private native void initAnSps49();
    private native void initSpy1d();
    private native void initAnMpq65();

    public native void move(Position pos, float basePhi, float baseTheta, float timestamp);

    public native void setGridParams(float range, int level);

    public native void capablePowerDensityGrid();

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