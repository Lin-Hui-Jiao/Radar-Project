package cn.org.higis.radar;

public class Demo {
    public static void main(String[] args) {
        Radar.begin_mpi();
        var radar = new Radar("custom", 1, 0.48f, 10000.0f, 32.0f, 20.0f, -30f, 30f, -30f, 30f);
        var pos = new Position();
        radar.move(pos, 0, 0, 0);
        // radar.move(pos, 45, 45, 0);
        radar.setGridParams(50000, 16);
        radar.capablePowerDensityGrid();
        radar.dispose();
        Radar.end_mpi();
    }
}