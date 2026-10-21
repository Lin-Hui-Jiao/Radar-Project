package cn.org.higis.radar;

import cn.org.higis.radar.Position;

public class GeoSOT3D {

    static {
        System.loadLibrary("geosot3djni");
    }

    public static native long encode(Position center, int level);

    public static Position decode(long code, int level) {
        Position pos = new Position();
        decode(pos, code, level);
        return pos;
    }

    private static native void decode(Position pos, long code, int level);
}