package cn.org.higis.radar;

import cn.org.higis.radar.GeoSOT3D;
import cn.org.higis.radar.Position;

public class Demo {
    public static void main(String[] args) {
        Position pos = new Position();
        pos.lon[0] = 145;
        pos.lon[1] = 12;
        pos.lon[2] = 5;
        pos.lat[0] = 75;
        pos.lat[1] = 50;
        pos.lat[2] = 10;
        pos.alt = 10000;
        System.out.println("original pos");
        System.out.printf("%d:%d:%f\n", (int) pos.lon[0], (int) pos.lon[1], pos.lon[2]);
        System.out.printf("%d:%d:%f\n", (int) pos.lat[0], (int) pos.lat[1], pos.lat[2]);
        System.out.println(pos.alt);
        long code = GeoSOT3D.encode(pos, 9);
        Position pos_inv = GeoSOT3D.decode(code, 9);
        System.out.println("encode to code:");
        System.out.println(code);
        System.out.println("decode to grid center:");
        System.out.printf("%d:%d:%f\n", (int) pos_inv.lon[0], (int) pos_inv.lon[1], pos_inv.lon[2]);
        System.out.printf("%d:%d:%f\n", (int) pos_inv.lat[0], (int) pos_inv.lat[1], pos_inv.lat[2]);
        System.out.println(pos_inv.alt);
    }
}
