#include <iostream>
#include "cn_org_higis_radar_GeoSOT3D.h"
#include "geosot3d.hpp"

struct Position
{
    float lon[3] = {0};
    float lat[3] = {0};
    float alt = 0;
};

jfieldID GetFid(JNIEnv *env, jobject obj, const char *name, const char *type)
{
    auto cls = env->GetObjectClass(obj);
    return env->GetFieldID(cls, name, type);
}

Position JPositionToPosition(JNIEnv *env, jobject pos)
{
    auto lon_obj = env->GetObjectField(pos, GetFid(env, pos, "lon", "[F"));
    auto lat_obj = env->GetObjectField(pos, GetFid(env, pos, "lat", "[F"));
    auto alt = env->GetFloatField(pos, GetFid(env, pos, "alt", "F"));
    auto lon = reinterpret_cast<jfloatArray *>(&lon_obj);
    auto lat = reinterpret_cast<jfloatArray *>(&lat_obj);
    auto lon_dms = env->GetFloatArrayElements(*lon, 0);
    auto lat_dms = env->GetFloatArrayElements(*lat, 0);
    Position position = {{lon_dms[0], lon_dms[1], lon_dms[2]},
                         {lat_dms[0], lat_dms[1], lat_dms[2]},
                         alt};
    env->ReleaseFloatArrayElements(*lon, lon_dms, 0);
    env->ReleaseFloatArrayElements(*lat, lat_dms, 0);
    return position;
}

JNIEXPORT jlong JNICALL Java_cn_org_higis_radar_GeoSOT3D_encode(JNIEnv *env, jclass, jobject jpos, jint level)
{
    auto pos = JPositionToPosition(env, jpos);
    auto res = GeoSOT3D::Encode(pos.lon, pos.lat, pos.alt, level);
    std::cout << pos.lon[0] << ":" << pos.lon[1] << ":" << pos.lon[2] << ",";
    std::cout << pos.lat[0] << ":" << pos.lat[1] << ":" << pos.lat[2] << "," << pos.alt << std::endl;
    std::cout << res << std::endl;
    return res;
}

JNIEXPORT void JNICALL Java_cn_org_higis_radar_GeoSOT3D_decode(JNIEnv *env, jclass, jobject jpos, jlong code, jint level)
{
    auto lon_obj = env->GetObjectField(jpos, GetFid(env, jpos, "lon", "[F"));
    auto lat_obj = env->GetObjectField(jpos, GetFid(env, jpos, "lat", "[F"));
    auto lon = reinterpret_cast<jfloatArray *>(&lon_obj);
    auto lat = reinterpret_cast<jfloatArray *>(&lat_obj);
    auto lon_dms = env->GetFloatArrayElements(*lon, 0);
    auto lat_dms = env->GetFloatArrayElements(*lat, 0);
    float alt;
    GeoSOT3D::Decode(lon_dms, lat_dms, alt, code, level);
    env->SetFloatField(jpos, GetFid(env, jpos, "alt", "F"), alt);
    env->ReleaseFloatArrayElements(*lon, lon_dms, 0);
    env->ReleaseFloatArrayElements(*lat, lat_dms, 0);
}
