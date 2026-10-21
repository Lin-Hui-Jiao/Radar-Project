#include <iostream>
#include <string>
#include <mpi.h>
#include "cn_org_higis_radar_Radar.h"
#include "hiradar.hpp"
#include "hiradar/radar_pool.hpp"

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

    return position;
}

jfieldID GetDbHostFid(JNIEnv *env, jobject obj)
{
    auto cls = env->GetObjectClass(obj);
    auto host_fid = env->GetFieldID(cls, "dbHost", "Ljava/lang/String;");
    return host_fid;
}

jfieldID GetPtrFieldId(JNIEnv *env, jobject obj)
{
    static jfieldID ptrFieldId = 0;

    if (!ptrFieldId)
    {
        jclass c = env->GetObjectClass(obj);
        ptrFieldId = env->GetFieldID(c, "radarImpl", "J");
        env->DeleteLocalRef(c);
    }

    return ptrFieldId;
}

Radar *GetRadarImpl(JNIEnv *env, jobject obj)
{
    return (Radar *)env->GetLongField(obj, GetPtrFieldId(env, obj));
}

void SetRadarImpl(JNIEnv *env, jobject obj, IRadar *radar_impl)
{
    env->SetLongField(obj, GetPtrFieldId(env, obj), (jlong)radar_impl);
    auto host_fid = GetDbHostFid(env, obj);
    auto db_host = (jstring)env->GetObjectField(obj, host_fid);
    auto db_host_cstr = env->GetStringUTFChars(db_host, 0);
    auto db_port = env->GetIntField(obj, GetFid(env, obj, "dbPort", "I"));
    auto level = env->GetIntField(obj, GetFid(env, obj, "level", "I"));
    auto id = env->GetIntField(obj, GetFid(env, obj, "id", "I"));
    radar_impl->BindWriter(new GeoSotChPowerDensityWriter(std::string(db_host_cstr), db_port, "radar", id, level));
    env->ReleaseStringUTFChars(db_host, db_host_cstr);
}

JNIEXPORT void JNICALL Java_cn_org_higis_radar_Radar_initAnApg66(JNIEnv *env, jobject obj)
{
    Position pos;
    SetRadarImpl(env, obj, ANAPG66(pos, 0, 0));
}

JNIEXPORT void JNICALL Java_cn_org_higis_radar_Radar_initAnAps145(JNIEnv *env, jobject obj)
{
    Position pos;
    SetRadarImpl(env, obj, ANAPS145(pos, 0, 0));
}

JNIEXPORT void JNICALL Java_cn_org_higis_radar_Radar_initAnAps115(JNIEnv *env, jobject obj)
{
    Position pos;
    SetRadarImpl(env, obj, ANAPS115(pos, 0, 0));
}

JNIEXPORT void JNICALL Java_cn_org_higis_radar_Radar_initAnApy1(JNIEnv *env, jobject obj)
{
    Position pos;
    SetRadarImpl(env, obj, ANAPY12(pos, 0, 0));
}

JNIEXPORT void JNICALL Java_cn_org_higis_radar_Radar_initAnApy2(JNIEnv *env, jobject obj)
{
    Position pos;
    SetRadarImpl(env, obj, ANAPY12(pos, 0, 0));
}

JNIEXPORT void JNICALL Java_cn_org_higis_radar_Radar_initAnSps49(JNIEnv *env, jobject obj)
{
    Position pos;
    SetRadarImpl(env, obj, ANSPS49(pos, 0, 0));
}

JNIEXPORT void JNICALL Java_cn_org_higis_radar_Radar_initSpy1d(JNIEnv *env, jobject obj)
{
    Position pos;
    SetRadarImpl(env, obj, SPY1D(pos, 0, 0));
}

JNIEXPORT void JNICALL Java_cn_org_higis_radar_Radar_initAnMpq65(JNIEnv *env, jobject obj)
{
    Position pos;
    SetRadarImpl(env, obj, ANMPQ65(pos, 0, 0));
}

JNIEXPORT void JNICALL Java_cn_org_higis_radar_Radar_move(JNIEnv *env, jobject obj, jobject pos, jfloat base_phi, jfloat base_theta, jfloat timestamp)
{
    Radar *radar_impl = GetRadarImpl(env, obj);
    auto position = JPositionToPosition(env, pos);
    radar_impl->Move(position, base_phi, base_theta, timestamp);
}

JNIEXPORT void JNICALL Java_cn_org_higis_radar_Radar_setGridParams(JNIEnv *env, jobject obj, jfloat range, jint level)
{
    Radar *radar_impl = GetRadarImpl(env, obj);
    env->SetIntField(obj, GetFid(env, obj, "level", "I"), level);
    env->SetFloatField(obj, GetFid(env, obj, "range", "F"), range);
    ((GeoSotChPowerDensityWriter *)radar_impl->GetWriter())->SetLevel(level);
}

JNIEXPORT void JNICALL Java_cn_org_higis_radar_Radar_capablePowerDensityGrid(JNIEnv *env, jobject obj)
{
    Radar *radar_impl = GetRadarImpl(env, obj);
    auto level = env->GetIntField(obj, GetFid(env, obj, "level", "I"));
    auto range = env->GetFloatField(obj, GetFid(env, obj, "range", "F"));
    auto pos_list = CreateGrid(radar_impl->GetPosition(), range, level);
    radar_impl->CapablePowerDensity(pos_list.first, pos_list.second);
    free(pos_list.first);
}

JNIEXPORT void JNICALL Java_cn_org_higis_radar_Radar_dispose(JNIEnv *env, jobject obj)
{
    auto radar_impl = GetRadarImpl(env, obj);
    if (radar_impl->GetWriter() != NULL)
    {
        delete radar_impl->GetWriter();
    }
    delete radar_impl;
}

JNIEXPORT void JNICALL Java_cn_org_higis_radar_Radar_begin_1mpi(JNIEnv *, jclass)
{
    MPI_Init(NULL, NULL);
}

JNIEXPORT void JNICALL Java_cn_org_higis_radar_Radar_end_1mpi(JNIEnv *, jclass)
{
    MPI_Finalize();
}