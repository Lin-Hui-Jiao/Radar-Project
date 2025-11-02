#pragma once

#include <inttypes.h>

namespace GeoSOT3D
{
	uint64_t Encode(float lon[3], float lat[3], float alt_val, short level);
	void Decode(float out_lon[3], float out_lat[3], float &out_alt_val, uint64_t code, short level);
	uint64_t Encode(float lon, float lat, float alt, short level);
	void Decode(float &out_lon, float &out_lat, float &out_alt, uint64_t code, short level);
	// __uint128_t Encode_96(const float lon[3], const float lat[3], float alt_val, short level);
	// void Decode_96(float out_lon[3], float out_lat[3], float &out_alt_val, __uint128_t code, short level);
	// __uint128_t Encode_96(double lon_deg, double lat_deg, double alt_m, short level);
	// void Decode_96(double& out_lon_deg, double& out_lat_deg, double& out_alt_m, __uint128_t code, short level);
};