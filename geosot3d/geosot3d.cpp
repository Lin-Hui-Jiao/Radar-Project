#include "geosot3d.hpp"
#include "libmorton/morton.h"
#include <cstring>
#include <iostream>
#include <cmath>

#include <cstdint>
#define EARTH_R 6371.393

// this lib relies on https://github.com/Forceflow/libmorton
namespace GeoSOT3D
{
	//保留一个整数的低20位，并清除所有高于20位的比特,二进制表示为 0 1111 1111 1111 1111 1111
	static const uint64_t ClearMask = (1 << 20) - 1;
	/**
	 * 这个表有22行，每一行代表一个剖分层级（从第1级到第22级）。
	 * 每一行有3列，分别代表该层级在度、分、秒三个单位上的步长。
	 * 例如，第1行表示第1级，步长为256度，0分，0秒。
	 * 第2行表示第2级，步长为256度，0分，0秒。
	 * 第3行表示第3级，步长为128度，0分，0秒。
	 * 以此类推。
	 * 
	 * 这个表的目的是为了快速计算出每个层级在度、分、秒三个单位上的步长。
	 */
	static const float LevelStepMap[][3] = {
		{256, 0, 0},
		{256, 0, 0},
		{128, 0, 0},
		{64, 0, 0},
		{32, 0, 0},
		{16, 0, 0},
		{8, 0, 0},
		{4, 0, 0},
		{2, 0, 0},
		{1, 0, 0},
		{0, 32, 0},
		{0, 16, 0},
		{0, 8, 0},
		{0, 4, 0},
		{0, 2, 0},
		{0, 1, 0},
		{0, 0, 32},
		{0, 0, 16},
		{0, 0, 8},
		{0, 0, 4},
		{0, 0, 2},
		{0, 0, 1}};

	/**
	 * 设置一个整数中特定比特位的值。它可以在不影响其他位的情况下，
	 * 精确地将指定位置的比特（bit）设置为0或1。
	 */
	uint32_t SetBit(uint32_t num, size_t n, int bit)
	{
		return (num & ~(1UL << n)) | (bit << n);
	}
	/**
	 * 将一个整数角度值添加到一个已有的度分秒（DMS）坐标上，并对结果进行规范化（Normalization）处理。
	 */
	void AddDegree(uint32_t degree, float dms[3])
	{
		dms[0] += degree;
		if (dms[0] > 0)
		{
			if (dms[2] < 0)
			{
				dms[1]--;
				dms[2] += 60;
			}
			if (dms[1] < 0)
			{
				dms[0]--;
				dms[1] += 60;
			}
		}
	}

	//将一个表示总秒数的浮点值，转换为度分秒 (DMS) 格式。
	void SecondToDMS(float dms[3], float sec)
	{
		int sign = sec > 0 ? 1 : -1;  //有可能秒数值是负的
		sec = abs(sec);
		auto raw_m = int(sec) / 60;
		dms[2] = sec - raw_m * 60;
		dms[1] = raw_m % 60;
		dms[0] = raw_m / 60;
		for (size_t i = 0; i < 3; i++)
			dms[i] *= sign;
	}

	//这是上一个函数的逆操作，将度分秒 (DMS) 格式转换为总秒数。
	// float DMSToSecond(float dms[3])
	// {
	// 	return dms[2] * 3600 + dms[1] * 60 + dms[0];
	// }

	//将十进制度数转换为度分秒 (DMS) 格式。例如输入 116.5，输出的dms数组将约等于 {116.0, 30.0, 0.0}。
	void DegreeToDMS(float dms[3], float deg)
	{
		int sign = deg > 0 ? 1 : -1;
		deg = abs(deg);
		dms[0] = (int)deg;
		dms[1] = (int)((deg - (int)deg) * 60);
		dms[2] = (deg - dms[0] - dms[1] / 60) * 3600;
		for (size_t i = 0; i < 3; i++)
			dms[i] *= sign;
	}
	//这是上一个函数的逆操作，将度分秒 (DMS) 格式转换回十进制度数。
	float DMSToDegree(float dms[3])
	{
		return dms[0] + dms[1] / 60.0 + dms[2] / 3600.0;
	}




// 	inline __uint128_t spread_bits_32_to_96(uint32_t x) {
//         __uint128_t val = x;
//         // 以下每一步都将位的间距拉大，直到间隔2位
//         val = (val | (val << 32)) & 0x00000000FFFFFFFF00000000FFFFFFFF; // a..p -> a..p0..0a..p
//         val = (val | (val << 16)) & 0x0000FFFF0000FFFF0000FFFF0000FFFF; // a..h0..0a..h -> a..h0..0a..h0..0a..h
//         val = (val | (val << 8))  & 0x00FF00FF00FF00FF00FF00FF00FF00FF;
//         val = (val | (val << 4))  & 0x0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F;
//         val = (val | (val << 2))  & 0x33333333333333333333333333333333;
//         val = (val | (val << 1))  & 0x55555555555555555555555555555555;
//         return val;
//     }

//     // 辅助函数：将一个散开的96位数值，压缩回32位整数
//     inline uint32_t compact_bits_96_to_32(__uint128_t code) {
//         code &= 0x55555555555555555555555555555555;
//         code = (code | (code >> 1))  & 0x33333333333333333333333333333333;
//         code = (code | (code >> 2))  & 0x0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F;
//         code = (code | (code >> 4))  & 0x00FF00FF00FF00FF00FF00FF00FF00FF;
//         code = (code | (code >> 8))  & 0x0000FFFF0000FFFF0000FFFF0000FFFF;
//         code = (code | (code >> 16)) & 0x00000000FFFFFFFF00000000FFFFFFFF;
//         code = (code | (code >> 32)) & 0x0000000000000000FFFFFFFFFFFFFFFF;
//         return static_cast<uint32_t>(code);
//     }
    
//     // 编码：将3个32位整数交错成一个96位编码
//     // 论文中顺序为 高、纬、经，对应 z, y, x
//     __uint128_t mortonEncode_96bit(uint32_t x, uint32_t y, uint32_t z) {
//         __uint128_t z_spread = spread_bits_32_to_96(z);
//         __uint128_t y_spread = spread_bits_32_to_96(y);
//         __uint128_t x_spread = spread_bits_32_to_96(x);
//         return (z_spread << 2) | (y_spread << 1) | x_spread;
//     }

//     // 解码：将一个96位编码分解为3个32位整数
//     void mortonDecode_96bit(__uint128_t code, uint32_t& x, uint32_t& y, uint32_t& z) {
//         x = compact_bits_96_to_32(code);
//         y = compact_bits_96_to_32(code >> 1);
//         z = compact_bits_96_to_32(code >> 2);
//     }



// 	/**
//  * @brief 【新版本】将一个三维坐标（DMS格式）编码为96位的GeoSOT-3D Morton码。
//  * @param lon 输入的经度，DMS格式数组 float[3]
//  * @param lat 输入的纬度，DMS格式数组 float[3]
//  * @param alt_val 输入的海拔高度，单位：米
//  * @param level 剖分层级 (1-32)
//  * @return __uint128_t 96位编码值
//  */
// __uint128_t Encode_96(const float lon[3], const float lat[3], float alt_val, short level)
// {
//     // 先将DMS格式转换为更易于处理的十进制度数
//     double lon_deg = DMSToDegree(const_cast<float*>(lon));
//     double lat_deg = DMSToDegree(const_cast<float*>(lat));

//     // 调用我们新的、基于十进制度数的核心编码函数
//     return Encode_96(lon_deg, lat_deg, alt_val, level);
// }

// /**
//  * @brief 【新版本】将一个96位GeoSOT-3D Morton码解码为三维坐标（DMS格式）。
//  * @param out_lon 输出的经度，DMS格式数组 float[3]
//  * @param out_lat 输出的纬度，DMS格式数组 float[3]
//  * @param out_alt_val 输出的海拔高度，单位：米
//  * @param code 96位编码值
//  * @param level 剖分层级 (1-32)
//  */
// void Decode_96(float out_lon[3], float out_lat[3], float &out_alt_val, __uint128_t code, short level)
// {
//     double lon_deg, lat_deg,alt_val;
    
//     // 调用我们新的、基于十进制度数的核心解码函数
//     Decode_96(lon_deg, lat_deg, alt_val, code, level);

//     // 将解码出的十进制度数转换回DMS格式
//     DegreeToDMS(out_lon, lon_deg);
//     DegreeToDMS(out_lat, lat_deg);
//     out_alt_val = alt_val;
// }


// 	/**
// 	 * 将一个三维坐标（经度、纬度、高程）转换为GeoSOT-3D系统内部使用的64位整数编码。
// 	 * 它首先将高程值转换为0-512“度”单位，然后将其转换为度分秒 (DMS) 格式。
// 	 * 接着，它将经度、纬度、高程的DMS值转换为整数编码，并进行位交错处理。
// 	 * 最后，它调用第三方库 libmorton，将处理好的三个维度的整数编码，通过位交错的方式，“编织”成一个最终的64位一维整数。这就是降维的过程。
// 	 */
// 	/**
//  * @brief 【新版本-重载】将一个三维坐标（十进制度数）编码为96位GeoSOT-3D Morton码。
//  * 这是现在主要的编码函数。
//  */
// __uint128_t Encode_96(double lon_deg, double lat_deg, double alt_m, short level)
// {
//     // 这个函数就是我们之前编写的 Encode_96，现在作为主函数
//     if (level < 1 || level > 32) {
//         throw std::out_of_range("Level must be between 1 and 32.");
//     }

//     double dims_deg[3];
//     dims_deg[0] = (alt_m + EARTH_R * 1000.0) * 180.0 / (M_PI * EARTH_R * 1000.0);
//     dims_deg[1] = lat_deg;
//     dims_deg[2] = lon_deg;

//     uint32_t codes_32bit[3] = {0};

//     for (int i = 0; i < 3; ++i) {
//         double current_deg = dims_deg[i];
//         if (i > 0) { current_deg += 256.0; } // 经纬度偏移

//         uint32_t d_code = static_cast<uint32_t>(floor(current_deg));
//         double minutes_total = (current_deg - d_code) * 60.0;
//         uint32_t m_code = static_cast<uint32_t>(floor(minutes_total));
//         double seconds_total = (minutes_total - m_code) * 60.0;
//         uint32_t s_code = static_cast<uint32_t>(floor(seconds_total));
//         double fractional_second = seconds_total - s_code;
//         uint32_t fs_code = static_cast<uint32_t>(fractional_second * (1 << 11));

//         codes_32bit[i] = ((d_code & 0x1FF) << 23) | ((m_code & 0x3F) << 17) | ((s_code & 0x3F) << 11) | (fs_code & 0x7FF);
//     }

//     uint32_t alt_bits = codes_32bit[0] >> (32 - level);
//     uint32_t lat_bits = codes_32bit[1] >> (32 - level);
//     uint32_t lon_bits = codes_32bit[2] >> (32 - level);
    
//     return mortonEncode_96bit(lon_bits, lat_bits, alt_bits);
// }

// 	/**
//  * @brief 【新版本-重载】将一个96位GeoSOT-3D Morton码解码为三维坐标（十进制度数）。
//  * 这是现在主要的解码函数。
//  */
// void Decode_96(double& out_lon_deg, double& out_lat_deg, double& out_alt_m, __uint128_t code, short level)
// {
//     // 这个函数就是我们之前编写的 Decode_96，现在作为主函数
//     if (level < 1 || level > 32) {
//         throw std::out_of_range("Level must be between 1 and 32.");
//     }

//     uint32_t alt_bits, lat_bits, lon_bits;
//     mortonDecode_96bit(code, lon_bits, lat_bits, alt_bits);

//     uint32_t codes_32bit[3];
//     codes_32bit[0] = alt_bits << (32 - level);
//     codes_32bit[1] = lat_bits << (32 - level);
//     codes_32bit[2] = lon_bits << (32 - level);

//     double dims_deg[3];

//     for (int i = 0; i < 3; ++i) {
//         uint32_t d_code = (codes_32bit[i] >> 23) & 0x1FF;
//         uint32_t m_code = (codes_32bit[i] >> 17) & 0x3F;
//         uint32_t s_code = (codes_32bit[i] >> 11) & 0x3F;
//         uint32_t fs_code = codes_32bit[i] & 0x7FF;

//         dims_deg[i] = static_cast<double>(d_code) +
//                       static_cast<double>(m_code) / 60.0 +
//                       static_cast<double>(s_code) / 3600.0 +
//                       static_cast<double>(fs_code) / (3600.0 * (1 << 11));
        
//         if (i > 0) { dims_deg[i] -= 256.0; } // 恢复偏移
//     }
    
//     out_alt_m = (dims_deg[0] * M_PI * EARTH_R * 1000.0) / 180.0 - (EARTH_R * 1000.0);
//     out_lat_deg = dims_deg[1];
//     out_lon_deg = dims_deg[2];
// }




	uint64_t Encode(float lon[3], float lat[3], float alt_val, short level)
	{
		auto alt_deg = ((alt_val + EARTH_R) / (50000 + EARTH_R)) * 512; //将以米为单位的绝对高程（alt_val）通过线性映射，转换为GeoSOT-3D系统内部使用的0-512“度”单位。
		float alt[3];
		DegreeToDMS(alt, alt_deg); //将十进制转换为度分秒DMS的形式！

		float *dims[3] = {lon, lat, alt};
		uint32_t codes[3];
		for (size_t i = 0; i < 3; i++)
		{
			uint32_t d_code = abs(dims[i][0]);
			uint32_t m_code = abs(dims[i][1]);
			uint32_t s_code = abs((int)dims[i][2]);
			codes[i] = (d_code << 12) | (m_code << 6) | s_code;
			// for altitude, the code is simply increased from the lower
			if (i < 2)
			{
				bool pos = dims[i][0] >= 0 && dims[i][1] >= 0 && dims[i][2] >= 0;
				auto bit = pos ? 0 : 1;
				codes[i] = SetBit(codes[i], 20, bit);
			}
			codes[i] >>= (20 - level + 1);
		}
		//空间填充曲线编码 (Z序曲线)。它调用了第三方库 libmorton，将处理好的三个维度的整数编码，通过位交错的方式，“编织”成一个最终的64位一维整数。这就是降维的过程。
		return libmorton::morton3D_64_encode(codes[0], codes[1], codes[2]);
	}

	/**
	 * 将一个64位整数编码转换回三维坐标（经度、纬度、高程）。
	 * 它首先调用第三方库 libmorton，将64位整数编码分解为三个20位整数。
	 * 接着，它将这三个20位整数转换为度分秒 (DMS) 格式，并进行位交错处理。
	 * 最后，它将经度、纬度、高程的DMS值转换为十进制度数，并进行规范化处理。
	 */
	void Decode(float out_lon[3], float out_lat[3], float &out_alt_val, uint64_t code, short level)
	{
		uint_fast32_t lon_code;
		uint_fast32_t lat_code;
		uint_fast32_t alt_code;
		libmorton::morton3D_64_decode(code, lon_code, lat_code, alt_code);
		float out_alt[3];
		float *dims_out[3] = {out_lon, out_lat, out_alt};
		uint_fast32_t dims_code[3] = {lon_code, lat_code, alt_code};
		for (size_t i = 0; i < 3; i++)
		{
			// if first bit is 1, then negative for lon & lat
			auto dim_code = dims_code[i];
			dim_code <<= (20 - level + 1);
			bool one_leading = dim_code >= (1 << 20);
			int sign = 1;
			// if not altitude
			if (i < 2)
			{
				dim_code = SetBit(dim_code, 20, 0);
				sign = (one_leading ? -1 : 1);
			}
			auto code_d = (dim_code >> 12);
			auto code_m = (dim_code >> 6) & 0x3f;
			auto code_s = dim_code & 0x3f;
			dims_out[i][0] = *reinterpret_cast<int32_t *>(&code_d) * sign;
			dims_out[i][1] = *reinterpret_cast<int32_t *>(&code_m) * sign;
			dims_out[i][2] = *reinterpret_cast<int32_t *>(&code_s) * sign;
		}
		auto alt_deg = DMSToDegree(out_alt);
		out_alt_val = alt_deg / 512.0 * (50000 + EARTH_R) - EARTH_R;
	}

	uint64_t Encode(float lon, float lat, float alt, short level)
	{
		float lon_dms[3];
		float lat_dms[3];
		DegreeToDMS(lon_dms, lon);
		DegreeToDMS(lat_dms, lat);
		return Encode(lon_dms, lat_dms, alt, level);
	}

	void Decode(float &out_lon, float &out_lat, float &out_alt, uint64_t code, short level)
	{
		float lon_dms[3];
		float lat_dms[3];
		Decode(lon_dms, lat_dms, out_alt, code, level);
		out_lon = DMSToDegree(lon_dms);
		out_lat = DMSToDegree(lat_dms);
	}
}