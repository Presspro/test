#include "weather/calc_loss_by_weather_id.hpp"

int main() {
	int weather_id = 1;					//case 1:clear;case 2:rain;case 3:snow;case 4:fog;default::none;
	double f_Hz = 1000;						//载频
	double r_km = 1000;						//距离
	//met参数
	loss::Met met;

   met.p_hPa = 1013.25;       // 大气压强 (百帕, hPa) - 标准海平面大气压
   met.T_C = 15.0;            // 环境温度 (摄氏度, ℃)
   met.rho_w_gm3 = 7.5;       // 绝对水汽密度 (克/立方米, g/m³) - 影响大气的氧气和水汽吸收衰减

   met.R_mmph = 10.0;         // 降雨率 (毫米/小时, mm/h) - 10mm/h 属于中雨级别

   met.rhoL_gm3 = 0.1;        // 云雾的液态水含量密度 (克/立方米, g/m³) - 典型雾天取值 0.05~0.5
   met.cloudT_K = 273.15;     // 云层物理温度 (开尔文, K) - 273.15K 即 0℃

   met.S_liqeq_mmph = 1.0;    // 降雪的等效液态水降水率 (毫米/小时, mm/h) - 将雪融化成水后的降水率
   met.snowType = "wet";      // 雪的物理状态 - "wet" 代表湿雪（衰减较大），"dry" 代表干雪

	//opt参数
	loss::Options opt;

   // [常规极化与几何参数]
   opt.pol = "H";             // 雷达发射/接收极化方式 - "H" 为水平极化，"V" 为垂直极化
   opt.elev_deg = 0.0;        // 雷达波束的仰角/俯仰角 (度, degree) - 0度代表水平扫视
   opt.drySnowScale = 0.2;    // 干雪衰减缩放系数 - 经验常数，干雪对电磁波的衰减通常远小于湿雪

   // [高级极化倾角参数]
   opt.has_tau_deg = false;   // 是否强制覆盖计算中的极化倾角 (tau) - false 则根据 pol 和 elev_deg 自动推导
   opt.tau_deg = 0.0;         // 强制指定的极化倾角 (度, degree) - 仅当 has_tau_deg 为 true 时生效

 	loss::Output output;
	output = loss::calcLossByWeatherId(weather_id, f_Hz, r_km, &met, &opt);

   std::string weather_name = "未知";
   switch (output.weatherType) {
   case loss::WeatherType::none:  weather_name = "无损耗/默认"; break;
   case loss::WeatherType::clear: weather_name = "晴天 (大吸收)"; break;
   case loss::WeatherType::rain:  weather_name = "降雨"; break;
   case loss::WeatherType::fog:   weather_name = "云雾"; break;
   case loss::WeatherType::snow:  weather_name = "降雪"; break;
   }
   std::cout << "当前的频率 (Hz): " << output.r_km << std::endl;
   std::cout << "当前距离 (dB): " << output.f_Hz << std::endl;
   std::cout << "天气类型: " << weather_name << std::endl;
   std::cout << "单程天气衰减 (dB): " << output.A_weather_dB << std::endl;
   std::cout << "自由空间损耗 (dB): " << output.FSPL_dB << std::endl;
   std::cout << "总线性损耗倍数 (单程幅度): " << output.totalLoss_linear_amp << std::endl;

   return 0;
}