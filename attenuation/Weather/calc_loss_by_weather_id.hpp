#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <algorithm>
#include <iostream>

namespace loss {

    // -------------------- helpers --------------------
    inline double deg2rad(double deg) { return deg * M_PI / 180.0; }
    inline double clamp_min(double x, double mn) { return (x < mn) ? mn : x; }

    // -------------------- user structs --------------------
    struct Met {
        // clear (gas)
        double p_hPa = 1013.0;
        double T_C = 15.0;
        double rho_w_gm3 = 7.5;

        // rain
        double R_mmph = 10.0;

        // fog/cloud
        double rhoL_gm3 = 0.1;
        double cloudT_K = 273.75;

        // snow
        double S_liqeq_mmph = 1.0;
        std::string snowType = "wet"; // "wet" | "dry"
    };

    struct Options {
        std::string pol = "H";   // "H"/"V"/others
        double elev_deg = 0.0;
        double drySnowScale = 0.2;

        // 可选：如果你要强制 tau_deg，就把 has_tau_deg=true
        bool   has_tau_deg = false;
        double tau_deg = 0.0;
    };

    enum class WeatherType { none, clear, rain, snow, fog };

    struct OutputDetail {
        double clear_gas_dB = 0.0;
        double rain_dB = 0.0;
        double fog_cloud_dB = 0.0;
        double snow_dB = 0.0;
    };

    struct Output {
        int weather_id = 0;
        WeatherType weatherType = WeatherType::none;
        double f_Hz = 0.0;
        double r_km = 0.0;

        double FSPL_dB = 0.0;
        double A_weather_dB = 0.0;
        double totalLoss_linear_amp = 0.0;

        OutputDetail detail;
    };

    // -------------------- default mapping --------------------
    inline WeatherType weatherId2Type(int weather_id) {
        // 默认：0=none, 1=clear, 2=rain, 3=snow, 4=fog
        switch (weather_id) {
        case 1: return WeatherType::clear;
        case 2: return WeatherType::rain;
        case 3: return WeatherType::snow;
        case 4: return WeatherType::fog;
        default: return WeatherType::none;
        }
    }

    inline double getTauDeg(const Options& opts) {
        if (opts.has_tau_deg) return opts.tau_deg;

        std::string pol = opts.pol;
        for (auto& c : pol) c = (char)std::toupper((unsigned char)c);

        if (pol == "H" || pol == "HH") return 0.0;
        if (pol == "V" || pol == "VV") return 90.0;
        return 45.0;
    }

    // ================= P.838：k, alpha (H/V) =================
    inline void p838_k_alpha_H(double f_GHz, double& kH, double& aH) {
        const double x = std::log10(f_GHz);

        const double aj_k[] = { -5.33980, -0.35351, -0.23789, -0.94158 };
        const double bj_k[] = { -0.10008,  1.26970,  0.86036,  0.64552 };
        const double cj_k[] = { 1.13098,  0.45400,  0.15354,  0.16817 };
        const double mk = -0.18961, ck = 0.71147;

        double sumk = 0.0;
        for (int i = 0; i < 4; ++i) {
            const double t = (x - bj_k[i]) / cj_k[i];
            sumk += aj_k[i] * std::exp(-(t * t));
        }
        kH = std::pow(10.0, sumk + mk * x + ck);

        const double aj_a[] = { -0.14318, 0.29591, 0.32177, -5.37610, 16.1721 };
        const double bj_a[] = { 1.82442, 0.77564, 0.63773, -0.96230, -3.29980 };
        const double cj_a[] = { -0.55187, 0.19822, 0.13164,  1.47828,  3.43990 };
        const double ma = 0.67849, ca = -1.95537;

        double suma = 0.0;
        for (int i = 0; i < 5; ++i) {
            const double t = (x - bj_a[i]) / cj_a[i];
            suma += aj_a[i] * std::exp(-(t * t));
        }
        aH = (suma + ma * x + ca);
    }

    inline void p838_k_alpha_V(double f_GHz, double& kV, double& aV) {
        const double x = std::log10(f_GHz);

        const double aj_k[] = { -3.80595, -3.44965, -0.39902, 0.50167 };
        const double bj_k[] = { 0.56934, -0.22911,  0.73042, 1.07319 };
        const double cj_k[] = { 0.81061,  0.51059,  0.11899, 0.27195 };
        const double mk = -0.16398, ck = 0.63297;

        double sumk = 0.0;
        for (int i = 0; i < 4; ++i) {
            const double t = (x - bj_k[i]) / cj_k[i];
            sumk += aj_k[i] * std::exp(-(t * t));
        }
        kV = std::pow(10.0, sumk + mk * x + ck);

        const double aj_a[] = { -0.07771, 0.56727, -0.20238, -48.2991, 48.5833 };
        const double bj_a[] = { 2.33840, 0.95545,  1.14520,  0.791669, 0.791459 };
        const double cj_a[] = { -0.76284, 0.54039,  0.26809,  0.116226, 0.116479 };
        const double ma = -0.053739, ca = 0.83433;

        double suma = 0.0;
        for (int i = 0; i < 5; ++i) {
            const double t = (x - bj_a[i]) / cj_a[i];
            suma += aj_a[i] * std::exp(-(t * t));
        }
        aV = (suma + ma * x + ca);
    }

    // ================= P.838：雨衰 gammaR(dB/km) =================
    inline double p838_rain_specific_dBpkm(double f_GHz, double R_mmph, double elev_deg, double tau_deg) {
        const double f_eff = (std::max)(f_GHz, 1.0);

        double kH, aH, kV, aV;
        p838_k_alpha_H(f_eff, kH, aH);
        p838_k_alpha_V(f_eff, kV, aV);

        const double theta = deg2rad(elev_deg);
        const double tau = deg2rad(tau_deg);

        const double ctheta2 = std::cos(theta) * std::cos(theta);
        const double cos2tau = std::cos(2.0 * tau);

        const double k = (kH + kV + (kH - kV) * ctheta2 * cos2tau) / 2.0;

        // 防止除 0
        const double k_safe = (std::abs(k) < 1e-18) ? 1e-18 : k;

        const double alpha =
            (kH * aH + kV * aV + (kH * aH - kV * aV) * ctheta2 * cos2tau) / (2.0 * k_safe);

        return k_safe * std::pow(R_mmph, alpha);
    }

    // ================= P.840：Kl (dB/km per g/m^3) =================
    inline double p840_Kl_dBpkm_per_gm3(double f_GHz, double T_K) {
        const double f = f_GHz;
        const double T = T_K;

        const double theta = 300.0 / T;
        const double eps0 = 77.66 + 103.3 * (theta - 1.0);
        const double eps1 = 0.0671 * eps0;
        const double eps2 = 3.52;
        const double fp = 20.20 - 146.0 * (theta - 1.0) + 316.0 * (theta - 1.0) * (theta - 1.0);
        const double fs = 39.8 * fp;

        const double fp2 = fp * fp;
        const double fs2 = fs * fs;

        const double termp = 1.0 + (f * f) / fp2;
        const double terms = 1.0 + (f * f) / fs2;

        const double epspp =
            f * (eps0 - eps1) / (fp * termp) +
            f * (eps1 - eps2) / (fs * terms);

        const double epsp =
            (eps0 - eps1) / termp +
            (eps1 - eps2) / terms +
            eps2;

        return 0.819 * f * epspp / ((epsp + 2.0) * (epsp + 2.0) + epspp * epspp);
    }

    // ================= P.676：气体吸收（近似） gamma_o, gamma_w (dB/km) =================
    inline void p676_gas_specific_dBpkm(double f_GHz,
        double p_hPa,
        double T_C,
        double rho_gm3,
        double& gamma_o,
        double& gamma_w)
    {
        gamma_o = 0.0;
        gamma_w = 0.0;

        if (!(f_GHz >= 1.0 && f_GHz <= 350.0)) return;

        const double p = p_hPa;
        const double t = T_C;
        const double rho = rho_gm3;

        const double rp = p / 1013.0;
        const double rt = 288.0 / (273.0 + t);

        auto phi = [&](double a, double b, double c, double d) {
            return std::pow(rp, a) * std::pow(rt, b) * std::exp(c * (1.0 - rp) + d * (1.0 - rt));
            };

        const double xi1 = phi(0.0717, -1.8132, 0.0156, -1.6515);
        const double xi2 = phi(0.5146, -4.6368, -0.1921, -5.7416);
        const double xi3 = phi(0.3414, -6.5851, 0.2130, -8.5854);
        const double xi4 = phi(-0.0112, 0.0092, -0.1033, -0.0009);
        const double xi5 = phi(0.2705, -2.7192, -0.3016, -4.1033);
        const double xi6 = phi(0.2445, -5.9191, 0.0422, -8.0719);
        const double xi7 = phi(-0.1833, 6.5589, -0.2402, 6.131);

        const double g54 = 2.192 * phi(1.8286, -1.9487, 0.4051, -2.8509);
        const double g58 = 12.59 * phi(1.0045, 3.5610, 0.1588, 1.2834);
        const double g60 = 15.0 * phi(0.9003, 4.1335, 0.0427, 1.6088);
        const double g62 = 14.28 * phi(0.9886, 3.4176, 0.1827, 1.3429);
        const double g64 = 6.819 * phi(1.4320, 0.6258, 0.3177, -0.5914);
        const double g66 = 1.908 * phi(2.0717, -4.1404, 0.4910, -4.8718);
        const double delta = -0.00306 * phi(3.211, -14.94, 1.583, -16.37);

        // ---- gamma_o ----
        const double f = f_GHz;
        double go = 0.0;

        if (f <= 54.0) {
            const double denom1 = (f * f + 0.34 * rp * rp * std::pow(rt, 1.6));
            const double term1 = 7.2 * std::pow(rt, 2.8) / denom1;

            const double denom2 = std::pow(54.0 - f, 1.16 * xi1) + 0.83 * xi2;
            const double term2 = 0.62 * xi3 / denom2;

            go = (term1 + term2) * (f * f) * (rp * rp) * 1e-3;
        }
        else if (f > 54.0 && f <= 60.0) {
            go = std::exp((std::log(g54) / 24.0) * (f - 58.0) * (f - 60.0)
                - (std::log(g58) / 8.0) * (f - 54.0) * (f - 60.0)
                + (std::log(g60) / 12.0) * (f - 54.0) * (f - 58.0));
        }
        else if (f > 60.0 && f <= 62.0) {
            go = g60 + (g62 - g60) * (f - 60.0) / 2.0;
        }
        else if (f > 62.0 && f <= 66.0) {
            go = std::exp((std::log(g62) / 8.0) * (f - 64.0) * (f - 66.0)
                - (std::log(g64) / 4.0) * (f - 62.0) * (f - 66.0)
                + (std::log(g66) / 8.0) * (f - 62.0) * (f - 64.0));
        }
        else if (f > 66.0 && f <= 120.0) {
            const double denomA = ((f - 118.75) * (f - 118.75) + 2.91 * rp * rp * std::pow(rt, 1.6));
            const double termA = 0.283 * std::pow(rt, 3.8) / denomA;

            const double denomB = (std::pow(f - 66.0, 1.4346 * xi4) + 1.15 * xi5);
            const double termB = 0.502 * xi6 * (1.0 - 0.0163 * xi7 * (f - 66.0)) / denomB;

            const double termC = 3.02e-4 * std::pow(rt, 3.5);

            go = (termC + termA + termB) * (f * f) * (rp * rp) * 1e-3;
        }
        else { // (120, 350]
            const double denomA = (1.0 + 1.9e-5 * std::pow(f, 1.5));
            const double term1 = 3.02e-4 / denomA;

            const double denom2 = ((f - 118.75) * (f - 118.75) + 2.91 * rp * rp * std::pow(rt, 1.6));
            const double term2 = 0.283 * std::pow(rt, 0.3) / denom2;

            go = (term1 + term2) * (f * f) * (rp * rp) * std::pow(rt, 3.5) * 1e-3 + delta;
        }

        gamma_o = go;

        // ---- gamma_w ----
        const double eta1 = 0.955 * rp * std::pow(rt, 0.68) + 0.006 * rho;
        // eta2 在 MATLAB 里算了但没用，这里也不需要

        auto gfun = [&](double ff, double fi) {
            const double num = (ff - fi);
            const double den = (ff + fi);
            const double q = (den == 0.0) ? 0.0 : (num / den);
            return 1.0 + q * q;
            };

        const double fw = f;
        const double part =
            (3.98 * eta1 * std::exp(2.23 * (1.0 - rt)) / ((fw - 22.235) * (fw - 22.235) + 9.42 * eta1 * eta1)) * gfun(fw, 22.0) +
            (11.96 * eta1 * std::exp(0.70 * (1.0 - rt)) / ((fw - 183.31) * (fw - 183.31) + 11.14 * eta1 * eta1)) +
            (0.081 * eta1 * std::exp(6.44 * (1.0 - rt)) / ((fw - 321.226) * (fw - 321.226) + 6.29 * eta1 * eta1)) +
            (3.66 * eta1 * std::exp(1.60 * (1.0 - rt)) / ((fw - 325.153) * (fw - 325.153) + 9.22 * eta1 * eta1));

        gamma_w = part * (fw * fw) * std::pow(rt, 2.5) * rho * p * 1e-4;
    }

    // ================= 主函数：calcLossByWeatherId =================
    inline Output calcLossByWeatherId(int weather_id = 1,
        double f_Hz = 1000,
        double r_km = 1000,
        const Met* met_in = nullptr,
        const Options* opts_in = nullptr)
    {
        Met met = met_in ? *met_in : Met{};
        Options opts = opts_in ? *opts_in : Options{};

        Output out;
        out.weather_id = weather_id;
        out.weatherType = weatherId2Type(weather_id);
        out.f_Hz = f_Hz;
        out.r_km = r_km;

        const double f_GHz = f_Hz / 1e9;

        // 1) FSPL (与 MATLAB 相同：92.45 + 20log10(f_GHz) + 20log10(r_km))
        const double fG = clamp_min(f_GHz, 1e-12);
        const double rk = clamp_min(r_km, 1e-9);
        out.FSPL_dB = 92.45 + 20.0 * std::log10(fG) + 20.0 * std::log10(rk);

        // 2) weather attenuation: gamma(dB/km) * r_km
        const double tau_deg = getTauDeg(opts);

        double A_clear = 0.0, A_rain = 0.0, A_fog = 0.0, A_snow = 0.0;

        switch (out.weatherType) {
        case WeatherType::none:
            break;

        case WeatherType::clear: {
            double g_o = 0.0, g_w = 0.0;
            p676_gas_specific_dBpkm(f_GHz, met.p_hPa, met.T_C, met.rho_w_gm3, g_o, g_w);
            A_clear = (g_o + g_w) * r_km;
        } break;

        case WeatherType::rain: {
            const double gR = p838_rain_specific_dBpkm(f_GHz, met.R_mmph, opts.elev_deg, tau_deg);
            A_rain = gR * r_km;
        } break;

        case WeatherType::fog: {
            const double Kl = p840_Kl_dBpkm_per_gm3(f_GHz, met.cloudT_K);
            A_fog = (Kl * met.rhoL_gm3) * r_km;
        } break;

        case WeatherType::snow: {
            const std::string st = [&] {
                std::string s = met.snowType;
                for (auto& c : s) c = (char)std::tolower((unsigned char)c);
                return s;
                }();

            double R_eq = met.S_liqeq_mmph;
            if (st == "dry") R_eq = opts.drySnowScale * met.S_liqeq_mmph;

            const double gS = p838_rain_specific_dBpkm(f_GHz, R_eq, opts.elev_deg, tau_deg);
            A_snow = gS * r_km;
        } break;
        }

        out.detail.clear_gas_dB = A_clear;
        out.detail.rain_dB = A_rain;
        out.detail.fog_cloud_dB = A_fog;
        out.detail.snow_dB = A_snow;

        out.A_weather_dB = A_clear + A_rain + A_fog + A_snow;
        // 计算雷达传播模型的单程幅度线性衰减倍数，而非对数dB！
        out.totalLoss_linear_amp = std::pow(10.0, -(out.FSPL_dB + out.A_weather_dB) / 20.0);
        return out;
    }

    // 你之前 C++ 调用习惯：直接要 totalLoss
    inline double calcTotalLoss_linear_amp(int weather_id,
        double f_Hz,
        double r_km,
        const Met* met = nullptr,
        const Options* opt = nullptr)
    {
        return calcLossByWeatherId(weather_id, f_Hz, r_km, met, opt).totalLoss_linear_amp;
    }

} // namespace loss
