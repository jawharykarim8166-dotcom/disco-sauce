#pragma once
#include "Dsp.h"
#include "Presets.h"

namespace ds
{
// Real-time pitch correction.
//  - Pitch detection: YIN (de Cheveigne & Kawahara, 2002) on a decimated mono copy of the input.
//  - Correction: nearest allowed note in the chosen key/scale, with hysteresis and a retune-speed smoother.
//  - Shifting: two-tap crossfading delay-line shifter whose window is kept a multiple of the detected
//    pitch period (pitch-synchronous), which avoids the comb/flutter of a plain delay-line shifter.
class AutoTune
{
public:
    struct Settings
    {
        float amount = 1.0f;     // 0 = off, 1 = full correction
        float retuneMs = 20.0f;  // 0 = instant (hard tune)
        int key = 0;             // 0..11 (C..B)
        int scale = 0;           // index into scale list
    };

    void prepare(double sr)
    {
        sampleRate = sr;
        bufSize = nextPow2((int) (sr * 0.16));
        bufMask = bufSize - 1;
        buf[0].assign((size_t) bufSize, 0.f);
        buf[1].assign((size_t) bufSize, 0.f);
        wi = 0;
        dnom = (int) std::round(0.013 * sr);
        W = 2.0 * dnom;
        Wtarget = W;
        phase = 0.0;

        dec = std::max(1, (int) std::round(sr / 12000.0));
        anaRate = sr / dec;
        anaW = (int) std::round(anaRate * 0.024);
        tauMax = (int) std::round(anaRate / 70.0);
        tauMin = std::max(2, (int) std::round(anaRate / 1000.0));
        anaHop = std::max(16, (int) (anaRate * 0.008));
        ring.assign(2048, 0.f);
        ringW = 0;
        tmp.assign((size_t) (anaW + tauMax + 4), 0.f);
        diff.assign((size_t) (tauMax + 4), 0.0);
        cmnd.assign((size_t) (tauMax + 4), 0.0);
        reset();
    }

    void reset()
    {
        std::fill(buf[0].begin(), buf[0].end(), 0.f);
        std::fill(buf[1].begin(), buf[1].end(), 0.f);
        std::fill(ring.begin(), ring.end(), 0.f);
        accum = 0.f; accCount = 0; hopCount = 0;
        corr = 0.f; desired = 0.f; voicedG = 0.f; engage = 0.f;
        voicedHold = 0; haveDet = false; haveTarget = false;
        lastTarget = 0; detSm = 60.f; lastF = 0.f;
    }

    int getLatency() const { return dnom; }
    bool isVoiced() const { return voicedHold > 0; }
    float getDetectedMidi() const { return detSm; }

    // L/R are processed in place; dryL/dryR receive the input delayed by the reported latency.
    void process(float* L, float* R, float* dryL, float* dryR, int n, const Settings& s)
    {
        const bool active = s.amount > 0.001f;
        const double engC = 1.0 - std::exp(-1.0 / (0.02 * sampleRate));
        const double voicC = 1.0 - std::exp(-1.0 / (0.01 * sampleRate));
        const double corrC = (s.retuneMs <= 0.5f) ? 1.0
                           : 1.0 - std::exp(-1.0 / (0.001 * (double) s.retuneMs * sampleRate));
        const double releaseC = 1.0 - std::exp(-1.0 / (0.03 * sampleRate));
        curSettings = s;

        for (int i = 0; i < n; ++i)
        {
            const float inL = L[i], inR = R[i];
            buf[0][(size_t) wi] = inL;
            buf[1][(size_t) wi] = inR;

            const int di = (wi - dnom) & bufMask;
            const float dL = buf[0][(size_t) di];
            const float dR = buf[1][(size_t) di];
            dryL[i] = dL;
            dryR[i] = dR;

            // analysis path (decimate by averaging)
            accum += 0.5f * (inL + inR);
            if (++accCount >= dec)
            {
                ring[(size_t) (ringW & 2047)] = accum / (float) dec;
                ++ringW;
                accum = 0.f; accCount = 0;
                if (++hopCount >= anaHop) { hopCount = 0; analyse(); }
            }

            engage += (float) (((active ? 1.0f : 0.0f) - engage) * engC);
            voicedG += (float) (((voicedHold > 0 ? 1.0f : 0.0f) - voicedG) * voicC);

            if (voicedHold > 0) corr += (float) ((desired - corr) * corrC);
            else corr += (float) ((0.0f - corr) * std::max(corrC, releaseC));

            if (engage < 0.0005f && !active)
            {
                L[i] = dL; R[i] = dR;
                W += std::max(-0.01, std::min(0.01, Wtarget - W));
                wi = (wi + 1) & bufMask;
                continue;
            }

            const double r = std::exp2((double) corr / 12.0);
            phase += (1.0 - r) / W;
            phase -= std::floor(phase);
            double p2 = phase + 0.5;
            p2 -= std::floor(p2);
            const float w1 = (float) (0.5 - 0.5 * std::cos(2.0 * kPi * phase));
            const float w2 = (float) (0.5 - 0.5 * std::cos(2.0 * kPi * p2));
            const double d1 = phase * W, d2 = p2 * W;
            const double dm = 0.5 * W;

            float out[2];
            for (int c = 0; c < 2; ++c)
            {
                const float sh = w1 * readFrac(buf[c], d1) + w2 * readFrac(buf[c], d2);
                const float tap = readFrac(buf[c], dm);
                const float tuned = voicedG * sh + (1.0f - voicedG) * tap;
                const float dry = (c == 0) ? dL : dR;
                out[c] = engage * tuned + (1.0f - engage) * dry;
            }
            L[i] = out[0];
            R[i] = out[1];

            W += std::max(-0.01, std::min(0.01, Wtarget - W));
            wi = (wi + 1) & bufMask;
        }
    }

private:
    inline float readFrac(const std::vector<float>& b, double d) const
    {
        const double rp = (double) wi - d;
        const int i1 = (int) std::floor(rp);
        const float f = (float) (rp - (double) i1);
        const float y0 = b[(size_t) ((i1 - 1) & bufMask)];
        const float y1 = b[(size_t) (i1 & bufMask)];
        const float y2 = b[(size_t) ((i1 + 1) & bufMask)];
        const float y3 = b[(size_t) ((i1 + 2) & bufMask)];
        const float c0 = y1;
        const float c1 = 0.5f * (y2 - y0);
        const float c2 = y0 - 2.5f * y1 + 2.0f * y2 - 0.5f * y3;
        const float c3 = 0.5f * (y3 - y0) + 1.5f * (y1 - y2);
        return ((c3 * f + c2) * f + c1) * f + c0;
    }

    bool isAllowed(int note) const
    {
        const int pc = (((note - curSettings.key) % 12) + 12) % 12;
        return ((scaleMask(curSettings.scale) >> pc) & 1) != 0;
    }

    int findTarget(float midi)
    {
        if (haveTarget && isAllowed(lastTarget) && std::abs(midi - (float) lastTarget) < 0.65f)
            return lastTarget;
        const int c = (int) std::lround(midi);
        int best = c;
        float bestD = 1.0e9f;
        for (int n = c - 7; n <= c + 7; ++n)
        {
            if (!isAllowed(n)) continue;
            const float d = std::abs(midi - (float) n);
            if (d < bestD) { bestD = d; best = n; }
        }
        haveTarget = true;
        lastTarget = best;
        return best;
    }

    // YIN pitch detection on the decimated ring buffer.
    bool detect(float& freq)
    {
        const int N = anaW + tauMax;
        for (int j = 0; j < N; ++j)
            tmp[(size_t) j] = ring[(size_t) ((ringW - N + j) & 2047)];

        double e = 0.0;
        for (int j = tauMax; j < N; ++j) e += (double) tmp[(size_t) j] * tmp[(size_t) j];
        if (std::sqrt(e / anaW) < 0.003) return false;

        for (int tau = 1; tau <= tauMax; ++tau)
        {
            double d = 0.0;
            for (int j = 0; j < anaW; ++j)
            {
                const double v = (double) tmp[(size_t) j] - (double) tmp[(size_t) (j + tau)];
                d += v * v;
            }
            diff[(size_t) tau] = d;
        }
        cmnd[0] = 1.0;
        double run = 0.0;
        for (int tau = 1; tau <= tauMax; ++tau)
        {
            run += diff[(size_t) tau];
            cmnd[(size_t) tau] = run > 1.0e-12 ? diff[(size_t) tau] * tau / run : 1.0;
        }

        int est = -1;
        for (int tau = tauMin; tau < tauMax; ++tau)
        {
            if (cmnd[(size_t) tau] < 0.15)
            {
                while (tau + 1 < tauMax && cmnd[(size_t) (tau + 1)] < cmnd[(size_t) tau]) ++tau;
                est = tau;
                break;
            }
        }
        if (est < 2) return false;

        const double s0 = cmnd[(size_t) (est - 1)], s1 = cmnd[(size_t) est], s2 = cmnd[(size_t) (est + 1)];
        const double den = s0 - 2.0 * s1 + s2;
        const double shift = std::abs(den) > 1.0e-9 ? 0.5 * (s0 - s2) / den : 0.0;
        const double better = (double) est + std::max(-1.0, std::min(1.0, shift));
        freq = (float) (anaRate / better);
        return freq > 60.f && freq < 1200.f;
    }

    void analyse()
    {
        float f = 0.f;
        if (detect(f))
        {
            const float midi = 69.0f + 12.0f * std::log2(f / 440.0f);
            if (!haveDet || std::abs(midi - detSm) > 1.0f) detSm = midi;
            else detSm += 0.5f * (midi - detSm);
            haveDet = true;

            const int target = findTarget(detSm);
            desired = ((float) target - detSm) * curSettings.amount;
            voicedHold = 5;

            // pitch-synchronous window: half-window = whole number of periods, ~13 ms
            const double P = sampleRate / (double) f;
            const int K = std::max(1, (int) std::lround((double) dnom / P));
            Wtarget = std::min(2.0 * K * P, (double) (bufSize / 2 - 8));
            lastF = f;
        }
        else
        {
            if (voicedHold > 0) --voicedHold;
            if (voicedHold == 0) { haveDet = false; haveTarget = false; }
        }
    }

    double sampleRate = 48000.0, anaRate = 12000.0;
    int bufSize = 1, bufMask = 0, wi = 0, dnom = 624;
    std::vector<float> buf[2];
    double W = 1248.0, Wtarget = 1248.0, phase = 0.0;

    int dec = 4, anaW = 288, tauMax = 171, tauMin = 12, anaHop = 96;
    std::vector<float> ring, tmp;
    std::vector<double> diff, cmnd;
    int ringW = 0, accCount = 0, hopCount = 0;
    float accum = 0.f;

    float corr = 0.f, desired = 0.f, voicedG = 0.f, engage = 0.f, detSm = 60.f, lastF = 0.f;
    int voicedHold = 0, lastTarget = 0;
    bool haveDet = false, haveTarget = false;
    Settings curSettings;
};

} // namespace ds
